######################## imports
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)
library(glmnet)


#' LASSO with Differential Shrinkage and Stability Selection for GRN inference
#'
#' @param counts expression matrix with gene IDs as rownames and conditions in columns
#' @param genes list of genes used as inputs for GRN inference
#' @param tfs list of TFs used as predictors for GRN inference
#' @param alpha integration strength, value should be between 0 and 1.
#' @param scale weather or not to scale expression data to z-scores.
#' @param pwm_occurrence matrix of prior data Pi containing TFBS scores between
#' TFs and genes
#' @param int_pwm_noise value of noise added to differential shrinkage 
#' during stability selection
#' @param N Number of stability selection iterations
#' @param score scoring of the interactions "pval" or "freq"
#' @param robustness threshold of frequency selection required for a 
#' TF to be included in unpenalized regressions
#' @param nCores Number of cores for multithreading. (Not supported on Windows)
#'
#' @return list of regulatory interactions weighted by p-value
LASSO.D3S_inference <- function(counts, genes, tfs, alpha=0.25, 
                                   pwm_occurrence, int_pwm_noise = 0.1,
                                   N = 100, 
                                   score = "pval", robustness = 0.1,
                                   nCores = ifelse(is.na(detectCores()),1,
                                                   max(detectCores() - 1, 1))){
  
  
  # to avoid convergence issues : "inner loop 3; cannot correct step size"
  maxit = 1e+05
  if(alpha==1){
    alpha=1-1e-8
    maxit = 1e+07
  }
  
  
  counts <- round(counts, 0)
  x <- t(counts[tfs,])
  
  # pwm scores to bias variable selecttion toward pairs supported by a TFBS
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  # parallel computing of the lasso
  registerDoParallel(cores = nCores)
  message(paste("\nUsing LASSO-D3S with", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(result.reg <-
                                   doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = rbind, .inorder = T),
                                                    {
                                                      # getting rid of the TF variable on which the regression is made
                                                      # if needed
                                                      target_tfs <- setdiff(tfs, target)
                                                      x_target <- x[, target_tfs]
                                                      y <- t(counts[target, ])
                                                      
                                                      # weights for differential shrinkage
                                                      penalty_factor <- 1 - pwm_imputed[target, target_tfs] * alpha
                                                      
                                                      importances <- setNames(rep(0, length(tfs)), tfs)
                                                      
                                                      ######### Stability Selection
                                                      n_actual = 0
                                                      for(n in 1:N) {
                                                        # bootstrapping observations
                                                        sampled <- sample(1:nrow(x), replace = T, size = nrow(x))
                                                        oob <- setdiff(1:nrow(x), sampled)
                                                        
                                                        if(length(oob)>1){ #ensures a test set of sufficient size
                                                          
                                                          # perturbating differential shrinkage
                                                          noisy_penalty_factor <- pmax(pmin(penalty_factor + 
                                                                                              runif(n=length(penalty_factor),
                                                                                                    min = -int_pwm_noise*alpha, 
                                                                                                    max = int_pwm_noise*alpha), 1),0)
                                                          
                                                          # models are try-catched because in rare cases glmnet
                                                          # crashes for convergence issues
                                                          tryCatch(
                                                            error = function(cnd) "cv.glmnet internal error due to convergence issues",
                                                            
                                                            # model training on sampled observations
                                                            {mymodels_pen = glmnet(
                                                              x_target[sampled,],
                                                              y[sampled],
                                                              maxit=maxit,
                                                              family = "poisson",
                                                              penalty.factor = noisy_penalty_factor)
                                                            
                                                            # model predictions on OOB observations
                                                            y_hat <- predict.glmnet(mymodels_pen, newx = x_target[oob,],
                                                                                    type = "response")
                                                            
                                                            # getting lambda value granting minimal MSE 
                                                            # on OOB observations
                                                            mses_oob <- c()
                                                            for(s in colnames(y_hat)){
                                                              mses_oob <- c(mses_oob, mean((y[oob] - y_hat[,s])^2))
                                                            }
                                                            iLambdaMin = which(mses_oob == min(mses_oob))
                                                            
                                                            # feature selection for optimal lambda
                                                            selected_tfs <- names(which(mymodels_pen$beta[, iLambdaMin] != 0))
                                                            
                                                            # selection frequency
                                                            importances[selected_tfs] <- importances[selected_tfs]+1
                                                            n_actual = n_actual+1
                                                            }
                                                          )
                                                        }
                                                      }
                                                      
                                                      
                                                      
                                                      robust_tfs  = NULL
                                                      if(score=="freq")
                                                        to_return <- data.frame(from =names(importances),
                                                                                to = target,
                                                                                importance = importances/n_actual)
                                                      ###### Unpenalized regression on robust TFs
                                                      if(score=="pval"){
                                                        to_return <- data.frame(from =names(importances),
                                                                                to = target,
                                                                                importance = 1)
                                                        tryCatch(
                                                          error = function(cnd) "a problem occurred in the unpenalized regressions for this gene.",
                                                          
                                                          # model training on sampled observations
                                                          {robust_tfs <- names(importances[importances>=n_actual*robustness])
                                                          # non penalized models with only robustly selected TFs
                                                          if(length(robust_tfs)>0){
                                                            lm_target <- glm(
                                                              formula = paste(paste0("`", target, "`"), '~', 
                                                                              paste(paste0("`", robust_tfs, "`"), 
                                                                                    collapse = '+')), 
                                                              data = data.frame(t(counts), check.names = F), 
                                                              family = "poisson")
                                                            
                                                            pvals <- c(setNames(rep(1, length(setdiff(tfs, robust_tfs))), 
                                                                                nm = setdiff(tfs, robust_tfs)), 
                                                                       setNames(summary(lm_target)$coefficients[
                                                                         2:nrow(summary(lm_target)$coefficients),"Pr(>|z|)"], 
                                                                         robust_tfs))
                                                            #TF-target interactions are weighted by the TF pvalue in the glm
                                                            to_return <- data.frame(from =names(pvals),
                                                                                    to = target,
                                                                                    importance = pvals)
                                                          }
                                                          })
                                                      }
                                                      if(score=="freq" | score=="pval" & length(robust_tfs)>0)
                                                        to_return
                                                      
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  edges <- result.reg
  toc()
  return(edges)
}


#' Threshold LASSO-D3S GRN to a desired density
#'
#' @param mat result of LASSO.D3S_inference function
#' @param density desired network density
#' @param pwm_occurrence matrix of prior data Pi containing TFBS scores between
#' TFs and genes
#' @param genes list of genes used as inputs for GRN inference
#' @param tfs list of TFs used as predictors for GRN inference
#'
#' @return dataframe of oriented edges, and their prior value in pwm_occurrence
LASSO.D3S_network <- function(mat, density, pwm_occurrence, genes, 
                              tfs, freq=FALSE, with_score = FALSE){
  edges_ <- mat[order(mat$importance, decreasing = freq),] 
  
  nEdges = round(density * (length(genes) - 1) * length(tfs), 0)
  edges_ <- edges_[1:nEdges,]
  
  if(!freq)
    edges_ <- edges_[edges_$importance < 1,]
  
  if(nrow(edges_) < nEdges)
    warning(paste("The desired density can not be reached : too few interactions 
            were pre-selected by the robustness threshold. The actual density is",
                  round(nrow(edges_)/((length(genes) - 1) * length(tfs)), 6), 
                  ", which amounts to", nrow(edges_), "edges."))
  
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  edges_$pwm <- pwm_imputed[cbind(edges_$to, edges_$from)]
  if (!with_score)
    return(edges_[,c("from", "to", "pwm")])
  else
    return(edges_[,c("from", "to", "pwm", "importance")])
    
}


