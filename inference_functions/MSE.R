######################## imports for parallel version
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)
library(igraph)
library(boot)
library(randomForest)
library(stringr)


#' weightedRF GRN inference returns MSE
#'
#' @param counts expression matrix with gene IDs as rownames and conditions in columns
#' @param genes list of genes used as inputs for GRN inference
#' @param tfs list of TFs used as predictors for GRN inference
#' @param alpha integration strength, value should be between 0 and 1.
#' @param scale weather or not to scale expression data to z-scores.
#' @param pwm_occurrence matrix of prior data Pi containing TFBS scores between
#' TFs and genes
#' @param nTrees Number of trees in Random Forests
#' @param importance Importance metric in Random Forests.
#' MDA ("%IncMSE") is the default value. MDI can be used via "IncNodePurity".
#' @param nCores Number of cores for multithreading. (Not supported on Windows)
#'
#' @return The weighted list of regulatory interactions between genes and TFs
weightedRF_inference_MSE <- function(counts, genes, tfs, alpha=0.25, scale = FALSE,
                          pwm_occurrence, nTrees=500, importance="%IncMSE",
                          tf_expression_permutation = FALSE,
                          seed =sample(1:10000, 1),
                          nCores = ifelse(is.na(detectCores()),1,
                                          max(detectCores() - 1, 1))){
  
  # counts must be normalized so that genes have comparable node purities
  if(importance=="IncNodePurity") scale = TRUE
  
  # z-score if scaling is required
  if(scale){
    counts <- (counts - rowMeans(counts))/genefilter::rowSds(counts)
  }
  
  x <- t(counts[tfs,])
  
  # attributing 0.5 for prior value for PWM with unknown PWM
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  # the regressions for each genes are done in parallel
  registerDoParallel(cores = nCores)
  message(paste("\nbRF is running using", foreach::getDoParWorkers(), "cores. alpha =", alpha))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(result.reg <-
                                   doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = c,
                                                                     .final = function(x) {names(x) <- genes; x}, 
                                                                     .inorder = TRUE),
                                                    {
                                                      target_tfs <- setdiff(tfs, target)
                                                      x_target <- x[, target_tfs]
                                                      if(tf_expression_permutation){
                                                        # randomises the expression rows of TFs but not their ID
                                                        set.seed(sample(1:10000, 1))
                                                        x_target <- x_target[,sample(target_tfs, replace = F, 
                                                                                     size = length(target_tfs))]
                                                        colnames(x_target) <- target_tfs
                                                      }
                                                        
                                                      p = length(target_tfs)
                                                      y <- as.numeric(t(counts[target, ]))
                                                      
                                                      # stronger version
                                                      weights <- ifelse(pwm_imputed[target, target_tfs] == 1, sqrt(1-(alpha-1)^2)+1,
                                                                        ifelse(pwm_imputed[target, target_tfs] == 0.5, 1-alpha,
                                                                               -sqrt(1-(alpha-1)^2)+1))
                                                      
                                                      
                                                      
                                                      if(sum(weights)>0){ 
                                                        weights <- weights/sum(weights)
                                                        
                                                        set.seed(seed)
                                                        rf_weighted <-irafnet_onetarget(x_target,y=y,importance=TRUE,
                                                                                        mtry=round(sqrt(p)),
                                                                                        ntree=nTrees,
                                                                                        sw=weights)

                                                        # mean(rf_weighted$mse)/var(y)
                                                        mean((rf_weighted$predicted - y)^2)/var(y)
                                                      }
                                                      else{
                                                        # case when all weights are equal to 0 (no pwm found and alpha = 1)
                                                        # the mse is the error from predicting the gene mean
                                                        sampled <- sample(1:ncol(counts), replace = T, size = ncol(counts))
                                                        oob <- setdiff(1:ncol(counts), sampled)
                                                        mean((counts[target,oob] - mean(counts[target, sampled]))^2)/var(y)
                                                      }
                                                      
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  mat <- result.reg
  toc()
  return(mat)
}




#' weightedLASSO GRN inference returns MSE
#'
#' @param counts Expression matrix (genes in rownames, conditions in columns)
#' @param genes Vector of genes (in the rownames of counts) to be used in GRN inference as target genes
#' @param tfs vector of genes (in the rownames of counts) that are transcriptional regulators
#' to be used a predictors in the regressions for GRN inference
#' @param alpha The strength of data integration.
#' Numeric value (e.g 0, 0.5, 1) 
#' @param pwm_occurrence Prior matrix Pi, giving PWM presence scores for TFs in rows
#' and genes in columns. Can contain NAs for TFs that do not have a PWM available.
#' @param int_pwm_noise Random perturbation applied to PWM priors in pwm_occurrence.
#' Defalut is none, experimental.
#' @param N Number of iterations of Stability selection
#' @param mda_type value between "shuffle" or "zero" (weather to randomize a TF or put it to zero in 
#' feature importance estimation)
#' @param tf_expression_permutation weather or not to shuffle the expression of TFs between each other.
#' @param nCores Number of cores for multithreading
#'
#' @return a matrix of feature importances for each TF-target pairs
#' @export
#'
#' @examples
weightedLASSO_inference_MSE <- function(counts, genes, tfs, alpha=0.25, 
                                   pwm_occurrence, int_pwm_noise = 0,
                                   N = 100, mda_type= "shuffle", 
                                   tf_expression_permutation = FALSE,
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
  message(paste("\nUsing weightedLASSO with", foreach::getDoParWorkers(), "cores. alpha =", alpha))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(result.reg <-
                                   doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = c, 
                                                                     .final = function(x) {names(x) <- genes; x}, 
                                                                     .inorder = TRUE),
                                                    {
                                                      # getting rid of the TF variable on which the regression is made
                                                      # if needed
                                                      target_tfs <- setdiff(tfs, target)
                                                      x_target <- x[, target_tfs]
                                                      if(tf_expression_permutation){
                                                        # randomises the expression rows of TFs but not their ID
                                                        x_target <- x_target[,sample(target_tfs, replace = F, 
                                                                                     size = length(target_tfs))]
                                                        colnames(x_target) <- target_tfs
                                                      }
                                                      y <- t(counts[target, ])
                                                      # y_norm = (y-mean(y))/sd(y)
                                                      
                                                      # weights for differential shrinkage
                                                      penalty_factor <- 1 - pwm_imputed[target, target_tfs] * alpha
                                                      
                                                      importances <- setNames(rep(0, length(tfs)), tfs)
                                                      
                                                      ######### Stability Selection
                                                      mse_gene = 0
                                                      n_actual = 0
                                                      for(n in 1:N) {
                                                        # bootstrapping observations
                                                        sampled <- sample(1:nrow(x), replace = T, size = nrow(x))
                                                        oob <- setdiff(1:nrow(x), sampled)
                                                        
                                                        if(length(oob)>1){ #ensures a test set of sufficient size
                                                          
                                                          # perturbating differential shrinkage
                                                          # noisy_penalty_factor <- pmax(pmin(penalty_factor + 
                                                          #                                     runif(n=length(penalty_factor),
                                                          #                                           min = -int_pwm_noise*alpha, 
                                                          #                                           max = int_pwm_noise*alpha), 1),0)
                                                          
                                                          # models are try-catched because in rare cases glmnet
                                                          # crashes for convergence issues
                                                          tryCatch(
                                                            error = function(cnd) "cv.glmnet internal error due to convergence issues",
                                                            
                                                            # model training on sampled observations
                                                            # estimating lambda by a 5 fold CV
                                                            {mymodels_pen = cv.glmnet(
                                                              x_target[sampled,],
                                                              y[sampled],
                                                              maxit=maxit,
                                                              family = "poisson",
                                                              nfolds = 5,
                                                              penalty.factor = penalty_factor)
                                                            
                                                            # value of lambda 1se
                                                            ilambda.1se <- which(mymodels_pen$lambda == mymodels_pen$lambda.1se)
                                                            
                                                            # model predictions on OOB observations
                                                            y_hat <- exp(predict(mymodels_pen, newx = x_target[oob,],
                                                                                        type = "link", s= ilambda.1se))
                                                            # prediction MSE on OOB data
                                                            mse_oob <- mean((y_hat - y[oob])^2)
                                                            mse_gene <- mse_gene + mse_oob
                                                            
                                                            # selection frequency
                                                            n_actual = n_actual+1
                                                            }
                                                          )
                                                        }
                                                      }
                                                      # the mse is normalized by the gene variance
                                                      mse_gene/(n_actual*sd(y)^2)
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  edges <- result.reg
  toc()
  return(edges)
}






# 
# get_MSE_baseline <- function(counts, genes, nCores=ifelse(is.na(detectCores()),1,
#                                                           max(detectCores() - 1, 1))){
#   registerDoParallel(cores = nCores)
#   tic()
#   suppressPackageStartupMessages(mmse <-
#                                    doRNG::"%dorng%"(foreach::foreach(
#                                      target = genes,  
#                                      .inorder = TRUE),
#                                      {
#                                        sampled <- sample(1:ncol(counts), replace = T, size = ncol(counts))
#                                        oob <- setdiff(1:ncol(counts), sampled)
#                                        
#                                        mean((counts[target,oob] - mean(counts[target, sampled]))^2)
#                                      }
#                                    ))
#   attr(mmse, "rng") <- NULL
#   mmmse <- unlist(mmse)
#   return(mmmse)
# }
