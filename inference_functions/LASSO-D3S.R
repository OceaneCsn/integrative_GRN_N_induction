######################## imports for parallel version
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)


library(glmnet)


LASSO.D3S_inference <- function(counts, genes, tfs, alpha=0.25, 
                                pwm_occurrence, 
                               resampling_prop = 1, int_pwm_noise = 0.1,
                                nfolds.cv = 5, N = 100, lambda = "min", 
                               score = "pval", robustness = 0.8,
                                nCores = ifelse(is.na(detectCores()),1,
                                                max(detectCores() - 1, 1))){
  
  
  # to avoid convergence issues : "inner loop 3; cannot correct step size"
  if(alpha==1)
    alpha=1-1e-8
  
 counts <- round(counts, 0)
 x <- t(counts[tfs,])
  
  # pwm scores to bias variable selecttion toward pairs supported by a TFBS
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  # parallel computing of the lasso
  registerDoParallel(cores = nCores)
  message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
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
                                                      
                                                      for(n in 1:N) {
                                                        # bootstrapping observations
                                                        sampled <-sample(1:nrow(x), replace = T,
                                                                        size = round(resampling_prop*nrow(x), 0))
                                                        # perturbating PWM penalties
                                                        noisy_penalty_factor <- pmax(pmin(penalty_factor + 
                                                          runif(n=length(penalty_factor),
                                                                min = -int_pwm_noise*alpha, 
                                                                max = int_pwm_noise*alpha), 1),0)
                                                        # models are try-catched because in rare cases glmnet
                                                        # crashes for some obscure convergence issues
                                                        tryCatch(
                                                          error = function(cnd) "cv.glmnet internal error due to convergence issues",
                                                          
                                                          {mymodels_pen = cv.glmnet(
                                                          x_target[sampled,],
                                                          y[sampled],
                                                          family = "poisson",
                                                          nfolds = nfolds.cv,
                                                          penalty.factor = noisy_penalty_factor,
                                                          keep = TRUE)
                                                        
                                                        if(lambda=="min")
                                                          penalty = mymodels_pen$lambda.min
                                                        if(lambda=="1se")
                                                          penalty = mymodels_pen$lambda.1se
                                                        
                                                        iLambdaMin = which(mymodels_pen$lambda == penalty)
                                                        selected_tfs <- names(which(mymodels_pen$glmnet.fit$beta[, iLambdaMin] != 0))
                                                        # selection frequency
                                                        importances[selected_tfs] <- importances[selected_tfs]+1}
                                                        )
                                                      }
                                                      robust_tfs  = NULL
                                                      if(score=="freq")
                                                        to_return <- data.frame(from =names(importances),
                                                                   to = target,
                                                                   importance = importances/N)
                                                      
                                                      if(score=="pval"){
                                                        robust_tfs <- names(importances[importances>=N*robustness])
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
                                                      }
                                                      if(score=="freq" | score=="pval" & length(robust_tfs)>0)
                                                        to_return
                                                      
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  edges <- result.reg
  toc()
  return(edges)
}





LASSO.D3S_network <- function(mat, density, pwm_occurrence, genes, tfs){
  edges_ <- mat[order(mat$importance, decreasing = F),] 
  
  nEdges = round(density * (length(genes) - 1) * length(tfs), 0)
  edges_ <- edges_[1:nEdges,]
  
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  edges_$pwm <- pwm_imputed[cbind(edges_$to, edges_$from)]
  return(edges_[,c("from", "to", "pwm")])
}
