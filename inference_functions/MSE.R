######################## imports for parallel version
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)

library(igraph)
library(boot)
library(randomForest)
library(stringr)

get_MSE <- function(net, method, nCores=ifelse(is.na(detectCores()),1,
                                               max(detectCores() - 1, 1))){
  registerDoParallel(cores = nCores)
  message(paste("MSE for method", method))
  tic()
  if(method == "bRF" ){
    x <- t(counts[tfs,])
    suppressPackageStartupMessages(mmse <-
                                     doRNG::"%dorng%"(foreach::foreach(
                                       target = genes,  
                                       .inorder = TRUE),
                                       {
                                         
                                         target_tfs <- net[net$to == target,"from"]
                                         if(length(target_tfs) > 0){
                                           x_target <- x[, target_tfs]
                                           y <- as.numeric(t(counts[target, ]))
                                           rf <- randomForest(x=data.frame(x_target),y=y,importance=TRUE,
                                                              mtry=round(sqrt(length(target_tfs))),ntree=1000)
                                           mean(rf$mse)
                                         }
                                       }))
    attr(mmse, "rng") <- NULL
    mmmse <- mean(log(unlist(mmse)))
  }
  if(method == "LASSO-D3S" ){
    x <- round(t(counts[tfs,]),0)
    
    suppressPackageStartupMessages(mmse <-
                                     doRNG::"%dorng%"(foreach::foreach(
                                       target = genes,  
                                       .inorder = TRUE),
                                       {
                                         target_tfs <- net[net$to == target,"from"]
                                         x_ <- round(t(counts[c(tfs, target),]),0)
                                         if(length(target_tfs) > 0){
                                           x_target <- x[, target_tfs]
                                           y <- as.numeric(t(counts[target, ]))
                                           target_tfs <- str_replace_all(target_tfs, '-', '.')
                                           target <- str_replace_all(target, '-', '.')
                                           lm_target <- glm(
                                             formula = paste(paste0("`", target, "`"), '~', 
                                                             paste(paste0("`", target_tfs, "`"), 
                                                                   collapse = '+')), 
                                             data = data.frame(x_), 
                                             family = "poisson")
                                           
                                           cv.glm(glmfit = lm_target, K=length(y), 
                                                  data = data.frame(t(round(counts,0))))$delta[2]
                                         }
                                       }))
    
    attr(mmse, "rng") <- NULL
    mmmse <- mean(log(unlist(mmse)), na.rm=TRUE)
  }
  toc()
  return(mmmse)
}


get_MSE_baseline <- function(counts, genes, nCores=ifelse(is.na(detectCores()),1,
                                               max(detectCores() - 1, 1))){
  registerDoParallel(cores = nCores)
  tic()
  suppressPackageStartupMessages(mmse <-
                                   doRNG::"%dorng%"(foreach::foreach(
                                     target = genes,  
                                     .inorder = TRUE),
                                     {
                                       sampled <- sample(1:ncol(counts), replace = T, size = ncol(counts))
                                       oob <- setdiff(1:ncol(counts), sampled)
                                       
                                       mean((counts[target,oob] - mean(counts[target, sampled]))^2)
                                       }
                                     ))
    attr(mmse, "rng") <- NULL
    mmmse <- mean(log(unlist(mmse)))
  return(mmmse)
}


#' bRF GRN inference returns MSE
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
bRF_inference_MSE <- function(counts, genes, tfs, alpha=0.25, scale = FALSE,
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

                                                        mean(rf_weighted$mse)/var(y)
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





LASSO.D3S_inference_MSE <- function(counts, genes, tfs, alpha=0.25, 
                                   pwm_occurrence, int_pwm_noise = 0,
                                   N = 100, mda_type= "shuffle", 
                                   seed =sample(1:10000, 1),
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
  message(paste("\nUsing LASSO-D3S with", foreach::getDoParWorkers(), "cores. alpha =", alpha))
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
                                                        set.seed(sample(1:10000, 1))
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
                                                      set.seed(seed)
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
                                                            y_hat <- exp(predict.glmnet(mymodels_pen, newx = x_target[oob,],
                                                                                        type = "link"))
                                                            
                                                            # getting lambda value granting minimal MSE 
                                                            # on OOB observations
                                                            mses_oob <- c()
                                                            for(s in colnames(y_hat)){
                                                              mses_oob <- c(mses_oob, mean((y[oob] - y_hat[,s])^2))
                                                            }
                                                            iLambdaMin = which(mses_oob == min(mses_oob))
                                                            
                                                            # y_hat_norm <- (y_hat[,iLambdaMin] - mean(y_hat[,iLambdaMin]))/
                                                            #   sd(y_hat[,iLambdaMin])
                                                            # 
                                                            # if(normalized)
                                                            #   mse_gene <- mse_gene + mean((y_norm[oob] - y_hat_norm[,iLambdaMin])^2)
                                                            # else
                                                            mse_gene <- mse_gene + min(mses_oob)
                                                            
                                                            # selection frequency
                                                            n_actual = n_actual+1
                                                            }
                                                          )
                                                        }
                                                      }
                                                      mse_gene/(n_actual*sd(y)^2)
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  edges <- result.reg
  toc()
  return(edges)
}

