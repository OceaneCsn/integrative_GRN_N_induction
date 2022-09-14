######################## imports for parallel version
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)

library(igraph)
library(boot)
library(randomForest)
library(stringr)

get_MSE <- function(net, method, nCores=30){
  registerDoParallel(cores = nCores)
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
                                         
                                         if(length(target_tfs) > 0){
                                           x_target <- x[, target_tfs]
                                           y <- as.numeric(t(counts[target, ]))
                                           target_tfs <- str_replace_all(target_tfs, '-', '.')
                                           target <- str_replace_all(target, '-', '.')
                                           lm_target <- glm(
                                             formula = paste(paste0("`", target, "`"), '~', 
                                                             paste(paste0("`", target_tfs, "`"), 
                                                                   collapse = '+')), 
                                             data = data.frame(t(round(counts,0))), 
                                             family = "poisson")
                                           
                                           cv.glm(glmfit = lm_target, K=5, 
                                                  data = data.frame(t(round(counts,0))))$delta[2]
                                         }
                                       }))
    attr(mmse, "rng") <- NULL
    mmmse <- mean(log(unlist(mmse)), na.rm=TRUE)
  }
  return(mmmse)
}