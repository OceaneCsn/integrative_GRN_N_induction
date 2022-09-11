source('inference_functions/utils.R')
source('inference_functions/iRafNet_utils.R')


######################## imports for parallel version
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)

###################### imports to use C++ code
dyn.load("inference_functions/Cpp_dependencies/rfutils.so")
dyn.load("inference_functions/Cpp_dependencies/regTree.so")
dyn.load("inference_functions/Cpp_dependencies/regrf.so")
dyn.load("inference_functions/Cpp_dependencies/regrf_mse.so")


########## other imports
library(tictoc)
library(GENIE3)
library(igraph)

bRF_inference <- function(counts, genes, tfs, alpha=0.25, scale = FALSE,
                          prior_strength = 5,
                          pwm_occurrence, nTrees=500, importance="%IncMSE",
                          nCores = ifelse(is.na(detectCores()),1,
                                          max(detectCores() - 1, 1))){
  
  # counts must be normalized so that genes to have comparable node purities
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
  message(paste("\nCustom iRafNet Using", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(result.reg <-
                                   doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = cbind, 
                                                                     .final = function(x) {colnames(x) <- genes; x}, 
                                                                     .inorder = TRUE),
                                                    {
                                                      target_tfs <- setdiff(tfs, target)
                                                      x_target <- x[, target_tfs]
                                                      p = length(target_tfs)
                                                      y <- as.numeric(t(counts[target, ]))
                                                      weights <- 10^(prior_strength* pwm_imputed[target, target_tfs]*alpha)
                                                      weights <- weights/sum(weights)
                                                      
                                                      
                                                      rf_weighted <-irafnet_onetarget(x_target,y=y,importance=TRUE,
                                                                                      mtry=round(sqrt(p)),
                                                                                      ntree=nTrees,
                                                                                      sw=weights)
                                                      
                                                      im <- rf_weighted$importance[,importance]
                                                      
                                                      c(setNames(0, target), setNames(im, names(im)))[tfs]
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  mat <- result.reg
  toc()
  return(mat)
}






bRF_network <- function(mat, density, pwm_occurrence, genes, tfs){
  # DIANE's function network_thresholding
  
  # getting the number of genes for a desired density
  nEdges = round(density * (length(genes) - 1) * length(tfs), 0)
  
  # getting the ranked list of edges
  links <- getLinkList(mat, reportMax = nEdges)
  network <- graph_from_data_frame(links, directed = T)
  
  
  edges <- as_long_data_frame(network)[c(4,5)]
  colnames(edges) <- c('from', 'to')
  
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  edges$pwm <- pwm_imputed[cbind(edges$to, edges$from)]
  return(edges[,c("from", "to", "pwm")])
}

