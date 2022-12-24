######################## imports for parallel version
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)

###################### imports to use C++ code
# compiled files depend on the OS ...
if(.Platform$OS.type=="windows"){
  dyn.load("inference_functions/Cpp_dependencies/rfutils.dll")
  dyn.load("inference_functions/Cpp_dependencies/regTree.dll")
  dyn.load("inference_functions/Cpp_dependencies/regrf_mse.dll")
} else {
  dyn.load("inference_functions/Cpp_dependencies/rfutils.so")
  dyn.load("inference_functions/Cpp_dependencies/regTree.so")
  dyn.load("inference_functions/Cpp_dependencies/regrf_mse.so")
}


source('inference_functions/iRafNet_utils.R')

########## other imports
library(tictoc)
library(reshape2)
library(igraph)

#' bRF GRN inference
#'
#' @param counts expression matrix with gene IDs as rownames and conditions in columns
#' @param genes list of genes used as inputs for GRN inference
#' @param tfs list of TFs used as predictors for GRN inference
#' @param alpha integration strength, value should be between 0 and 1.
#' @param scale weather or not to scale expression data to z-scores.
#' @param prior_strength value of k in the link function
#' @param pwm_occurrence matrix of prior data Pi containing TFBS scores between
#' TFs and genes
#' @param nTrees Number of trees in Random Forests
#' @param importance Importance metric in Random Forests.
#' MDA ("%IncMSE") is the default value. MDI can be used via "IncNodePurity".
#' @param nCores Number of cores for multithreading. (Not supported on Windows)
#'
#' @return The weighted list of regulatory interactions between genes and TFs
bRF_inference <- function(counts, genes, tfs, alpha=0.25, scale = FALSE,
                          prior_strength = 2,
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
  message(paste("\nbRF is running using", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  force(irafnet_onetarget)
  suppressPackageStartupMessages(result.reg <-
                                   doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = cbind, 
                                                                     .final = function(x) {colnames(x) <- genes; x}, 
                                                                     .inorder = TRUE),
                                                    {
                                                      target_tfs <- setdiff(tfs, target)
                                                      x_target <- x[, target_tfs]
                                                      p = length(target_tfs)
                                                      y <- as.numeric(t(counts[target, ]))
                                                      # first version for Phd manuscript
                                                      #weights <- 10^(prior_strength* pwm_imputed[target, target_tfs]*alpha)
                                                      
                                                      # new version for a more similar behavior to the lasso for 
                                                      # large alphas and no need for k
                                                      weights <- ifelse(pwm_imputed[target, target_tfs] == 1, 1,
                                                                        ifelse(pwm_imputed[target, target_tfs] == 0.5, 1-alpha,
                                                                               -sqrt(1-(alpha-1)^2)+1))
                                                      # sophie's version
                                                      # weights <- ifelse(pwm_imputed[target, target_tfs] == 1, 1,
                                                      #                   ifelse(pwm_imputed[target, target_tfs] == 0, 1-alpha,
                                                      #                          sqrt(1-alpha^2)))
                                                      # do not include genes in which no known PWMs are found
                                                      # at alpha = 1
                                                      if(sum(weights)>0){
                                                        weights <- weights/sum(weights)
                                                        
                                                        
                                                        rf_weighted <-irafnet_onetarget(x_target,y=y,importance=TRUE,
                                                                                        mtry=round(sqrt(p)),
                                                                                        ntree=nTrees,
                                                                                        sw=weights)
                                                        
                                                        im <- rf_weighted$importance[,importance]
                                                        # get relative increase in MSE instead of absolute increase
                                                        if(importance == "%IncMSE"){
                                                          im <- im/mean(rf_weighted$mse)
                                                        }
                                                        c(setNames(0, target), setNames(im, names(im)))[tfs]
                                                      } else {
                                                        setNames(rep(0, length(tfs)), tfs)
                                                      }
                                                      
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  mat <- result.reg
  toc()
  return(mat)
}






#' Threshold bRF GRN to a desired density
#'
#' @param mat result of bRF_inference function
#' @param density desired network density
#' @param pwm_occurrence matrix of prior data Pi containing TFBS scores between
#' TFs and genes
#' @param genes list of genes used as inputs for GRN inference
#' @param tfs list of TFs used as predictors for GRN inference
#'
#' @return dataframe of oriented edges, and their prior value in pwm_occurrence
bRF_network <- function(mat, density, pwm_occurrence, genes, tfs){
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



getLinkList <- function (weightMatrix, reportMax = NULL, threshold = 0) 
{
  if (!is.numeric(threshold)) {
    stop("threshold must be a number.")
  }
  regulatorsInTargets <- rownames(weightMatrix)[rownames(weightMatrix) %in% 
                                                  colnames(weightMatrix)]
  if (length(regulatorsInTargets) == 1) 
    weightMatrix[regulatorsInTargets, regulatorsInTargets] <- NA
  if (length(regulatorsInTargets) > 1) 
    diag(weightMatrix[regulatorsInTargets, regulatorsInTargets]) <- NA
  linkList <- reshape2::melt(weightMatrix, na.rm = TRUE)
  colnames(linkList) <- c("regulatoryGene", "targetGene", "weight")
  linkList <- linkList[linkList$weight >= threshold, ]
  linkList <- linkList[order(linkList$weight, decreasing = TRUE), 
  ]
  if (!is.null(reportMax)) {
    linkList <- linkList[1:min(nrow(linkList), reportMax), 
    ]
  }
  rownames(linkList) <- NULL
  uniquePairs <- nrow(unique(linkList[, c("regulatoryGene", 
                                          "targetGene")]))
  if (uniquePairs < nrow(linkList)) 
    warning("There might be duplicated regulator-target (gene id/name) pairs.")
  return(linkList)
}
