######################## imports
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)
library(glmnet)
library(tictoc)
library(ppcor)

# most faithful implementation of the method in this paper: https://genome.cshlp.org/content/early/2019/02/19/gr.238253.118.full.pdf+html
# lasso.stars had to be modified in order to allow for differential shrinkage
# include the modified stars stability selection function:
source("inference_functions/lasso.stars.R")


#' mLASSO-STARS GRN inference (closest R implementation possible!)
#'
#' @param counts Expression matrix (genes in rownames, conditions in columns)
#' @param genes Vector of genes (in the rownames of counts) to be used in GRN inference as target genes
#' @param tfs vector of genes (in the rownames of counts) that are transcriptional regulators
#' to be used a predictors in the regressions for GRN inference
#' @param alpha The strength of data integration.
#' Numeric value (e.g 0, 1) or a named vector giving the value of alpha for each target gene
#' @param pwm_occurrence Prior matrix Pi, giving PWM presence scores for TFs in rows
#' and genes in columns. Can contain NAs for TFs that do not have a PWM available.
#' @param N Number of iterations of Stability selection in Stars
#' @param tf_expression_permutation weather or not to shuffle the expression of TFs between each other.
#' @param nCores Number of cores for multithreading
#' @param family Type of distribution for the glm ("gaussian" or "poisson" (default))
#' @return a matrix of feature importances for each TF-target pairs
#' @export
#'
#' @examples
mLASSO_stars_inference <- function(counts, genes, tfs, alpha=0.25, 
                                    pwm_occurrence, stars.thresh = 0.05,
                                    N_boot = 5, N_stars_ss = 10, 
                                    family = "gaussian",
                                    nCores = ifelse(is.na(detectCores()),1,
                                                    max(detectCores() - 1, 1))){
  
  # for a gaussian lasso, data is log transformed
  if(family == "gaussian")
    counts <- log(counts+0.5)
  # else # for poisson glm, data needs to be integers
  #   counts <- round(counts, 0)
  
  # expression of regulators
  x <- t(counts[tfs,])
  
  # weather or not this is a gene-specific alpha model or not
  gene_specific = length(alpha) > 1
  
  # pwm scores to bias variable selection toward pairs supported by a TFBS
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  # partial correlation between genes, useful for feature importance
  partial_cor <- pcor(t(counts), method = "spearman")$estimate
  colnames(partial_cor) <- genes
  rownames(partial_cor) <- genes
  
  # parallel computing of the lasso
  registerDoParallel(cores = nCores)
  message(paste("\n weightedLASSO is running using", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(result.reg <-
                                   doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = cbind, 
                                                                     .final = function(x) {colnames(x) <- genes; x}, 
                                                                     .inorder = TRUE),
                                                    {
                                                      # getting rid of the TF variable on which the regression is made
                                                      # if needed
                                                      target_tfs <- setdiff(tfs, target)
                                                      x_target <- x[, target_tfs]
                                                      
                                                      y <- t(counts[target, ])
                                                      
                                                      if(gene_specific)
                                                        alpha_gene = alpha[target]
                                                      else
                                                        alpha_gene = alpha
                                                      
                                                      # to avoid convergence issues : "inner loop 3; cannot correct step size"
                                                      maxit = 1e+05
                                                      
                                                      if(alpha_gene==1){
                                                        alpha_gene=1-1e-4
                                                        # maxit = 1e+07
                                                      }
                                                      
                                                      # weights for differential shrinkage
                                                      penalty_factor <- 1 - pwm_imputed[target, target_tfs] * alpha_gene
                                                      
                                                      importances <- setNames(rep(0, length(tfs)), tfs)
                                                      selections <- setNames(rep(0, length(tfs)), tfs)
                                                      
                                                      ######### Stability Selection
                                                      n_actual = 0
                                                      mse_gene = c()
                                                      for(n in 1:N_boot) {
                                                        # bootstrapping observations
                                                        # sampled <- sample(1:nrow(x), replace = T, size = nrow(x))
                                                        
                                                        # bootstrapping observations and controlling that
                                                        # duplicated observations are in the same fold
                                                        nfolds.cv = 5
                                                        idx <- sample(1:length(y), replace = F)
                                                        folds_bg <- split(idx, ceiling(seq_along(idx)/(length(y)/nfolds.cv)))
                                                        breaks <- c(0,cumsum(lengths(folds_bg)))
                                                        
                                                        sampled <- rep(0, length(y))
                                                        foldid <- rep(0, length(y))
                                                        for(fold in 1:length(folds_bg)){
                                                          sampled[(breaks[fold]+1):(breaks[fold+1])] <- 
                                                            sample(folds_bg[[fold]], replace = T, size = lengths(folds_bg)[fold])
                                                          foldid[(breaks[fold]+1):(breaks[fold+1])] <- fold
                                                        }
                                                        oob <- setdiff(1:nrow(x), sampled)
                                                        
                                                        if(length(oob)>1){ #ensures a test set of sufficient size
                                                          
                                                          # models are try-catched because in rare cases glmnet
                                                          # crashes for convergence issues
                                                          tryCatch(
                                                            error = function(cnd) "cv.glmnet internal error due to convergence issues",
                                                            
                                                            # model training on sampled observations
                                                            {
                                                              mymodels_pen = lasso.stars(
                                                              x = scale(x_target[sampled,]),
                                                              y = scale(y[sampled]),
                                                              maxit=maxit,stars.thresh = stars.thresh,
                                                              rep.num = N_stars_ss,
                                                              penalty.factor = penalty_factor)
                                                            
                                                            # feature selection for optimal lambda
                                                            selected_tfs <- names(which(mymodels_pen$opt.beta != 0))
                                                            
                                                            
                                                            # model predictions on OOB observations
                                                            # home-made predict()
                                                            y_hat <- rowSums(t(t(scale(x_target[oob,]))*mymodels_pen$opt.beta)) + 
                                                              mymodels_pen$opt.intercept
                                                            
                                                            # prediction of MSE on OOB data
                                                            mse_oob <- mean((y_hat - scale(y[oob]))^2)
                                                            mse_gene <- c(mse_gene, mse_oob)
                                                            n_actual = n_actual+1
                                                            selections[selected_tfs] <- selections[selected_tfs]+1
                                                            
                                                            })
                                                        }
                                                      }
                                                      # feature importance is computed as the number of times a TF is
                                                      # selected + its partial correlation to y
                                                      for(tf in tfs){
                                                        importances[tf] <- selections[tf] + abs(partial_cor[target, tf])
                                                      }
                                                      # return value for one target gene
                                                      c(importances[tfs], setNames(median(mse_gene), "mse"))
                                                    }))
  attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
  edges <- result.reg
  toc()
  return(edges)
}





#' Threshold weighted LASSO GRN to a desired density
#'
#' @param mat result of weightedLASSO_inference function
#' @param density desired network density
#' @param pwm_occurrence matrix of prior data Pi containing TFBS scores between
#' TFs and genes
#' @param genes list of genes used as inputs for GRN inference
#' @param tfs list of TFs used as predictors for GRN inference
#'
#' @return dataframe of oriented edges, and their prior value in pwm_occurrence
mLASSO_stars_network <- function(mat, density, pwm_occurrence, genes, tfs, decreasing=FALSE){
  # getting the number of genes for a desired density
  nEdges = round(density * (length(genes) - 1) * length(tfs), 0)
  
  mat <- mat[!str_detect(rownames(mat), 'mse'),]
  # getting the ranked list of edges
  links <- getLinkListLasso(mat, reportMax = nEdges, decreasing=decreasing)
  network <- graph_from_data_frame(links, directed = T)
  edges <- as_long_data_frame(network)[c(4,5)]
  colnames(edges) <- c('from', 'to')
  
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
  edges$pwm <- pwm_imputed[cbind(edges$to, edges$from)]
  return(edges[,c("from", "to", "pwm")])
}



getLinkListLasso <- function (weightMatrix, reportMax = NULL, threshold = 0, decreasing=FALSE) 
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
  linkList <- linkList[order(linkList$weight, decreasing = decreasing), 
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