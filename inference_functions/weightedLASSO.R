######################## imports
library(doParallel)
library(parallel)
library(foreach)
library(doRNG)
library(glmnet)
library(tictoc)




#' weightedLASSO GRN inference
#'
#' @param counts Expression matrix (genes in rownames, conditions in columns)
#' @param genes Vector of genes (in the rownames of counts) to be used in GRN inference as target genes
#' @param tfs vector of genes (in the rownames of counts) that are transcriptional regulators
#' to be used a predictors in the regressions for GRN inference
#' @param alpha The strength of data integration.
#' Numeric value (e.g 0, 1) or a named vector giving the value of alpha for each target gene
#' @param pwm_occurrence Prior matrix Pi, giving PWM presence scores for TFs in rows
#' and genes in columns. Can contain NAs for TFs that do not have a PWM available.
#' @param int_pwm_noise Random perturbation applied to PWM priors in pwm_occurrence.
#' Defalut is none, experimental.
#' @param N Number of iterations of Stability selection
#' @param mda_type value between "shuffle" or "zero" (weather to randomize a TF or put it to zero in 
#' feature importance estimation)
#' @param tf_expression_permutation weather or not to shuffle the expression of TFs between each other.
#' @param nCores Number of cores for multithreading
#' @param nfolds.cv Number of folds for cross validation
#' @param family Type of distribution for the glm ("gaussian" or "poisson" (default))
#' @param EN_param ElasticNet parameter. 1 (default) is LASSO, 0 is Rigde.
#' @param lambda "min" ou "1se" (default)
#' @return a matrix of feature importances for each TF-target pairs
#' @export
#'
#' @examples
weightedLASSO_inference <- function(counts, genes, tfs, alpha=0.25, 
                                pwm_occurrence, int_pwm_noise = 0,
                                N = 100, mda_type="shuffle",
                                family = "poisson", EN_param = 1,
                                tf_expression_permutation = FALSE,
                                nfolds.cv=5, lambda = "1se",
                                nCores = ifelse(is.na(detectCores()),1,
                                                max(detectCores() - 1, 1))){
  
  
  
  
  # for a gaussian lasso, data is log transformed
  if(family == "gaussian")
    counts <- log(counts+0.5)
  else # for poisson glm, data needs to be integers
    counts <- round(counts, 0)
  
  # expression of regulators
  x <- t(counts[tfs,])
  
  # weather or not this is a gene-specific alpha model or not
  gene_specific = length(alpha) > 1
  
  # pwm scores to bias variable selection toward pairs supported by a TFBS
  pwm_imputed <- pwm_occurrence
  pwm_imputed[is.na(pwm_imputed)] <- 0.5
  
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
                                                      if(tf_expression_permutation){
                                                        # randomises the expression rows of TFs but not their ID
                                                        # this is an in-silico null hypothesis 
                                                        x_target <- x_target[,sample(target_tfs, replace = F, 
                                                                                     size = length(target_tfs))]
                                                        colnames(x_target) <- target_tfs
                                                      }
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
                                                      
                                                      ######### Stability Selection / bootstraps
                                                      n_actual = 0
                                                      mse_gene = c()
                                                      for(n in 1:N) {

                                                        # bootstrapping observations and controlling that
                                                        # duplicated observations are in the same fold
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
                                                            {mymodels_pen = cv.glmnet(
                                                              x_target[sampled,],
                                                              y[sampled],
                                                              maxit=maxit,
                                                              alpha = EN_param,
                                                              family = family,
                                                              nfolds = nfolds.cv,
                                                              foldid = foldid,
                                                              penalty.factor = penalty_factor)
                                                            
                                                            # value of lambda
                                                            if(lambda == "1se")
                                                              ilambda <- which(mymodels_pen$lambda == mymodels_pen$lambda.1se)
                                                            else ilambda <- which(mymodels_pen$lambda == mymodels_pen$lambda.min)
                                                            
                                                            # feature selection for optimal lambda
                                                            selected_tfs <- names(which(mymodels_pen$glmnet.fit$beta[, ilambda] != 0))
                                                            
                                                            # model predictions on OOB observations
                                                            if(family == "poisson"){
                                                              y_hat <- exp(predict(mymodels_pen, newx = x_target[oob,],
                                                                                   type = "link", s= ilambda))
                                                            }
                                                            else{
                                                              y_hat <- predict(mymodels_pen, newx = x_target[oob,],type = "response",
                                                                                  s= paste0("lambda.", lambda))
                                                            }
                                                            
                                                            # prediction of MSE on OOB data
                                                            mse_oob <- mean((y_hat - y[oob])^2)
                                                            mse_gene <- c(mse_gene, mse_oob)
                                                            n_actual = n_actual+1
                                                            
                                                            # computing importance on OOB conditions
                                                            for(sel_tf in selected_tfs){
                                                              x_target_rand <- x_target[oob,]

                                                              # randomizes the conditions to break the link
                                                              # between regulator and target
                                                              if(mda_type == "shuffle")
                                                                x_target_rand[,sel_tf] <- sample(x_target_rand[,sel_tf],
                                                                                                 size = nrow(x_target_rand),
                                                                                                 replace = F)
                                                              # alternative mda : removes the variable from predictors
                                                              if(mda_type == "zero")
                                                                x_target_rand[,sel_tf] <- 0
                                                              
                                                              
                                                              # makes predictions on the permuted data
                                                              if(family == "poisson"){
                                                                y_hat_rand <- exp(predict(mymodels_pen, newx = x_target_rand,
                                                                                          type = "link", s = ilambda))
                                                              }
                                                              else{
                                                                y_hat <- predict(mymodels_pen, newx = x_target_rand,
                                                                                 s= paste0("lambda.", lambda))
                                                              }
                                                              
                                                              # mse on this shuffled dataset
                                                              mse_rand <- mean((y_hat_rand - y[oob])^2)
                                                              
                                                              # importance is the normalized MSE increase induced by the shuffling of
                                                              # a given variable
                                                              importances[sel_tf] <- importances[sel_tf] +
                                                                max(0,(mse_rand - mse_oob)/mse_rand)
                                                            }
                                                            }
                                                          )
                                                        }
                                                      }
                                                      # returns importance values and median mse in the different bootstraps
                                                      c(importances[tfs]/N, setNames(median(mse_gene)/(sd(y)^2), "mse"))
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
weightedLASSO_network <- function(mat, density, pwm_occurrence, genes, tfs, decreasing=FALSE){
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




# older versions

# LASSO.D3S_inferenceMDA <- function(counts, genes, tfs, alpha=0.25, 
#                                    pwm_occurrence, int_pwm_noise = 0,
#                                    N = 100, mda_type= "shuffle",
#                                    tf_expression_permutation = FALSE,
#                                    nCores = ifelse(is.na(detectCores()),1,
#                                                    max(detectCores() - 1, 1))){
#   
#   
#   # to avoid convergence issues : "inner loop 3; cannot correct step size"
#   maxit = 1e+05
#   if(alpha==1){
#     alpha=1-1e-8
#     maxit = 1e+07
#   }
#   
#   
#   counts <- round(counts, 0)
#   x <- t(counts[tfs,])
#   
#   # pwm scores to bias variable selecttion toward pairs supported by a TFBS
#   pwm_imputed <- pwm_occurrence
#   pwm_imputed[is.na(pwm_imputed)] <- 0.5
#   
#   # parallel computing of the lasso
#   registerDoParallel(cores = nCores)
#   message(paste("\nUsing LASSO-D3S with", foreach::getDoParWorkers(), "cores."))
#   "%dopar%" <- foreach::"%dopar%"
#   tic()
#   suppressPackageStartupMessages(result.reg <-
#                                    doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = cbind, 
#                                                                      .final = function(x) {colnames(x) <- genes; x}, 
#                                                                      .inorder = TRUE),
#                                                     {
#                                                       # getting rid of the TF variable on which the regression is made
#                                                       # if needed
#                                                       target_tfs <- setdiff(tfs, target)
#                                                       x_target <- x[, target_tfs]
#                                                       if(tf_expression_permutation){
#                                                         # randomises the expression rows of TFs but not their ID
#                                                         x_target <- x_target[,sample(target_tfs, replace = F, 
#                                                                                      size = length(target_tfs))]
#                                                         colnames(x_target) <- target_tfs
#                                                       }
#                                                       y <- t(counts[target, ])
#                                                       
#                                                       # weights for differential shrinkage
#                                                       penalty_factor <- 1 - pwm_imputed[target, target_tfs] * alpha
#                                                       
#                                                       importances <- setNames(rep(0, length(tfs)), tfs)
#                                                       
#                                                       ######### Stability Selection
#                                                       n_actual = 0
#                                                       for(n in 1:N) {
#                                                         # bootstrapping observations
#                                                         sampled <- sample(1:nrow(x), replace = T, size = nrow(x))
#                                                         oob <- setdiff(1:nrow(x), sampled)
#                                                         
#                                                         if(length(oob)>1){ #ensures a test set of sufficient size
#                                                           
#                                                           # perturbating differential shrinkage
#                                                           noisy_penalty_factor <- pmax(pmin(penalty_factor + 
#                                                                                               runif(n=length(penalty_factor),
#                                                                                                     min = -int_pwm_noise*alpha, 
#                                                                                                     max = int_pwm_noise*alpha), 1),0)
#                                                           
#                                                           # models are try-catched because in rare cases glmnet
#                                                           # crashes for convergence issues
#                                                           tryCatch(
#                                                             error = function(cnd) "cv.glmnet internal error due to convergence issues",
#                                                             
#                                                             # model training on sampled observations
#                                                             {mymodels_pen = glmnet(
#                                                               x_target[sampled,],
#                                                               y[sampled],
#                                                               maxit=maxit,
#                                                               family = "poisson",
#                                                               penalty.factor = noisy_penalty_factor)
#                                                             
#                                                             # model predictions on OOB observations
#                                                             y_hat <- exp(predict.glmnet(mymodels_pen, newx = x_target[oob,],
#                                                                                         type = "link"))
#                                                             
#                                                             # getting lambda value granting minimal MSE 
#                                                             # on OOB observations
#                                                             mses_oob <- c()
#                                                             for(s in colnames(y_hat)){
#                                                               mses_oob <- c(mses_oob, mean((y[oob] - y_hat[,s])^2))
#                                                             }
#                                                             iLambdaMin = which(mses_oob == min(mses_oob))
#                                                             mse_initial <- min(mses_oob)
#                                                             
#                                                             # feature selection for optimal lambda
#                                                             selected_tfs <- names(which(mymodels_pen$beta[, iLambdaMin] != 0))
#                                                             
#                                                             # computing MDA on OOB data
#                                                             for(sel_tf in selected_tfs){
#                                                               x_target_rand <- x_target[oob,]
#                                                               
#                                                               # randomizes the conditions to break the link
#                                                               # between regulator and target
#                                                               if(mda_type == "shuffle")
#                                                                 x_target_rand[,sel_tf] <- sample(x_target_rand[,sel_tf], 
#                                                                                                  size = nrow(x_target_rand), 
#                                                                                                  replace = F)
#                                                               # alternative mda : removes the variable from predictors
#                                                               if(mda_type == "zero")
#                                                                 x_target_rand[,sel_tf] <- 0
#                                                               
#                                                               y_hat_rand <- exp(predict.glmnet(mymodels_pen, newx = x_target_rand,
#                                                                                                type = "link", s = iLambdaMin))
#                                                               
#                                                               mse_rand <- mean((y_hat_rand - y[oob])^2)
#                                                               
#                                                               importances[sel_tf] <- importances[sel_tf] + 
#                                                                 (mse_rand - mse_initial)/mse_initial
#                                                             }
#                                                             
#                                                             # selection frequency
#                                                             n_actual = n_actual+1
#                                                             }
#                                                           )
#                                                         }
#                                                       }
#                                                       importances / n_actual
#                                                     }))
#   attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
#   edges <- result.reg
#   toc()
#   return(edges)
# }




# 
# 
# LASSO.D3S_inference <- function(counts, genes, tfs, alpha=0.25, 
#                                 pwm_occurrence, int_pwm_noise = 0,
#                                 N = 100, 
#                                 tf_expression_permutation = FALSE,
#                                 score = "pval", robustness = 0.1,
#                                 nCores = ifelse(is.na(detectCores()),1,
#                                                 max(detectCores() - 1, 1))){
#   
#   
#   
#   
#   
#   counts <- round(counts, 0)
#   x <- t(counts[tfs,])
#   
#   gene_specific = length(alpha) > 1
#   
#   # pwm scores to bias variable selection toward pairs supported by a TFBS
#   pwm_imputed <- pwm_occurrence
#   pwm_imputed[is.na(pwm_imputed)] <- 0.5
#   
#   # parallel computing of the lasso
#   registerDoParallel(cores = nCores)
#   message(paste("\nUsing LASSO-D3S with", foreach::getDoParWorkers(), "cores."))
#   "%dopar%" <- foreach::"%dopar%"
#   tic()
#   suppressPackageStartupMessages(result.reg <-
#                                    doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = cbind, 
#                                                                      .final = function(x) {colnames(x) <- genes; x}, 
#                                                                      .inorder = TRUE),
#                                                     {
#                                                       # getting rid of the TF variable on which the regression is made
#                                                       # if needed
#                                                       target_tfs <- setdiff(tfs, target)
#                                                       x_target <- x[, target_tfs]
#                                                       if(tf_expression_permutation){
#                                                         # randomises the expression rows of TFs but not their ID
#                                                         x_target <- x_target[,sample(target_tfs, replace = F, 
#                                                                                      size = length(target_tfs))]
#                                                         colnames(x_target) <- target_tfs
#                                                       }
#                                                       y <- t(counts[target, ])
#                                                       
#                                                       if(gene_specific)
#                                                         alpha_gene = alpha[target]
#                                                       else
#                                                         alpha_gene = alpha
#                                                       
#                                                       # to avoid convergence issues : "inner loop 3; cannot correct step size"
#                                                       maxit = 1e+05
#                                                       
#                                                       if(alpha_gene==1){
#                                                         alpha=1-1e-8
#                                                         maxit = 1e+07
#                                                       }
#                                                       
#                                                       # weights for differential shrinkage
#                                                       penalty_factor <- 1 - pwm_imputed[target, target_tfs] * alpha_gene
#                                                       
#                                                       importances <- setNames(rep(0, length(tfs)), tfs)
#                                                       
#                                                       ######### Stability Selection
#                                                       n_actual = 0
#                                                       for(n in 1:N) {
#                                                         # bootstrapping observations
#                                                         sampled <- sample(1:nrow(x), replace = T, size = nrow(x))
#                                                         oob <- setdiff(1:nrow(x), sampled)
#                                                         
#                                                         if(length(oob)>1){ #ensures a test set of sufficient size
#                                                           
#                                                           # perturbating differential shrinkage
#                                                           noisy_penalty_factor <- pmax(pmin(penalty_factor + 
#                                                                                               runif(n=length(penalty_factor),
#                                                                                                     min = -int_pwm_noise*alpha, 
#                                                                                                     max = int_pwm_noise*alpha), 1),0)
#                                                           
#                                                           # models are try-catched because in rare cases glmnet
#                                                           # crashes for convergence issues
#                                                           tryCatch(
#                                                             error = function(cnd) "cv.glmnet internal error due to convergence issues",
#                                                             
#                                                             # model training on sampled observations
#                                                             {mymodels_pen = glmnet(
#                                                               x_target[sampled,],
#                                                               y[sampled],
#                                                               maxit=maxit,
#                                                               family = "poisson",
#                                                               penalty.factor = noisy_penalty_factor)
#                                                             
#                                                             # model predictions on OOB observations
#                                                             y_hat <- exp(predict.glmnet(mymodels_pen, newx = x_target[oob,],
#                                                                                     type = "link"))
#                                                             
#                                                             # getting lambda value granting minimal MSE 
#                                                             # on OOB observations
#                                                             mses_oob <- c()
#                                                             for(s in colnames(y_hat)){
#                                                               mses_oob <- c(mses_oob, mean((y[oob] - y_hat[,s])^2))
#                                                             }
#                                                             iLambdaMin = which(mses_oob == min(mses_oob))
#                                                             
#                                                             # feature selection for optimal lambda
#                                                             selected_tfs <- names(which(mymodels_pen$beta[, iLambdaMin] != 0))
#                                                             
#                                                             # selection frequency
#                                                             importances[selected_tfs] <- importances[selected_tfs]+1
#                                                             n_actual = n_actual+1
#                                                             }
#                                                           )
#                                                         }
#                                                       }
#                                                       
#                                                       
#                                                       
#                                                       robust_tfs  = NULL
#                                                       if(score=="freq")
#                                                         to_return <- importances/n_actual
#                                                       ###### Unpenalized regression on robust TFs
#                                                       if(score=="pval"){
#                                                         to_return <- setNames(rep(1, length(tfs)), tfs)
#                                                         tryCatch(
#                                                           error = function(cnd) "a problem occurred in the unpenalized regressions for this gene.",
#                                                           
#                                                           # model training on sampled observations
#                                                           {robust_tfs <- names(importances[importances>=n_actual*robustness])
#                                                           # non penalized models with only robustly selected TFs
#                                                           if(length(robust_tfs)>0){
#                                                             lm_target <- glm(
#                                                               formula = paste(paste0("`", target, "`"), '~', 
#                                                                               paste(paste0("`", robust_tfs, "`"), 
#                                                                                     collapse = '+')), 
#                                                               data = data.frame(t(counts), check.names = F), 
#                                                               family = "poisson")
#                                                             
#                                                             pvals <- c(setNames(rep(1, length(setdiff(tfs, robust_tfs))), 
#                                                                                 nm = setdiff(tfs, robust_tfs)), 
#                                                                        setNames(summary(lm_target)$coefficients[
#                                                                          2:nrow(summary(lm_target)$coefficients),"Pr(>|z|)"], 
#                                                                          robust_tfs))
#                                                             #TF-target interactions are weighted by the TF pvalue in the glm
#                                                             to_return <-  pvals
#                                                           }
#                                                           })
#                                                       }
#                                                       if(score=="freq" | score=="pval" & length(robust_tfs)>0)
#                                                         to_return[tfs]
#                                                       
#                                                     }))
#   attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
#   edges <- result.reg
#   toc()
#   return(edges)
# }

# 
# 


#' 
#' 
#' #' weightedLASSO GRN inference
#' #'
#' #' @param counts Expression matrix (genes in rownames, conditions in columns)
#' #' @param genes Vector of genes (in the rownames of counts) to be used in GRN inference as target genes
#' #' @param tfs vector of genes (in the rownames of counts) that are transcriptional regulators
#' #' to be used a predictors in the regressions for GRN inference
#' #' @param alpha The strength of data integration.
#' #' Numeric value (e.g 0, 1) or a named vector giving the value of alpha for each target gene
#' #' @param pwm_occurrence Prior matrix Pi, giving PWM presence scores for TFs in rows
#' #' and genes in columns. Can contain NAs for TFs that do not have a PWM available.
#' #' @param int_pwm_noise Random perturbation applied to PWM priors in pwm_occurrence.
#' #' Defalut is none, experimental.
#' #' @param N Number of iterations of Stability selection
#' #' @param mda_type value between "shuffle" or "zero" (weather to randomize a TF or put it to zero in 
#' #' feature importance estimation)
#' #' @param tf_expression_permutation weather or not to shuffle the expression of TFs between each other.
#' #' @param robustness Rate of selection for a TF to be considered for feature importance estimation
#' #' @param nCores Number of cores for multithreading
#' #'
#' #' @return a matrix of feature importances for each TF-target pairs
#' #' @export
#' #'
#' #' @examples
#' weightedLASSO_inference_stable <- function(counts, genes, tfs, alpha=0.25, 
#'                                            pwm_occurrence, int_pwm_noise = 0,
#'                                            N = 100, mda_type="shuffle",
#'                                            tf_expression_permutation = FALSE,
#'                                            robustness = 0.2,
#'                                            nCores = ifelse(is.na(detectCores()),1,
#'                                                            max(detectCores() - 1, 1))){
#'   
#'   
#'   counts <- round(counts, 0)
#'   x <- t(counts[tfs,])
#'   
#'   gene_specific = length(alpha) > 1
#'   
#'   # pwm scores to bias variable selection toward pairs supported by a TFBS
#'   pwm_imputed <- pwm_occurrence
#'   pwm_imputed[is.na(pwm_imputed)] <- 0.5
#'   
#'   # parallel computing of the lasso
#'   registerDoParallel(cores = nCores)
#'   message(paste("\n weightedLASSO is running using", foreach::getDoParWorkers(), "cores."))
#'   "%dopar%" <- foreach::"%dopar%"
#'   tic()
#'   suppressPackageStartupMessages(result.reg <-
#'                                    doRNG::"%dorng%"(foreach::foreach(target = genes, .combine = cbind, 
#'                                                                      .final = function(x) {colnames(x) <- genes; x}, 
#'                                                                      .inorder = TRUE),
#'                                                     {
#'                                                       # getting rid of the TF variable on which the regression is made
#'                                                       # if needed
#'                                                       target_tfs <- setdiff(tfs, target)
#'                                                       x_target <- x[, target_tfs]
#'                                                       if(tf_expression_permutation){
#'                                                         # randomises the expression rows of TFs but not their ID
#'                                                         x_target <- x_target[,sample(target_tfs, replace = F, 
#'                                                                                      size = length(target_tfs))]
#'                                                         colnames(x_target) <- target_tfs
#'                                                       }
#'                                                       y <- t(counts[target, ])
#'                                                       
#'                                                       if(gene_specific)
#'                                                         alpha_gene = alpha[target]
#'                                                       else
#'                                                         alpha_gene = alpha
#'                                                       
#'                                                       # to avoid convergence issues : "inner loop 3; cannot correct step size"
#'                                                       maxit = 1e+05
#'                                                       
#'                                                       if(alpha_gene==1){
#'                                                         alpha=1-1e-8
#'                                                         maxit = 1e+07
#'                                                       }
#'                                                       
#'                                                       # weights for differential shrinkage
#'                                                       penalty_factor <- 1 - pwm_imputed[target, target_tfs] * alpha_gene
#'                                                       
#'                                                       importances <- setNames(rep(0, length(tfs)), tfs)
#'                                                       
#'                                                       ######### Stability Selection
#'                                                       n_actual = 0
#'                                                       for(n in 1:N) {
#'                                                         # bootstrapping observations
#'                                                         sampled <- sample(1:nrow(x), replace = T, size = nrow(x))
#'                                                         oob <- setdiff(1:nrow(x), sampled)
#'                                                         
#'                                                         if(length(oob)>1){ #ensures a test set of sufficient size
#'                                                           
#'                                                           # perturbating differential shrinkage
#'                                                           # noisy_penalty_factor <- pmax(pmin(penalty_factor + 
#'                                                           #                                     runif(n=length(penalty_factor),
#'                                                           #                                           min = -int_pwm_noise*alpha, 
#'                                                           #                                           max = int_pwm_noise*alpha), 1),0)
#'                                                           
#'                                                           # models are try-catched because in rare cases glmnet
#'                                                           # crashes for convergence issues
#'                                                           tryCatch(
#'                                                             error = function(cnd) "cv.glmnet internal error due to convergence issues",
#'                                                             
#'                                                             # model training on sampled observations
#'                                                             {mymodels_pen = cv.glmnet(
#'                                                               x_target[sampled,],
#'                                                               y[sampled],
#'                                                               maxit=maxit,
#'                                                               family = "poisson",
#'                                                               nfolds = 5,
#'                                                               penalty.factor = penalty_factor)
#'                                                             
#'                                                             # value of lambda 1se
#'                                                             ilambda.1se <- which(mymodels_pen$lambda == mymodels_pen$lambda.1se)
#'                                                             
#'                                                             # feature selection for optimal lambda
#'                                                             selected_tfs <- names(which(mymodels_pen$glmnet.fit$beta[, ilambda.1se] != 0))
#'                                                             
#'                                                             # selection frequency
#'                                                             # importances[selected_tfs] <- importances[selected_tfs]+1
#'                                                             
#'                                                             mse_initial <- mean((exp(predict(mymodels_pen, newx = x_target[oob,],
#'                                                                                              type = "link", s = ilambda.1se))-y[oob])^2)
#'                                                             
#'                                                             
#'                                                             
#'                                                             n_actual = n_actual+1
#'                                                             }
#'                                                             
#'                                                           )
#'                                                           
#'                                                           
#'                                                         }
#'                                                       }
#'                                                       
#'                                                       ###### Unpenalized regression on robust TFs
#'                                                       to_return <- setNames(rep(0, length(tfs)), tfs)
#'                                                       tryCatch(
#'                                                         error = function(cnd) "a problem occurred in the unpenalized regressions for this gene.",
#'                                                         
#'                                                         # model training on sampled observations
#'                                                         {
#'                                                           robust_tfs <- names(importances[importances>=n_actual*robustness])
#'                                                           # non penalized models with only robustly selected TFs
#'                                                           if(length(robust_tfs)>0){
#'                                                             lm_target <- glm(
#'                                                               formula = paste(paste0("`", target, "`"), '~',
#'                                                                               paste(paste0("`", robust_tfs, "`"),
#'                                                                                     collapse = '+')),
#'                                                               data = data.frame(t(counts), check.names = F),
#'                                                               family = "poisson")
#'                                                             
#'                                                             mse <- mean((lm_target$fitted.values - y)^2)
#'                                                             
#'                                                             for(sel_tf in robust_tfs){
#'                                                               x_target_rand <- x_target
#'                                                               
#'                                                               # randomizes the conditions to break the link
#'                                                               # between regulator and target
#'                                                               if(mda_type == "shuffle")
#'                                                                 x_target_rand[,sel_tf] <- sample(x_target_rand[,sel_tf],
#'                                                                                                  size = nrow(x_target_rand),
#'                                                                                                  replace = F)
#'                                                               # alternative mda : removes the variable from predictors
#'                                                               if(mda_type == "zero")
#'                                                                 x_target_rand[,sel_tf] <- 0
#'                                                               
#'                                                               y_hat_rand <- exp(predict.glm(lm_target, newdata = data.frame(x_target_rand),
#'                                                                                             type = "link"))
#'                                                               
#'                                                               mse_rand <- mean((y_hat_rand - y)^2)
#'                                                               
#'                                                               to_return[sel_tf] <- to_return[sel_tf] +
#'                                                                 (mse_rand - mse )/ mse_rand
#'                                                               
#'                                                             }
#'                                                           }
#'                                                         })
#'                                                       
#'                                                       importances[tfs]/n_actual
#'                                                     }))
#'   attr(result.reg, "rng") <- NULL # It contains the whole sequence of RNG seeds
#'   edges <- result.reg
#'   toc()
#'   return(edges)
#' }

