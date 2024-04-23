source('inference_functions/weightedRF.R')
source('inference_functions/weightedLASSO.R')
source('inference_functions/data_integration_optimization.R')


#' DIOgene : Data Integration Optimization for gene networks
#' 
#' This function allows the inference of integrative GRNs with a data driven gene-specific
#' contribution of the prior data.
#' 
#' It optimizes the contribution of prior data (here TFBMs)
#' to gene expression, using the maximal improvement of MSE over a 
#' simulated null hypothesis, for each gene.
#'
#' @param counts Expression matrix (genes in rownames, conditions in columns)
#' @param genes Vector of genes (must be present in the rownames of counts) 
#' to be used in GRN inference as target genes
#' @param tfs vector of genes (must be present in the rownames of counts) 
#' that are transcriptional regulators to be used as predictors in the regressions for GRN inference
#' @param prior_matrix Prior matrix Pi, a score (0, 0.5, 1) for each TF-target interaction.
#' Can contain NAs for TFs that do not have a PWM available. NAs will be turned into priors of 1/2.
#' @param model "linear" or "non-linear": the type of regression that must be performed.
#' "linear is for weightedLASSO, "non-linear" is for weightedRF.
#' @param nrep Number of repetitions of each regression model at each alpha value, 
#' to estimate MSE mean and dispersion on the true and null datasets. Must be
#' large enough (100 'default) works well on our case study).
#' @param nCores Number of cores for multithreading. (Not supported on Windows)
#' @param quiet Print individual model estimation times? 
#' @param ALPHAS Set of alpha values to be explored. Default value:
#' from 0 to 1 with a step of 0.1.
#' @param S Number of bootstrapped samples for robust regression models estimation
#'
#' @return A list containing:
#' - mats: the list of importance matrices for each alpha value and each replicate,
#' on the true and shuffled datasets
#' - mse: the MSE values for each gene, alpha value and replicate on the true and shuffled datasets
#' - edges: the list of inferred edges for each alpha value, replicate, and three densities
#' on the true and shuffled datasets
#' - alphas_opt: the list of optimal alpha values per gene.
#' - diogene: the importance matrix of the gene-specific models with optimal alpha.
DIOgene <- function(counts, genes, tfs, prior_matrix,
                    model = "non-linear",
                    nrep = 100,  nCores = 45,
                    quiet = T,
                    ALPHAS = seq(0,1, by = 0.1),
                    S = NULL){
  
  # number of bootstrapped models depends on the model type
  if(is.null(S)){
    if(model == "linear") S = 50
    if(model == "non-linear") S = 2000
  }
  
  # looping over alpha and replicates to estimate the MSE and EDI curves as
  # a function of alpha
  mats <- list()
  for(alpha in ALPHAS){ # exploring PWM integration strength
    print(paste("------------------- Running models for alpha =", alpha))
    for(rep in 1:nrep){ # exploring inherent variability
      # weightedLASSO
      if(model == "linear"){
        mat_lasso <- weightedLASSO_inference(counts = counts, genes = genes, tfs = tfs,
                                             alpha = alpha, N = S, 
                                             tf_expression_permutation = FALSE,
                                             mda_type = "shuffle",
                                             quiet = quiet,
                                             pwm_occurrence = prior_matrix,
                                             nCores = nCores)
        
        mat_lasso_perm <- weightedLASSO_inference(counts, genes, tfs,
                                                  alpha = alpha, N = S, 
                                                  tf_expression_permutation = TRUE,
                                                  mda_type = "shuffle",
                                                  quiet = quiet,
                                                  pwm_occurrence = prior_matrix,
                                                  nCores = nCores)
        
        mats[[paste0("LASSO_", as.character(alpha),  '_trueData_', rep)]] <- mat_lasso
        mats[[paste0("LASSO_", as.character(alpha),  '_shuffled_', rep)]] <- mat_lasso_perm
      }
      # weightedRF
      if(model == "non-linear"){
        mat_rf <- weightedRF_inference(counts, genes, tfs, nTrees = S,
                                       alpha = alpha, quiet = quiet,
                                       pwm_occurrence = prior_matrix,
                                       nCores = nCores,
                                       importance = "%IncMSE")
        
        mat_rf_perm <- weightedRF_inference(counts, genes, tfs, nTrees = S,
                                            alpha = alpha, tf_expression_permutation = TRUE,
                                            pwm_occurrence = prior_matrix,
                                            nCores = nCores,quiet = quiet,
                                            importance = "%IncMSE")
        
        mats[[paste0("RF_", as.character(alpha),  '_trueData_', rep)]] <- mat_rf
        mats[[paste0("RF_", as.character(alpha),  '_shuffled_', rep)]] <- mat_rf_perm
      }
    }
  }
  
  # thresholds the regulatory weights at certain densities to build GRNs
  edges <- list()
  lmses <- data.frame(row.names = genes)
  densities <- c(0.005, 0.01,0.05)
  for(name in names(mats)){ 
    for(density in densities){# exploring importance threshold stringency
      edges[[paste0(name, '_', density)]] <-
        weightedLASSO_network(mats[[name]], density = density, prior_matrix, 
                              genes, tfs, decreasing = TRUE)
      
      lmses[,name] <- mats[[name]]["mse",]
    }
  }
  
  
  # number of cores for evaluation and alpha opt computation
  # (lower as it requires more RAM)
  nCores_alphas_opt = max(round(nCores/2), 1)
  
  
  # applying DIOgene to compute the optimal alpha for each gene
  print("------------- Computing the optimal value of alpha")
  alphas_opt <-  mcsapply(genes, get_opt_alpha_per_gene, mats = mats, 
                          lmses = lmses, dev = "student",
                          return_alpha = T, 
                          mc.cores = nCores_alphas_opt, 
                          pval.adjust = "fdr")
  
  
  # running the appropriate model (weightedLASSO or weightedRF) on 
  # a target-gene specific basis:
  print("------------- Learning the final gene specific GRN")
  if(model == "linear"){
    diogene <- weightedLASSO_inference(counts = counts, genes = genes, tfs = tfs,
                                         alpha = alphas_opt, N = S, EN_param = 1,
                                         mda_type = "shuffle",
                                         quiet = quiet,
                                         pwm_occurrence = prior_matrix,
                                         nCores = nCores)
  }
  if(model == "non-linear"){
    diogene <- weightedRF_inference(counts, genes, tfs, nTrees = S,
                                   alpha = alphas_opt, quiet = quiet,
                                   pwm_occurrence = prior_matrix,
                                   nCores = nCores,
                                   importance = "%IncMSE")
  }
  
  # returns the results
  return(list(mats = mats, mse = lmses, edges = edges, 
              alphas_opt = alphas_opt, diogene = diogene))
  
}
