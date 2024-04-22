source('inference_functions/weightedRF.R')
source('inference_functions/weightedLASSO.R')
source('inference_functions/evaluateNetwork.R')
source('inference_functions/data_integration_optimization.R')


DIOgene <- function(counts, genes, tfs, pwm_occurrence,
                    model = "non-linear",
                    nrep = 100,  nCores = 45,
                    nTrees = 2000,
                    quiet = T,
                    ALPHAS = seq(0,1, by = 0.1),
                    S = NULL){
  
  # number of bootstrapped models depends on the model type
  if(is.null(S)){
    if(model == "linear") S = 50
    if(model == "non-linear") S = 2000
  }
  
  
  # looping over alpha and replicates to build the MSE and EDI curves
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
                                             pwm_occurrence = pwm_occurrence,
                                             nCores = nCores)
        
        mat_lasso_perm <- weightedLASSO_inference(counts, genes, tfs,
                                                  alpha = alpha, N = S, 
                                                  tf_expression_permutation = TRUE,
                                                  mda_type = "shuffle",
                                                  quiet = quiet,
                                                  pwm_occurrence = pwm_occurrence,
                                                  nCores = nCores)
        
        mats[[paste0("LASSO_", as.character(alpha),  '_trueData_', rep)]] <- mat_lasso
        mats[[paste0("LASSO_", as.character(alpha),  '_shuffled_', rep)]] <- mat_lasso_perm
      }
      # weightedRF
      if(model == "non-linear"){
        mat_rf <- weightedRF_inference(counts, genes, tfs, nTrees = S,
                                       alpha = alpha, quiet = quiet,
                                       pwm_occurrence = pwm_occurrence,
                                       nCores = nCores,
                                       importance = "%IncMSE")
        
        mat_rf_perm <- weightedRF_inference(counts, genes, tfs, nTrees = S,
                                            alpha = alpha, tf_expression_permutation = TRUE,
                                            pwm_occurrence = pwm_occurrence,
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
        weightedLASSO_network(mats[[name]], density = density, pwm_occurrence, 
                              genes, tfs, decreasing = TRUE)
      
      lmses[,name] <- mats[[name]]["mse",]
    }
  }
  
  
  # number of cores for evaluation and alpha opt computation
  # (lower as it requires more RAM)
  nCores_alphas_opt = max(round(nCores/2), 1)
  
  
  # validation of GRN edges against DAP-Seq
  settings <- c("model", "alpha", "dataset", "rep", "density")
  val_dap <-
    evaluate_networks(
      edges,
      validation = c("DAPSeq"),
      nCores = nCores_alphas_opt,
      input_genes = genes,
      input_tfs = tfs
    )
  val_dap[, settings] <-
    str_split_fixed(val_dap$network_name, '_', length(settings))
  
  
  
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
                                         pwm_occurrence = pwm_occurrence,
                                         nCores = nCores)
  }
  if(model == "non-linear"){
    diogene <- weightedRF_inference(counts, genes, tfs, nTrees = S,
                                   alpha = alphas_opt, quiet = quiet,
                                   pwm_occurrence = pwm_occurrence,
                                   nCores = nCores,
                                   importance = "%IncMSE")
  }
  
  # returns the results
  return(list(mats = mats, mse = lmses, edges = edges, val_dap = val_dap,
              alphas_opt = alphas_opt, diogene = diogene))
  
}
