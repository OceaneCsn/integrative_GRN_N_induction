load('rdata/connectf_N_responsive_genes.rdata')
library(tidyverse)
# library(ggh4x)
library(igraph)
library(pROC)
library(PRROC)

evaluate_grn <- function(mat, plot = F, pAUC_specificity=0.9, 
                         validation = "DAPSeq"){
  
  # validation edges restricted to genes of interest
  val <- validated_edges %>%
    filter(type == validation) %>%
    filter(from %in% rownames(mat),
           to %in% colnames(mat))
  
  # predicted complete grn is joined with validation edges
  val_data <-mat %>%
    reshape2::melt() %>%
    rename(from = Var1, to = Var2, importance = value) %>%
    filter(from != "mse") %>%
    filter(from %in% val$from) %>%
    full_join(val, by = c("from", "to")) %>%
    mutate(gold = ifelse(is.na(type), 0, 1))
  
  # evaluation metrics are computed:
  ##### AUC
  partial_auc <- val_data %>%
    roc(gold, importance, plot = F, direction = "<") %>%
    auc(partial.auc = c(1, pAUC_specificity), 
        partial.auc.correct=T, 
        partial.auc.focus = "specificity") %>%
    as.numeric()
  
  auc =  val_data %>%
    roc(gold, importance, plot = F, 
        direction = "<", plot.roc = F) %>%
    ci.auc() %>%
    as.numeric()
  
  curve <- val_data %>%
    roc(gold, importance, plot = F, direction = "<") %>%
    coords(transpose =F, ret = "all")  %>%
    top_frac(1,threshold)  %>%
    ggplot(aes(y=sensitivity, x = 1-specificity))+
    geom_point()+geom_line()+
    geom_abline(slope = 1, intercept = 0)+
    ggtitle(paste("AUC =", round(auc, 3), "- pAUC =", round(partial_auc, 3)))+
    geom_vline(xintercept = 1-pAUC_specificity)
  
  ####### AUPR
  pr <- pr.curve(scores.class0 = val_data$importance,
                   weights.class0 = val_data$gold, rand.compute = T)
  aupr <- pr$auc.integral
  aupr_rand <- pr$rand$auc.integral
  
  if(plot)
    print(curve)
  
  return(list(auc = auc[2], auc.lower = auc[1], auc.higher = auc[3],
              aupr = aupr, aupr_rand = aupr_rand,
              partial_auc = partial_auc, plot = curve))
}


evaluate_genes <- function(mat, plot = F, pAUC_specificity = 0.9, nCores = 34,
                           validation = "DAPSeq"){
  
  # validation edges restricted to genes of interest
  val <- validated_edges %>%
    filter(type == validation) %>%
    filter(from %in% rownames(mat),
           to %in% colnames(mat))
  
  # predicted complete grn is joined with validation edges
  get_auc_gene <- function(gene){
    val_data <-mat %>%
      reshape2::melt() %>%
      rename(from = Var1, to = Var2, importance = value) %>%
      filter(to == gene) %>%
      filter(from %in% val$from) %>%
      left_join(val, by = c("from", "to")) %>%
      mutate(gold = ifelse(is.na(type), 0, 1))
    
    if(sum(val_data$gold)==0){
      return(NA)
    }
    if(sum(val_data$gold)==nrow(val_data)){
      return(NA)
    }
    else{
      auc =  val_data %>%
        roc(gold, importance, plot = F, direction = "<", auc = F) %>%
        auc() %>%
        as.numeric()
      
      partial_auc = val_data %>%
        roc(gold, importance, plot = F, direction = "<", auc = F) %>%
        auc(partial.auc = c(1, pAUC_specificity), 
            #partial.auc.correct=T, allow.invalid.partial.auc.correct=T,
            partial.auc.focus = "specificity") %>%
        as.numeric()
      
      pr <- pr.curve(scores.class0 = val_data$importance,
                     weights.class0 = val_data$gold, rand.compute = T)
      aupr <- pr$auc.integral
      aupr_rand <- pr$rand$auc.integral
      
      
      return(c(auc, partial_auc, aupr, aupr_rand))
    }
  }
  aucs <- mcsapply(colnames(mat), FUN = get_auc_gene, mc.cores=nCores)
  
  return(t(data.frame(aucs)) %>%
           data.frame() %>%
           rename(auc = X1, p_auc = X2, aupr = X3, aupr_rand = X4))
}

#' Evaluates an inferred network against validated regulatory interactions
#'
#' @param net dataframe with a column 'from' (regulators) and a column 'to'
#' (targets), representing the inferred network of edges to evaluate
#' @param input_genes vector of gene AGIs used for network inference as input
#' (useful to accurately compute recall, so that genes in input but not in
#' predicted edges are counted as false negatives. If not given, will only 
#' count as false negatives genes missed in the predicted network)
#' @param input_tfs vector of TFs AGIs used for network inference as input
#' (useful to accurately compute recall, so that genes in input but not in
#' predicted edges are counted as false negatives. If not given, will only 
#' count as false negatives genes missed in the predicted network)
#' @param validation type of edge in the validation database that
#' should be considered to defined a true/supported prediction in the
#' evaluation process.
#' The validation type must be a vector of one or more
#' of the following values : CHIPSeq, DAPSeq, Litterature, TARGET
#' @return a list containing true positives, true positive rate (precision),
#' false positives, false positive rate, the false positives,
#' and the recall value. It also returns the input network dataframe
#' with an additional column to characterize how is the edge supported
#' by known interactions
evaluate_network <-
  function(net, input_genes = NULL, input_tfs = NULL,
           validation = c("CHIPSeq", "DAPSeq", "TARGET"),
           subset_validated_edges = NULL, tf_validation_subset = NULL) {
    
    if (nrow(net) == 0) {
      warning("Empty dataframe of edges \n")
      return(list(
        tp = NA,
        fp = NA,
        tpr = NA,
        fpr = NA,
        fn = NA,
        recall = NA
      ))
    }
    
    if (!is.null(subset_validated_edges))
      validated_edges <- subset_validated_edges
    
    if (!is.null(tf_validation_subset))
      validated_edges <- validated_edges[validated_edges$from %in% tf_validation_subset, ]
    
    # ----------------------------- CHECKS ---------------------------------- #
    
    if (any(!validation %in% c("CHIPSeq", "DAPSeq", "Litterature", "TARGET")))
      stop(
        "The validation type must be a vector of one or more of the following values :
           CHIPSeq, DAPSeq, Litterature, TARGET"
      )
    
    if (!("from" %in% colnames(net) &
          "to" %in% colnames(net)))
      stop("The network dataframe should have two columns named 'from' and 'to")
    
    from_DIANE <- FALSE
    
    # regex check for Arabidopsis AGIs
    matched <-
      sum(stringr::str_detect(c(net$from, net$to), 
                              pattern = "^AT[[:alnum:]]G[[:digit:]]{5}"))
    
    if (matched != 2 * nrow(net)) {
      if (matched > 0)
        stop("Some of the gene IDs do not match the expected regex for Arabidopsis AGIs")
      else
        stop("None of the gene IDs match the expected regex Arabidopsis AGIs")
    }
    
    # all genes studied in the network inference
    if(is.null(input_genes))
      input_genes <- unique(net$to)
    
    
    # all genes studied in the network inference
    if(is.null(input_tfs))
      input_tfs <- unique(net$from)
    
    # ------------------------------------------------------------------------ #
    
    
    # validation from specific type of validation
    validated_edges_specific <-
      validated_edges[validated_edges$type %in% validation,]
    
    # aggregate validation edges to have unique validated pairs
    validated_edges_specific_unique <-
      aggregate(. ~ from + to,
                data = validated_edges_specific,
                FUN = paste0,
                collapse = '+')
    
    # restricting validation edges to TFs and genes present in the input data
    # (useful to compute false negatives, missed edges)
    validated_edges_specific_unique <- validated_edges_specific_unique[
      validated_edges_specific_unique$from %in% input_tfs & 
        validated_edges_specific_unique$to %in% input_genes,]
    

        # TFs for which we possess validation data
    studied_tfs <- unique(validated_edges_specific$from)
    n_studied_interactions <- sum(net$from %in% studied_tfs)
    
    if (length(studied_tfs) == 0) {
      warning("No regulator present in your network was studied in the database \n")
      return(list(
        tp = NA,
        fp = NA,
        tpr = NA,
        fpr = NA,
        fn = NA,
        recall = NA
      ))
    }
    
    
    # validated edges
    val <-
      merge(net, validated_edges_specific,
            by = c('from', 'to'))
    
    if (nrow(val) == 0) {
      warning(
        "No predicted edge was found in the validation databse...Coffee to cheer you up? \n"
      )
      edges <- net
      edges$type = "Not supported"
      return(list(
        tp = 0,
        fp = n_studied_interactions,
        tpr = 0,
        fpr = 1,
        fn = nrow(validated_edges_specific_unique),
        recall = 0,
        edges = edges
      ))
    }
    
    val_unique <-
      aggregate(. ~ from + to,
                data = val,
                FUN = paste0,
                collapse = '+')
    
    # true positives
    tp <- nrow(val_unique)
    
    # true positive rate (misnommer) : PRECISION
    if (nrow(val) == 0) {
      tpr = 0
    }
    else{
      tpr <- tp / n_studied_interactions
    }
    
    # false positives
    fp <- n_studied_interactions - tp
    
    # false positive rate
    # negs = n_possible_val_interactions - val_interactions
    
    # number of TFs present in validation
    n_val_tfs <- length(unique(validated_edges_specific_unique$from))
    
    # number of genes present in validation
    n_val_genes <- length(unique(c(validated_edges_specific_unique$from, validated_edges_specific_unique$to)))
    
    # number of negative interactions
    negs <- n_val_tfs * (n_val_genes -1) - nrow(validated_edges_specific_unique)
    
    fpr <- fp / negs
    
    # false negatives
    fn <- nrow(validated_edges_specific_unique) - tp
    
    # recall (real TPR)
    recall <- tp / (tp+fn)
    
    
    # dataframe of edges with validation information
    edges <- merge(net, val_unique,
                   all = TRUE, by = c('from', 'to'))
    edges$type <- ifelse(edges$from %in% studied_tfs &
                           is.na(edges$type),
                         "Not supported",
                         edges$type)
    
    edges$type <-
      ifelse(is.na(edges$type), "TF not studied", edges$type)
    
    # sorting type of edge
    reorder_type <- function(type){
      types <- unlist(strsplit(type, '\\+'))
      return(paste(types[order(types)], collapse = '+'))
    }
    
    edges$type <- sapply(edges$type, reorder_type)
    edges <- edges[order(edges$type),] 
    
    results <- list(
      tp = tp,
      fp = fp,
      tpr = tpr,
      precision = tpr,
      fpr = fpr,
      fn = fn,
      recall = recall,
      edges = edges
    )
    return(results)
  }



# runs the evaluate network function in parallel for a list of networks
# returns the results as a dataframe with precision and recall for each network
# of the list
evaluate_networks <- function(edges_list, validation = c("TARGET", "CHIPSeq", "DAPSeq"), nCores = 32,
                              input_genes = NULL, input_tfs = NULL, tf_validation_subset = NULL){
  
  registerDoParallel(cores = nCores)
  message(paste("\nParallel network validation with AranetBench Using", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(
    result.validation <-
      doRNG::"%dorng%"(foreach::foreach(network_name = names(edges_list), .combine = rbind),
                       {
                         val <- evaluate_network(edges_list[[network_name]], 
                                                 validation = validation, 
                                                 input_genes = input_genes, 
                                                 input_tfs = input_tfs, 
                                                 tf_validation_subset = tf_validation_subset)
                         
                         data.frame("network_name" = c(network_name), 
                                    "precision" = c(val$tpr), 
                                    "recall" = c(val$recall))
                         
                       })
  )
  toc()
  attr(result.validation, "rng") <- NULL # It contains the whole sequence of RNG seeds
  return(result.validation)
}


# return a precision/recall curve for a range of different importance threshold (densities)
# this is useful to draw precision and recall curves
evaluate_fully_connected <- function(mat, pwm_occurrence = NULL, validation = c("TARGET", "CHIPSeq", "DAPSeq"), nCores = 32,
                              input_genes = NULL, input_tfs = NULL, tf_validation_subset = NULL){
  
  # reads the input matrix
  imp <-  reshape2::melt(mat) %>%
    rename(from = Var1, to = Var2, importance = value) %>%
    mutate(importance = ifelse(importance<0, 0, importance))
  
  if(!is.null(pwm_occurrence)){
    pwm_imputed <- pwm_occurrence
    pwm_imputed[is.na(pwm_imputed)] <- 0.5
    imp$pwm <- pwm_imputed[cbind(imp$to, imp$from)]
  }
  
  
  # probs <- seq(0, 1, length.out = N)
  # densities <- quantile(imp$importance, probs)
  # or
  # densities <- seq(min(imp$importance), max(imp$importance), length.out = N)
  
  top_edges <- c(1, 0.5,0.2,0.1,0.075, 0.05,0.025,0.01,0.0075,
                 0.005, 0.0025,0.002,0.0015, 0.001)
  
  edges_list <- list()
  for(t in top_edges){
    edges_list[[as.character(t)]] <- imp %>%
      top_frac(t, wt = importance)
      # filter(importance > 0)
  }
  
  registerDoParallel(cores = nCores)
  message(paste("\nParallel network validation using", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(
    result.validation <-
      doRNG::"%dorng%"(foreach::foreach(network_name = names(edges_list), .combine = rbind),
                       {
                         val <- evaluate_network(edges_list[[network_name]], 
                                                 validation = validation, 
                                                 input_genes = input_genes, 
                                                 input_tfs = input_tfs, 
                                                 tf_validation_subset = tf_validation_subset)
                         
                         # returns all the info
                         if(!is.null(pwm_occurrence)){
                           data.frame("network_name" = c(network_name), 
                                      "precision" = c(val$tpr), 
                                      "recall" = c(val$recall),
                                      "pwm_supprt" = mean(edges_list[[network_name]]$pwm))
                         }
                         else
                           data.frame("network_name" = c(network_name), 
                                      "precision" = c(val$tpr), 
                                      "recall" = c(val$recall))
                         
                       })
  )
  toc()
  attr(result.validation, "rng") <- NULL # It contains the whole sequence of RNG seeds
  return(result.validation)
}



# computes the degree of input networks (in a list)
get_networks_degrees <- function(nets, nCores = 30){
  
  n_genes <- read.csv("data/Ngenes.csv", header = T, sep = ';')%>%
    rename(genes = AGI)
  
  
  registerDoParallel(cores = nCores)
  message(paste("\nParallel network validation using", foreach::getDoParWorkers(), "cores."))
  "%dopar%" <- foreach::"%dopar%"
  tic()
  suppressPackageStartupMessages(
    result.validation <-
      doRNG::"%dorng%"(foreach::foreach(network_name = names(nets), .combine = rbind),
                       {
                         graph <- graph_from_data_frame(d = nets[[network_name]], directed = T)
                         data.frame(degree = igraph::degree(graph, mode = "total"), 
                                    genes = names(igraph::degree(graph, mode = "total"))) %>%
                           
                           mutate(network_name = c(network_name))
                         
                       })
  )
  toc()
  attr(result.validation, "rng") <- NULL # It contains the whole sequence of RNG seeds
  return(result.validation)
}


draw_validation <- function(validation, densities = c(0.005,0.01,0.05), 
                            precision_only = F, recall_only = F){
  data_val <- validation %>%
    filter(density %in% densities)%>%
    group_by(model, alpha, dataset, density) %>%
    mutate(mean_precision = mean(precision, na.rm = T),
           sd_precision = sd(precision, na.rm = T),
           mean_recall = mean(recall, na.rm = T),
           sd_recall = sd(recall, na.rm = T),
           density = paste("D =", density),
           model = str_replace(model, "bRF", "weightedRF"),
           model = str_replace(model, "LASSO", "weightedLASSO")) %>%
    dplyr::select(-network_name) %>%
    mutate(alpha = as.numeric(alpha))
  
  precision <- data_val %>%
    ggplot(aes(
      x = as.numeric(alpha),
      y = precision,
      color = dataset,
      fill = dataset
    ))+
    ggh4x::facet_nested_wrap(vars(density), ncol = 8, nest_line = T) + 
    geom_ribbon(aes(ymin = mean_precision - sd_precision , 
                    ymax = mean_precision + sd_precision  ), 
                alpha = .4) + theme_pubr(legend = "top")+
    geom_point(alpha = 0.1) +
    xlab(expression(alpha)) + ylab("Precision") +
    ggtitle(paste("Precision against DAP-Seq")) +
    geom_line(aes(group = dataset, y = mean_precision)) +
    theme(
      strip.background = element_blank(),
      axis.title.x = element_text(size = 22),
      title = element_text(size = 16),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.position = 'top'
    )+scale_color_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))+
    scale_fill_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData"))) +
    geom_hline(color = "#C55A11", yintercept = 0.331042, size= 1.5, show.legend = T)
  
  recall <- data_val %>%
    ggplot(aes(
      x = as.numeric(alpha),
      y = recall,
      color = dataset,
      fill = dataset
    ))+
    ggh4x::facet_nested_wrap(vars(density), ncol = 8, nest_line = T, scales = "free") + 
    geom_ribbon(aes(ymin = mean_recall - sd_recall , 
                    ymax = mean_recall + sd_recall  ), 
                alpha = .4)  +theme_pubr(legend = "top")+
    geom_point(alpha = 0.1) + xlab("alpha") +
    xlab(expression(alpha)) + ylab("Recall") +
    geom_line(aes(group = dataset, y = mean_recall)) +
    ggtitle(paste("Recall against DAPSeq")) +
    theme(
      strip.background = element_blank(),
      axis.title.x = element_text(size = 22),
      title = element_text(size = 16),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.position = 'top'
    )+scale_color_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))+
    scale_fill_manual(values = setNames(c("grey", "#70AD47"), c("shuffled", "trueData")))
  
  if(precision_only)
    return(precision)
  if(recall_only)
    return(recall)
  if(!precision_only & !precision_only)
    return(precision/recall)
}




# Plots the ranking of TFs based on their out degree and highlights
# the known nitrate regulators
get_hubs <- function(nets, metric = "degree", unique_regs = c()){
  g <- setNames(rep(0, length(genes)), genes)
  for(net in names(nets)){
    if(metric == "degree")
      g[names(table(nets[[net]]$from))] <- g[names(table(nets[[net]]$from))] + table(nets[[net]]$from)/length(nets)
    if(metric == "betweenness")
      g[names(table(nets[[net]]$from))]<- g[names(table(nets[[net]]$from))]+
        igraph::betweenness(igraph::graph_from_data_frame(nets[[net]], directed = T),
                            v = names(table(nets[[net]]$from)), directed = T)
  }
  df <- data.frame(genes = genes, average_degree = g[genes], 
                   label = annot[match(genes, rownames(annot)), "label"])
  
  df$regulator <- ifelse(df$label %in% c("DIV1", "TGA1", "TGA4", "HHO2", "HHO3", "HRS1", "BT1", "BT2"), 
                         "Nitrate regulator",
                         ifelse(df$label %in%  c("VRN1", "CRF4"), "Candidate nitrate\nregulator",
                                ifelse(df$label %in% unique_regs,
                                       "Retreived by gene\nspecific data integration\n optimisation", "No knowlege")))
  
  
  plt <- df %>%
    slice_max(average_degree, n = 25) %>% 
    ggplot(aes(x=average_degree, label = label, color = regulator,
               y=reorder(label, average_degree)))+ 
    geom_segment(x=0, color = "grey",
                 aes(xend=average_degree, 
                     yend = reorder(label, average_degree)))+
    geom_point()+
    geom_text(nudge_x = 7, show.legend = F)+ theme_bw() + 
    ylab("") + xlab("Average out-degree")+theme(axis.text.y = element_blank())+
    scale_color_manual(name = "", values = setNames(c("#C55A11", "orange", "#4670CD", "black"), 
                                                    c("Nitrate regulator", "Candidate nitrate\nregulator",
                                                      "Retreived by gene\nspecific data integration\n optimisation", "No knowledge")))+
    guides(color = guide_legend(nrow=4,byrow=TRUE))
  
  top20<- df %>%
    slice_max(average_degree, n = 25) %>%
    select(label)
  
  return(list(ranking = plt, top25 = top20$label))
}



# computes other topological metrics on a list of GRNs
get_transitivity <- function(nets, settings){
  res <- data.frame(trans = unlist(lapply(names(nets), 
                                          function(net){return(igraph::transitivity(igraph::graph_from_data_frame(nets[[net]], directed = T), 
                                                                                    type = "global"))})),
                    length = unlist(lapply(names(nets), 
                                           function(net){return(igraph::mean_distance(igraph::graph_from_data_frame(nets[[net]], directed = T)))})),
                    diam = unlist(lapply(names(nets), 
                                         function(net){return(igraph::diameter(igraph::graph_from_data_frame(nets[[net]], directed = T)))})))
  
  res$settings <- settings
  return(res)
}




# draw density of importances values of edges present in the gold standard,
# and compares two models (m1, m2), given as importance matrices.
draw_importances_comparison <- function(m1, m2, saved, lost, same, edge_type = "DAPSeq"){
  # m1 : min mse
  # m2 : dev criterion
  m2 <-  reshape2::melt(m2) %>%
    rename(from = Var1, to = Var2, baseline = value) 
  
  min <- 0.75*min(m2$baseline[m2$baseline>0])
  
  known_tfs <- unique(validated_edges[validated_edges$type=="DAPSeq",]$from)
  
  importances <- reshape2::melt(m1) %>%
    rename(from = Var1, to = Var2, min_MSE = value) %>%
    left_join(m2, by = c("from", "to")) %>%
    left_join(validated_edges, by = c("from", "to")) %>%
    left_join(prior, by = c("from", "to")) %>%
    mutate(importance_change = baseline - min_MSE)  %>%
    mutate(type = ifelse(is.na(type), 
                         ifelse(from %in% known_tfs, "unsupported by DAP-Seq", "unknown"), 
                         type)) %>%
    mutate(target_type = ifelse(to %in% saved, "2. saved",
                                ifelse(to %in% lost, "3. lost", 
                                       ifelse(to %in% same, "1. both", "4. none")))) %>%
    filter(type == edge_type) %>%
    mutate(pi = paste("pi =", pi, "-")) %>%
    mutate(baseline = ifelse(baseline <0, 0, baseline),
           min_MSE=ifelse(min_MSE<0, 0, min_MSE))
  
  baseline_ranking <- importances %>%
    top_n(1432, wt = baseline)
  thr_dev <- log(baseline_ranking[nrow(baseline_ranking),]$baseline+min)
  
  min_ranking <- importances %>%
    top_n(1432, wt = min_MSE)
  thr_min <- log(min_ranking[nrow(min_ranking),]$min_MSE+min)
  
  return(importances %>%
           # filter(baseline > 0) %>%
           ggplot(aes(x = log(baseline+min),
                      y = interaction(target_type), fill = interaction(target_type))) +
           stat_density_ridges(scale = 2, geom = "density_ridges_gradient", calc_ecdf = TRUE,
                               quantiles = 2, quantile_lines = TRUE, alpha = 0.5, rel_min_height = 0.01) + 
           scale_fill_viridis_d(name = "Target gene")+
           geom_vline(xintercept = thr_dev, color = "darkorange", linewidth = 1.5)
         +
           importances%>%
           # filter(min_MSE > 0) %>%
           ggplot(aes(x = log(min_MSE+min),
                      y = interaction(target_type), fill = interaction(target_type))) +
           stat_density_ridges(scale = 2, geom = "density_ridges_gradient", calc_ecdf = TRUE,
                               quantiles = 2, quantile_lines = TRUE, alpha = 0.5, rel_min_height = 0.01) +
           scale_fill_viridis_d(name = "Target gene")+
           geom_vline(xintercept = thr_min, color = "darkorange", linewidth = 1.5)+
           plot_layout(guides = "collect") & theme(legend.position = "top"))
}





N_go <- c("GO:0006807", "GO:0051171", "GO:0051173", "GO:0022622", "GO:0034641", "GO:0044271")

convert_from_agi <- function(ids){
  x <- org.At.tair.db::org.At.tairENTREZID
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  return(unlist(xx[as.vector(ids)]))
}


get_gsea <- function(nets, settings){
  g <- setNames(rep(0, length(genes)), genes)
  for(net in names(nets)){
    
    # g[names(table(nets[[net]]$to))]<-g[names(table(nets[[net]]$to))]+
    #   igraph::betweenness(igraph::graph_from_data_frame(nets[[net]], directed = T),
    #                       v = names(table(nets[[net]]$to)), directed = F)
    
    # g[names(table(nets[[net]]$to))]<-g[names(table(nets[[net]]$to))]+
    #   igraph::transitivity(igraph::graph_from_data_frame(nets[[net]], directed = T),
    #                       vids = names(table(nets[[net]]$to)), type = "local")
    #g[names(table(nets[[net]]$from))] <- g[names(table(nets[[net]]$from))] + table(nets[[net]]$from)/length(nets)
    
    # g[genes] <- g[genes] + (table(nets[[net]]$from)[genes]+table(nets[[net]]$to)[genes])/length(nets)
    ge <- intersect(genes, unique(c(nets[[net]]$from, nets[[net]]$to)))
    g[ge]<-g[ge]+
      igraph::degree(igraph::graph_from_data_frame(nets[[net]], directed = T),
                     v = ge, mode = "total")
  }
  df <- data.frame(genes = genes, average_degree = g[genes],
                   label = annot[match(genes, rownames(annot)), "label"])
  
  df <- df[df$average_degree>0,]
  degrees <- sort(setNames(df$average_degree,convert_from_agi(df$genes)),
                  decreasing = T)
  
  degree <- degrees[!is.na(names(degrees))]
  
  gse <- gseGO(geneList=degrees,
               ont ="BP",
               maxGSSize = length(degrees)-1, nPerm =1000,
               pvalueCutoff = 1,
               verbose = TRUE, scoreType = "pos",
               OrgDb = org.At.tair.db::org.At.tair.db,
               pAdjustMethod = "BH")
  
  res <- gse@result[gse@result$ID %in% N_go,c("Description", "enrichmentScore", "p.adjust")]
  res$settings = settings
  return(res)
}

