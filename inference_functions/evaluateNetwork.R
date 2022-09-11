load('rdata/connectf_N_responsive_genes.rdata')

#' Evaluates an inferred network against validated regulatory interactions
#'
#' @param net dataframe with a column 'from' (regulators) and a column 'to'
#' (targets), representing the inferred network of edges to evaluate
#' @param input_genes vector of gene AGIs used for network inference as input
#' (useful to accurately compute recall, so that genes in input but not in
#' predicted edges are counted as false negatives. If not given, will only 
#' count as false negatives genes missed in the predicted network)
#' @param validation type of edge in the validation database that
#' should be considered to defined a true/supported prediction in the
#' evaluation process.
#' The validation type must be a vector of one or more
#' of the following values : CHIPSeq, DAPSeq, Litterature, TARGET
#' @param subset_validated_edges potential dataframe of validated edges
#' (that should however contain all the nodes in the inferred network).
#' Made to speed up computation in ROC and test again random.
#'
#' @return a list containing true positives, true positive rate (precision),
#' false positives, false positive rate, the false positives,
#' and the recall value. It also returns the input network dataframe
#' with an additional column to characterize how is the edge supported
#' by known interactions
#'
#' @export
#' @examples 
#' #' data("abiotic_stress_Heat_genes_net")
#' set.seed(999)
#' results <- evaluate_network(abiotic_stress_Heat_genes_net)
#' results[c("tp", "fp", "tpr", "fpr", "fn", "recall")]
evaluate_network <-
  function(net, input_genes = NULL, input_tfs = NULL,
           validation = c("CHIPSeq", "DAPSeq", "TARGET"),
           subset_validated_edges = NULL, tf_validation_subset = NULL) {
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
    
    # restricting validation edges to TFs and genes present in predicted network
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
    
    # true positive rate
    if (nrow(val) == 0) {
      tpr = 0
    }
    else{
      tpr <- tp / n_studied_interactions
    }
    
    # false positives
    fp <- n_studied_interactions - tp
    
    # false positive rate
    fpr <- fp / n_studied_interactions
    
    # false negatives
    fn <- nrow(validated_edges_specific_unique) - tp
    
    # recall
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
      fpr = fpr,
      fn = fn,
      recall = recall,
      edges = edges
    )
    return(results)
  }




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
