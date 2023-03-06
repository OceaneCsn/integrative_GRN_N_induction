# What is the precision and recall on the different validation datasets of the prior PWM network?
source("inference_functions/evaluateNetwork.R")
source("inference_functions/bRF.R")

load('rdata/pwm_prom_jaspar_dap.rdata')
load(file = "rdata/inference_input_N_response_varala.rdata")
load("rdata/pwm_occurrences_N_response_varala.rdata")

genes <- input_data$grouped_genes
known_tfs <- unique(pwm_prom$TF)
known_tfs <- intersect(known_tfs, colnames(pwm_occurrence))
tfs <- input_data$grouped_regressors

dim(pwm_occurrence)

# number of edges of the PWM network
sum(pwm_occurrence[,known_tfs])

links <- getLinkList(t(pwm_occurrence[,known_tfs]), reportMax = sum(pwm_occurrence[,known_tfs]))
colnames(links)[1:2] <- c('from', 'to')
head(links)

evaluate_network(links, input_genes = genes, input_tfs = known_tfs)[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("CHIPSeq"))[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("TARGET"))[c("tpr", "recall")]

evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("DAPSeq"))[c("tpr", "recall", "fn", "tp")]




########### on positive genes

# groups of genes depending on their behavior toward data integration
load("results/clusters_mse_bRF_100permutations.rdata")
positive_genes_rf <- names(clusters_rf[clusters_rf==2])

links <- getLinkList(t(pwm_occurrence[positive_genes_rf,known_tfs]), reportMax = sum(pwm_occurrence[positive_genes_rf,known_tfs]))
colnames(links)[1:2] <- c('from', 'to')
head(links)

evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence), validation = c("CHIPSeq"))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence), validation = c("TARGET"))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_rf, input_tfs = known_tfs, validation = c("DAPSeq"))[c("tpr", "recall")]






load("results/clusters_mse_lasso_100permutations.rdata")
positive_genes_lasso <- names(clusters_lasso[clusters_lasso==1])

links <- getLinkList(t(pwm_occurrence[positive_genes_lasso,known_tfs]), reportMax = sum(pwm_occurrence[positive_genes_lasso,known_tfs]))
colnames(links)[1:2] <- c('from', 'to')
head(links)

evaluate_network(links, input_genes = positive_genes_lasso, input_tfs = colnames(pwm_occurrence))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_lasso, input_tfs = colnames(pwm_occurrence), validation = c("CHIPSeq"))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_lasso, input_tfs = colnames(pwm_occurrence), validation = c("TARGET"))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_lasso, input_tfs = colnames(pwm_occurrence), validation = c("DAPSeq"))[c("tpr", "recall")]





val_conn$tpr
val_conn$recall

val_chip <- evaluate_network(links, input_genes = genes, input_tfs = colnames(pwm_occurrence))
val_conn$tpr
val_conn$recall


network <- graph_from_data_frame(links, directed = T)
edges <- as_long_data_frame(network)[c(4,5)]

bRF_network(pwm_occurrence)


valdap <- validated_edges[stringr::str_detect(validated_edges$type, "DAPSeq"),]
intersect(genes, unique(valdap$to))
intersect(tfs, unique(valdap$from))
studied_tfs <- unique(valdap$from)

both <- intersect(known_tfs, studied_tfs)

setdiff(known_tfs, studied_tfs)

tf <- "AT4G31800"

validation_fraction <- function(tf){
  return(length(setdiff(valdap[valdap$from==tf,"to"],
                names(which(pwm_occurrence[,tf] == 1))))
         / length(names(which(pwm_occurrence[,tf] == 1))))
}

sapply(both, FUN = validation_fraction)

hist(sapply(both, FUN = validation_fraction), breaks = 100)

nrow(merge(links, valdap, by = c('from', 'to')))/sum(links$from %in% studied_tfs)

# loaded new connecTF with amplification
# new_val <- validated_edges
# new_val_n <- new_val[new_val$from %in% tfs & new_val$to %in% genes,]
# 
# nrow(merge(valdap, new_val_n, by = c("from", "to")))
# validated_edges <- new_val_n
# 
# unique(new_val_n[stringr::str_detect(new_val_n$type, "DAPSeq"),]$from)
# unique(valdap[stringr::str_detect(valdap$type, "DAPSeq"),]$from)
# 
# 
# save(validated_edges, file = "rdata/connectf_N_responsive_genes_amp.rdata")


# dirty tests to see if the amplification changes a lot the precision of inferred grns


load("results/100_permutations_bRF_edges.rdata")
load("results/gene_specific_alphas_grns_max_div_0others.rdata")

evaluate_network(networks$`rf-max_div-0-all-trueData-5-0.005`, 
                 input_genes = genes, input_tfs = tfs, validation = c("DAPSeq"))[c("tpr", "recall")]


evaluate_network(networks$`rf-max_div-0-all-shuffled-3-0.005`, 
                 input_genes = genes, input_tfs = tfs, validation = c("DAPSeq"))[c("tpr", "recall")]


evaluate_network(edges$bRF_0_trueData_1_0.005, 
                 input_genes = genes, input_tfs = tfs, validation = c("DAPSeq"))[c("tpr", "recall")]


evaluate_network(edges$bRF_0_shuffled_4_0.005, 
                 input_genes = genes, input_tfs = tfs, validation = c("DAPSeq"))[c("tpr", "recall")]



evaluate_network(edges$bRF_1_trueData_1_0.005, 
                 input_genes = genes, input_tfs = tfs, validation = c("DAPSeq"))[c("tpr", "recall")]


evaluate_network(edges$bRF_1_shuffled_4_0.005, 
                 input_genes = genes, input_tfs = tfs, validation = c("DAPSeq"))[c("tpr", "recall")]
