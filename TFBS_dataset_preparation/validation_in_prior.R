# What is the precision and recall on the different validation datasets of the prior PWM network?
source("inference_functions/evaluateNetwork.R")
source("inference_functions/bRF.R")

load('rdata/pwm_prom_jaspar_dap.rdata')
load(file = "rdata/inference_input_N_response_varala.rdata")
load("rdata/pwm_occurrences_N_response_varala.rdata")

genes <- input_data$grouped_genes
known_tfs <- unique(pwm_prom$TF)
known_tfs <- intersect(known_tfs, colnames(pwm_occurrence))


dim(pwm_occurrence)

# number of edges of the PWM network
sum(pwm_occurrence[,known_tfs])

links <- getLinkList(pwm_occurrence[,known_tfs], reportMax = sum(pwm_occurrence[,known_tfs]))
colnames(links)[1:2] <- c('from', 'to')
head(links)

evaluate_network(links, input_genes = genes, input_tfs = colnames(pwm_occurrence))[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = colnames(pwm_occurrence), validation = c("CHIPSeq"))[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = colnames(pwm_occurrence), validation = c("TARGET"))[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = colnames(pwm_occurrence), validation = c("DAPSeq"))[c("tpr", "recall")]













########### on positive genes

# groups of genes depending on their behavior toward data integration
load("results/clusters_mse_bRF_100permutations.rdata")
positive_genes_rf <- names(clusters_rf[clusters_rf==2])

links <- getLinkList(pwm_occurrence[positive_genes_rf,known_tfs], reportMax = sum(pwm_occurrence[positive_genes_rf,known_tfs]))
colnames(links)[1:2] <- c('from', 'to')
head(links)

evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence), validation = c("CHIPSeq"))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence), validation = c("TARGET"))[c("tpr", "recall")]
evaluate_network(links, input_genes = positive_genes_rf, input_tfs = colnames(pwm_occurrence), validation = c("DAPSeq"))[c("tpr", "recall")]






load("results/clusters_mse_lasso_100permutations.rdata")
positive_genes_lasso <- names(clusters_lasso[clusters_lasso==1])

links <- getLinkList(pwm_occurrence[positive_genes_lasso,known_tfs], reportMax = sum(pwm_occurrence[positive_genes_lasso,known_tfs]))
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