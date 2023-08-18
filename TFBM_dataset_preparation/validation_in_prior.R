# What is the precision and recall on the different validation datasets of the prior PWM network?
source("inference_functions/evaluateNetwork.R")

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

evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("DAPSeq"))[c("tpr", "recall", "fn", "tp")]