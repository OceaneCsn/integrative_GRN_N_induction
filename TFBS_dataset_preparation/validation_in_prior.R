# What is the precision and recall on the different validation datasets of the prior PWM network?
source("inference_functions/evaluateNetwork.R")
source("inference_functions/bRF.R")

load('rdata/pwm_prom_jaspar_dap.rdata')
load(file = "rdata/inference_input_N_response_varala.rdata")
load("rdata/pwm_occurrences_N_response_varala.rdata")

genes <- input_data$grouped_genes
tfs <- input_data$grouped_regressors
known_tfs <- unique(pwm_prom$TF)
known_tfs <- intersect(known_tfs, colnames(pwm_occurrence))


dim(pwm_occurrence)

# number of edges of the PWM network
sum(pwm_occurrence[,known_tfs])

links <- getLinkList(t(pwm_occurrence[,known_tfs]), reportMax = sum(pwm_occurrence[,known_tfs]))
colnames(links)[1:2] <- c('from', 'to')
head(links)

evaluate_network(links, input_genes = genes, input_tfs = known_tfs)[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("CHIPSeq"))[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("TARGET"))[c("tpr", "recall")]
evaluate_network(links, input_genes = genes, input_tfs = known_tfs, validation = c("DAPSeq"))[c("tpr", "recall")]



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


valdap[val_dap$from %in% known_tfs, ]
