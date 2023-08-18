########## PWM occurrences in N response dataset
load(file = "rdata/inference_input_N_response_varala.rdata")
load('rdata/pwm_prom_jaspar_dap.rdata')

genes <- input_data$grouped_genes
tfs <- input_data$grouped_regressors
counts <- input_data$counts

pwm_prom_ <- pwm_prom[pwm_prom$target %in% genes | 
                        pwm_prom$TF %in% tfs, ]
pwm_prom_ <- pwm_prom_[pwm_prom_$pval < 1e-4,]
known_tfs <- unique(pwm_prom_$TF)

# construction de la matrice de présence du score
pwm_occurrence <- matrix(NA, nrow = length(genes),
                         ncol = length(tfs),
                         dimnames = list(genes, tfs))

for(tf in tfs){
  print(tf)
  if(length(intersect(tf, known_tfs)) >=1){
    pwm_prom_tf <- pwm_prom_[pwm_prom_$TF %in% tf,]
    for(gene in genes){
      pwm_occurrence[gene, tf] <-
        nrow(pwm_prom_tf[pwm_prom_tf$target %in% gene,])/(length(gene)*length(tf))
    }
  }
}

save(pwm_occurrence, file = "rdata/pwm_occurrences_N_response_varala.rdata")



############## testing the prom+- 1000kb model

########## PWM occurrences in N response dataset
load(file = "rdata/inference_input_N_response_varala.rdata")
load('../PWM/pwm_prom_extended_jaspar_dap.rdata')

genes <- input_data$grouped_genes
tfs <- input_data$grouped_regressors
counts <- input_data$counts

pwm_prom_ <- pwm_prom[pwm_prom$target %in% genes | 
                        pwm_prom$TF %in% tfs, ]
pwm_prom_ <- pwm_prom_[pwm_prom_$pval < 1e-4,]
known_tfs <- unique(pwm_prom_$TF)

# construction de la matrice de présence du score
pwm_occurrence <- matrix(NA, nrow = length(genes),
                         ncol = length(tfs),
                         dimnames = list(genes, tfs))

for(tf in tfs){
  print(tf)
  if(length(intersect(tf, known_tfs)) >=1){
    pwm_prom_tf <- pwm_prom_[pwm_prom_$TF %in% tf,]
    for(gene in genes){
      pwm_occurrence[gene, tf] <-
        nrow(pwm_prom_tf[pwm_prom_tf$target %in% gene,])/(length(gene)*length(tf))
    }
  }
}

save(pwm_occurrence, file = "rdata/pwm_occurrences_N_response_varala_prom_extended.rdata")

sum(pwm_occurrence, na.rm = T)





########## PWM occurrences in N response dataset
load(file = "rdata/inference_input_N_response_varala.rdata")
load('../PWM/pwm_jaspar_dapseq.rdata')

genes <- input_data$grouped_genes
tfs <- input_data$grouped_regressors
counts <- input_data$counts

pwm_prom_ <- pwm[pwm$target %in% genes | 
                   pwm$TF %in% tfs, ]
pwm_prom_ <- pwm_prom_[pwm_prom_$pval < 1e-4,]
known_tfs <- unique(pwm_prom_$TF)

# construction de la matrice de présence du score
pwm_occurrence <- matrix(NA, nrow = length(genes),
                         ncol = length(tfs),
                         dimnames = list(genes, tfs))

for(tf in tfs){
  print(tf)
  if(length(intersect(tf, known_tfs)) >=1){
    pwm_prom_tf <- pwm_prom_[pwm_prom_$TF %in% tf,]
    for(gene in genes){
      pwm_occurrence[gene, tf] <-
        nrow(pwm_prom_tf[pwm_prom_tf$target %in% gene,])/(length(gene)*length(tf))
    }
  }
}

save(pwm_occurrence, file = "rdata/pwm_occurrences_N_response_varala_prom_introns.rdata")

sum(pwm_occurrence, na.rm = T)
