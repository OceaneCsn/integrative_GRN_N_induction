load(file = "rdata/inference_input_N_response_varala.rdata")
load('rdata/pwm_prom_jaspar_dap.rdata')

library(tidyverse)
library(ggpubr)

genes <- input_data$grouped_genes
tfs <- input_data$grouped_regressors
counts <- input_data$counts

# taking only interactions concerning Nitrate responsive genes
pwm_prom_ <- pwm_prom[pwm_prom$target %in% genes & 
                        pwm_prom$TF %in% tfs, ]

known_tfs <- intersect(pwm_prom$TF, tfs)


counts <- table(pwm_prom_$target)
mot_promN <- ggplot(data.frame(counts), aes(x =Freq)) + 
  geom_histogram(fill = "#8497B0", alpha=1) + ggtitle("Number of distinct TFBS per promoter") +
  labs(subtitle ="Among the 70 nitrate responsive TFs with a PWM")+ xlab("") + ylab("Count") + 
  geom_vline(xintercept = mean(counts), col = rgb(174/255, 90/255,30/255), size = 2)+
  theme_pubr()

counts <- table(pwm_prom_$TF)
cib_promN <- ggplot(data.frame(counts), aes(x =Freq)) + 
  geom_histogram(fill = "gray", alpha=1) + ggtitle("Number of distinct target TFBS per TF")+
  labs(subtitle = "In the promoters of the 1426 nitrate responsive genes") + xlab("") + ylab("Count") + 
  geom_vline(xintercept = mean(counts), col = rgb(174/255, 90/255,30/255), size = 2) +
  theme_pubr()

figure <- mot_promN + cib_promN + 
  patchwork::plot_annotation(tag_levels = "a"); figure

ggexport(figure, filename = "results/prior_network.pdf", width = 10, height = 4)






