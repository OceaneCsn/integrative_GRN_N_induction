load(file = "rdata/inference_input_N_response_varala.rdata")
load('rdata/pwm_prom_jaspar_dap.rdata')

source('utils/utils.R')

genes <- input_data$grouped_genes
tfs <- input_data$grouped_regressors
counts <- input_data$counts

tfs_spread <- ungroup_tfs(tfs)
genes_spread <- ungroup_tfs(genes)
pwm_prom_ <- pwm_prom[pwm_prom$target %in% genes_spread & 
                        pwm_prom$TF %in% tfs_spread, ]

known_tfs <- intersect(pwm_prom$TF, tfs_spread)


counts <- table(pwm_prom$target)
mot_prom <- ggplot(data.frame(counts), aes(x =Freq)) + 
  geom_histogram(fill = "#8497B0", alpha=0.65) + ggtitle("Number of distinct TFBS per promoter")+
  labs(subtitle ="Genome-wide")+ xlab("") + ylab("Count") + 
  geom_vline(xintercept = mean(counts), col = rgb(174/255, 90/255,30/255), size = 2)

counts <- table(pwm_prom$TF)
cib_prom <- ggplot(data.frame(counts), aes(x =Freq)) + 
  geom_histogram(fill = "gray", alpha=0.65) + ggtitle("Number of distinct TFBS per TF")+
  labs(subtitle ="Genome-wide")+ xlab("") + ylab("Count") + 
  geom_vline(xintercept = mean(counts), col = rgb(174/255, 90/255,30/255), size = 2)

counts <- table(pwm_prom_$target)
mot_promN <- ggplot(data.frame(counts), aes(x =Freq)) + 
  geom_histogram(fill = "#8497B0", alpha=1) + ggtitle("Number of distinct TFBS per promoter") +
  labs(subtitle ="Among the 1426 N responsive genes")+ xlab("") + ylab("Count") + 
  geom_vline(xintercept = mean(counts), col = rgb(174/255, 90/255,30/255), size = 2)

counts <- table(pwm_prom_$TF)
cib_promN <- ggplot(data.frame(counts), aes(x =Freq)) + 
  geom_histogram(fill = "gray", alpha=1) + ggtitle("Number of distinct TFBS per TF")+
  labs(subtitle = "Among the 1426 N responsive genes") + xlab("") + ylab("Count") + 
  geom_vline(xintercept = mean(counts), col = rgb(174/255, 90/255,30/255), size = 2)

figure <- mot_prom + cib_prom + mot_promN + cib_promN + 
  patchwork::plot_annotation(tag_levels = "a"); figure

ggexport(figure, filename = "results/prior_network.pdf", width = 9, height = 7)






