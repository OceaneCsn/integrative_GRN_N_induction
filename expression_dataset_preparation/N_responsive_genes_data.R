library(DIANE)
library(splines)
library(limma)
library(tidyverse)
library(ggpubr)

data <- read.csv("data/Arabidopsis_Varala_2018_nitrate_response_roots_with_control.txt",
                 sep = '\t', row.names = "Gene",
                 h = T, check.names = F)

colnames(data) <- str_replace(colnames(data), '-', '_')

# normalisation and filtering
tcc_object <-
  normalize(data, str_split_fixed(colnames(data), '_', 2)[, 1], iteration = FALSE)
threshold = 10 * ncol(data)
tcc_object <- DIANE::filter_low_counts(tcc_object, threshold)
normalized_counts <- TCC::getNormalizedData(tcc_object)
pca <- draw_PCA(normalized_counts)
ggexport(pca, filename = "results/supp_figures/PCA.pdf", width = 15, height = 12)


## samples 10min and 15min seem really odd, like generated in another experiment or something
normalized_counts <- normalized_counts[,!str_detect(colnames(normalized_counts), '10|15')]
draw_PCA(normalized_counts)

##### DEGs by N*control


N_treatment <- ifelse(str_detect(colnames(normalized_counts), "N"), "N", "C")
time <- as.numeric(substr(str_split_fixed(colnames(normalized_counts), '_', 2)[,1], start = 2, 
                          stop = length(str_split_fixed(colnames(normalized_counts), '_', 2)[,1])))

### Differential expression via a spline design
X<-ns(time, df=5)
design<-model.matrix(~X*N_treatment)
fit<-lmFit(normalized_counts,design)
fit<-eBayes(fit)
Qle0.01<-topTable(fit,coef=8:12,adjust="BH",p.value=0.001,number=20000)
degs<-rownames(Qle0.01)

# genes annotated as transcriptional regulators
load("rdata/regulators.rdata")
regulators_in_degs <- intersect(regulators, degs)

# saving all input data for network inference later
input_data <- list(counts = normalized_counts[degs,], grouped_genes = degs, 
                   grouped_regressors = regulators_in_degs)

write.table(normalized_counts[degs,], row.names = T,
            file = "results/supp_tables/normalized_expression_nitrate_genes.tsv", 
            quote = F, sep = '\t')

write.table(regulators_in_degs, row.names = F, col.names = F,
            file = "results/supp_tables/nitrate_responsive_regulators.tsv", 
            quote = F, sep = '\t')

save(input_data, file = "rdata/inference_input_N_response_varala.rdata")

