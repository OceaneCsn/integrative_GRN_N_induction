# this scripts gets the number of introns and TSSs per gene in Arabidopsis
# based on the TAIR10 GFF3

library(tidyverse)

gff <- read.table("data/TAIR10_GFF3_genes.gff", h=F, sep = '\t')
# downloaded from 
# https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff

# gets the number of introns per gene (nb exons -1)
introns <- gff %>%
  filter(V3 =="exon") %>%
  mutate(V9 = str_split_fixed(str_split_fixed(V9, "=", 2)[,2], "\\.", 2)[,1]) %>%
  group_by(V9) %>%
  summarise(n_introns = n()-1)

# gets the number of transcripts per gene (number of TSSs according to charles)
tsss <- gff %>%
  filter(V3 =="mRNA") %>%
  mutate(V9 = str_split_fixed(str_split_fixed(V9, "=", 2)[,2], "\\.", 2)[,1]) %>%
  group_by(V9) %>%
  summarise(n_transcripts = n())

gene_structure <- left_join(tsss, introns, by = "V9")
save(gene_structure, file = "rdata/gene_structure.rdata")
