## PWM search scores in the promoters of Arabidopsis

The aim is to run FIMO, a search for TFBS, in all Arabidopsis promoters.

Those commands should be run from a folder which contains the Arabidopsis GFF3 from TAIR10.

The following analyses require bedtools and a docker image of FIMO installed.

To create the file that gives the TSS of each AGI (start of the mRNA transcript in the GFF, taking into account the reading direction) : 

``` awk 'BEGIN{OFS="\t"}{if($3=="mRNA"){$3="tss";if($7 == "+"){$5=$4}else{$4=$5};print $0}}' TAIR10_GFF3_genes.gff > tss_file.gff``` 

Get promoter bed intervals :

``` awk 'BEGIN{OFS="\t"}{if($7 == "+"){$5=$5+200;$4=$4-1000}else{$5=$5+1000;$4=$4-200};print $0}' tss_file.gff > promoter_file.gff``` 

``` awk 'BEGIN{OFS="\t"}{$3=$9;print $0}' promoter_file.gff > promoter_file_renamed.gff``` 

Remove duplicate promoters from gene isoforms : 

``` sort -k1,1 -k4,4 -k5,5 promoter_file_renamed.gff > promoter_file_renamed_no_duplicates.gff``` 

Get promoter sequences :

``` bedtools getfasta -fi TAIR10_chr_all.fasta -fo TAIR10_promoters_no_duplicates.fasta -bed promoter_file_renamed_no_duplicates.gff -name``` 

Launch TFBS search from DAPSeq and JASPAR databases with FIMO, that is run from a docker container :

``` sudo docker run --rm -v /data/ocassan_data/PWM:/home/meme --user `id -u`:`id -g` memesuite/memesuite fimo --text --max-strand --max-stored-scores 10000000 AllJASPAR2022_PWM_MEME_Arabidopsis.txt TAIR10_promoters_no_duplicates.fasta > fimo_output_intron_prom_pv_max_strand.tsv``` 

``` sudo docker run --rm -v /data/ocassan_data/PWM:/home/meme --user `id -u`:`id -g` memesuite/memesuite fimo --text --max-strand --max-stored-scores 10000000 DAPSeqMEME.txt TAIR10_promoters_no_duplicates.fasta > fimo_output_intron_prom_pv_max_strand_DAPSeq.tsv``` 

Can also be done on not only promoters, but introns as well (concatenated fastas).

The csv outputs from those two commands were then imported in R and merged to form the TFBS matrix for all regulators and genes in Arabidopsis, saved as `pwm_prom_jaspar_dap.rdata`.
