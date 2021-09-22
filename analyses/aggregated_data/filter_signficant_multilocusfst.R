##########################################################
#             ADD MULTILOCUS-FST VALUES TO THE           # 
#                  SNP ANNOTATION TABLE                  #
##########################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

###
###
### ---- SETUP ANALYZES ENVIRONMENT ----
###
###

# Recover R-renv environment
setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()

# Remove last features
rm(list=ls())
ls()

# Load R libraries
library(tidyverse)

###
###
### ---- UPLOAD FILES ----
###
###

## Load table containing the SNPs positions, annotation and summary stats (from 'genetic_diversity_poolsnp.R' analyses)
snps_ann_path = 'results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/combined_summary_statistics_annotation_table.txt'

snps_ann <- read.table(file = snps_ann_path, header = T, sep = "\t")

## Load the bed file containing the significant 'multilocus-fst snp window limits' 
snp_fstBED_path      = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst.bed'
snp_fstBED_header <- c('CHROM', 'start', 'end')

fstsnps <- read.table(file = snp_fstBED_path, sep = '\t', col.names = snp_fstBED_header) 
dim(fstsnps)

## Load the file with multilocus-fst values
snp_fstVALUES_path = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst_values.txt'

fstVALUES <- read.table(file = snp_fstVALUES_path, header = T) 
dim(fstVALUES)

## Combine the Multilocus-FST coordinates with the values

# Sanity check if the two file are in the same order
sum(paste0(fstsnps$CHROM,"_",fstsnps$end) == paste0(fstVALUES$CHRM,"_",fstVALUES$POS))

fstsnps['multilocus_fst'] <- fstVALUES$FST

## Keep only SNPs in the SNP annotation file that felt within the 10 Kbp limits each side around the multilocus-FST (20Kbp)
dist2snp=10000
significant_snps_window = NULL
for(i in seq_along(fstsnps$CHROM)){
  
  fstsnps_ <- fstsnps[i, ]
  snps_ann_ <- snps_ann[snps_ann$CHROM %in% fstsnps_$CHROM, ]
  
  snps_in_window <- snps_ann_[snps_ann_$POS >= (fstsnps_$end - dist2snp) & 
                              snps_ann_$POS < (fstsnps_$end + dist2snp), ]
  
  d <- data.frame(snps_in_window, 
                  multfstSNP=paste0(fstsnps_$start,':', fstsnps_$end), 
                  window_size = ((fstsnps_$end + dist2snp)-(fstsnps_$end - dist2snp)), 
                  multilocus_fst=fstsnps_$multilocus_fst)
  
  significant_snps_window <- rbind(significant_snps_window, d)
  
}

dim(significant_snps_window) # 801  33
head(significant_snps_window)

## Save the preliminar Annotation
write.table(x=significant_snps_window, file = 'results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign_multilocus_fst_ann-sumstats.txt',
            quote = F, row.names = F, col.names = T)

## SAVE THE SNPeff-annR windows predicted genes
write(unique(significant_snps_window$GENE), file = 'results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/raw_outliers_genes_snpeff-annR-windows_predicted.pep.txt')

## Take only the NS AND SS mutations
#res_annotation_genes <- res_annotation[which(res_annotation$FITNES_EFFECT == "NS" | res_annotation$FITNES_EFFECT == "SS"), ]
#
### Calculat the pNpS for the genes
#genes_fitness = NULL
#for (i in seq_along(unique(res_annotation_genes$GENE)))
#{
#  
#  gene_set <- res_annotation_genes[which(res_annotation_genes$GENE == unique(res_annotation_genes$GENE)[i]), ]
#  ns = sum(gene_set$FITNES_EFFECT == "NS")
#  ss = sum(gene_set$FITNES_EFFECT == "SS")
#  pnps = ns/ss
#  
#  g <- data.frame(GENE=gene_set[1, 4], pNpS=pnps)
#  genes_fitness <- rbind(genes_fitness, g)
#  
#}
#
#genes_fitness[is.infinite(genes_fitness$pNpS), 'pNpS'] <- 1
#
### Save pNpS table
#write.table(x=genes_fitness, file = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_snp_gene_pnps.txt',
#            quote = F, row.names = F, col.names = T)