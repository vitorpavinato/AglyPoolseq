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
dim(snps_ann)

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

## FOLDER: signf_multilocusfst_poolsnp
## Save the preliminar Annotation
#write.table(x=significant_snps_window, file = 'results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign_multilocus_fst_ann-sumstats.txt',
#            quote = F, row.names = F, col.names = T)

## SAVE THE SNPeff-annR windows predicted genes
#write(unique(significant_snps_window$GENE), file = 'results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/raw_outliers_genes_snpeff-annR-windows_predicted.pep.txt')

###ADDED 26/JAN/2022
###
###
### --- ISOLATE WC_BETA OF SIGNIFICANT SCAFFOLDS
###
###

## GLOBAL VARIABLES FOR PLOTS AND FILTERING
## FROM poolfstat_poolsnp_analyses
dt.1.fst.scan_gFST = 0.004362626
dt.1.fst.scan_thresh = 0.1930718 
dt.1.fst.scan.025_thres = 0.25

sign_scaffolds <- snps_ann[snps_ann$CHROM %in% unique(fstsnps$CHROM), ]
table(sign_scaffolds$CHROM)

## PLOT 
## Multi-manhattam plots of W&C Beta's for all significant scaffolds
## FOLDER: signf_multilocusfst_poolsnp
## Blueish line: smoothed WC Beta
## Yellowish line: smoothed pi within
## Red dashed line: multilocus WC FST 0.25 (genome scan)
## Blue dashed line: multilocus WC FST 99% quantile (genome scan)
pdf("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/significant_intralocus_WCBetas_scaffolds.pdf",  # File name
    width = 8.5, height = 11,                                                                                           # Width and height in inches
    bg = "white",                                                                                                       # Background color
    colormodel = "cmyk",                                                                                                # Color model (cmyk is required for most publications)
)
par(mar=c(5,3.5,4,3.2)+.1, mfrow=c(4, 4)) # 16 scaffolds
for (s in seq_along(unique(sign_scaffolds$CHROM)))
{
  chr_name <- unique(sign_scaffolds$CHROM)[s]
  n <- sign_scaffolds[sign_scaffolds$CHROM %in% chr_name,]
  
  ##Plot W&C Beta points and lowess line
  plot(n$WC_Beta ~ as.numeric(n$POS/1e6), 
       xlab="Position (in Mbp)", ylab="", main=chr_name,
       col=adjustcolor("black", alpha.f = 0.5),
       pch=19, cex.lab=1.4, cex.axis=1.2)
  abline(h=c(dt.1.fst.scan.025_thres, dt.1.fst.scan_thresh), 
         lty=c(2,2), col=c("red","blue"), lwd=c(2,2))
  lines(lowess(n$POS/1e6, n$WC_Beta, f = 1/3), col="#71DFE7", lwd=3)
  mtext(expression(paste(italic("W&C")[beta])), side=2, line=1.5)
  
  
  ##Add pi_within lowess line
  par(new = TRUE)
  plot(n$pi_within ~ as.numeric(n$POS/1e6), 
       pch="", axes = FALSE, bty = "n", xlab = "", ylab = "", 
       cex.lab=1.4, cex.axis=1.2)
  axis(side=4, at = pretty(range(n$pi_within)))
  lines(lowess(n$POS/1e6, n$pi_within, f = 1/3), col="#FFC900", lwd=3)
  mtext(expression(paste(italic(pi)[within])), side=4, line=2)
  
}
dev.off()

## PLOT 
## Multi-manhattam plots of POOLFSTAT FSTs for all significant scaffolds
## FOLDER: signf_multilocusfst_poolsnp
## Blueish line: smoothed POOLFSTAT FSTs
## Yellowish line: smoothed pi within
## Red dashed line: multilocus WC FST 0.25 (genome scan)
## Blue dashed line: multilocus WC FST 99% quantile (genome scan)
pdf("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/significant_intralocus_NewFst_scaffolds.pdf",  # File name
    width = 8.5, height = 11,                                                                                           # Width and height in inches
    bg = "white",                                                                                                       # Background color
    colormodel = "cmyk",                                                                                                # Color model (cmyk is required for most publications)
)
par(mar=c(5,3.5,4,3.2)+.1, mfrow=c(4, 4)) # 16 scaffolds
for (s in seq_along(unique(sign_scaffolds$CHROM)))
{
  chr_name <- unique(sign_scaffolds$CHROM)[s]
  n <- sign_scaffolds[sign_scaffolds$CHROM %in% chr_name, c(1,2,26,30)]
  m <- n[complete.cases(n), ]
  
  ##Plot W&C Beta points and lowess line
  plot(m$POOLFSTAT_FST ~ as.numeric(m$POS/1e6), 
       xlab="Position (in Mbp)", ylab="", main=chr_name,
       col=adjustcolor("black", alpha.f = 0.5),
       pch=19, cex.lab=1.4, cex.axis=1.2)
  abline(h=c(dt.1.fst.scan.025_thres, dt.1.fst.scan_thresh, dt.1.fst.scan_gFST), 
         lty=c(2,2,2), col=c("red","blue","lightgreen"), lwd=c(2,2,2))
  lines(lowess(m$POS/1e6, m$POOLFSTAT_FST, f = 1/3), col="#71DFE7", lwd=3)
  mtext(expression(paste(italic(F)[ST])), side=2, line=1.5)
  
  
  ##Add pi_within lowess line
  par(new = TRUE)
  plot(n$pi_within ~ as.numeric(n$POS/1e6), 
       pch="", axes = FALSE, bty = "n", xlab = "", ylab = "", 
       cex.lab=1.4, cex.axis=1.2)
  axis(side=4, at = pretty(range(n$pi_within)))
  lines(lowess(n$POS/1e6, n$pi_within, f = 1/3), col="#FFC900", lwd=3)
  mtext(expression(paste(italic(pi)[within])), side=4, line=2)
  
}
dev.off()

## FILTER 1
## SNP-POOLFSTAT FST > 0.25 (multilocus fst)
sign_snp_fst <- sign_scaffolds[which(sign_scaffolds$POOLFSTAT_FST > dt.1.fst.scan.025_thres), c(1:2)]

# Save to a file
write.table(sign_snp_fst,
            file="results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/signf_snp_fst.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)

## FILTER 2
## SNP-W&C Beta > 0.25 (multilocus fst)
sign_WC_beta <- sign_scaffolds[which(sign_scaffolds$WC_Beta > dt.1.fst.scan.025_thres), c(1:2)]

write.table(sign_WC_beta,
            file="results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/signf_WC_beta.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)



##
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