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
setwd("/fs/project//PAS1554/aphidpool")
renv::restore()
#renv::snapshot()

# Remove last features
rm(list=ls())
ls()

# Load R libraries
library(tidyverse)
library(gplots)

##### PART 1: MULTILOCUS-FST SIGNIFICANT SNPS
#####----------------------------------------

###
###
### --- UPLOAD FILES ---
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
#write.table(x=significant_snps_window, file = 'results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign_multilocus_fst_annot_table.txt',
#            quote = F, row.names = F, col.names = T)

## SAVE THE SNPeff-annR windows predicted genes
#write(unique(significant_snps_window$GENE), file = 'results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/raw_outliers_genes_snpeff-annR-windows_predicted.pep.txt')

###ADDED 26/JAN/2022
###
###
### --- ISOLATE WC_BETA OF SIGNIFICANT SCAFFOLDS --- 
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
## ALL SNP-POOLFSTAT FST > 0.25 WITHIN SIGNIFICANT SCAFFOLDS IDENTIFIED WITH MULTILOCUS FST
sign_snp_fst <- sign_scaffolds[which(sign_scaffolds$POOLFSTAT_FST > dt.1.fst.scan.025_thres), c(1:2)]

# Save to a file
write.table(sign_snp_fst,
            file="results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/loci/signf_snp_fst.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)

## FILTER 2
## ALL  SNP-W&C Beta > 0.25  WITHIN SIGNIFICANT SCAFFOLDS IDENTIFIED WITH MULTILOCUS FST
sign_WC_beta <- sign_scaffolds[which(sign_scaffolds$WC_Beta > dt.1.fst.scan.025_thres), c(1:2)]

write.table(sign_WC_beta,
            file="results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/loci/signf_WC_beta.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)


## FILTER 3
## ONLY SNP-W&C Beta > 0.25 WITHIN 10Kbp WINDOWS WITHIN SIGNIFICANT SCAFFOLDS IDENTIFIED WITH MULTILOCUS FST
sign_WC_beta_window <- significant_snps_window[which(significant_snps_window$WC_Beta > dt.1.fst.scan.025_thres), c(1:2)]

write.table(sign_WC_beta_window,
            file="results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/loci/signf_WC_beta_window.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)

## Savepoint_1
#save.image("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")
#load("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")


##### PART 2: ADDITIONAL EVIDENCE TO SUPPORT SIGNIFICANT SNPS
#####-------------------------------------------------------
###
###
### --- SNPS IN COMMON AMONG GENOME SCAN METHODS ---
###
###
### METHODS::MULTILOCUS-FST; BAYPASS-XtX; BAYPASS-BFmc; SELESTIM-KLD

## Upload files
lociids_signf_WC_beta_window <- scan(file = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/loci/lociIDs_signf_WC_beta_window.txt",what = "character")
lociids_signf_snps_core0_xtx <- scan(file = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/lociIDs_signf_snps_core0_podsQuant99.txt", what = "character")
lociids_signf_snps_aux0_bfmc <- scan(file = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/lociIDs_signf_snps_aux0_podsQuant99.txt", what = "character")
lociids_signf_snps_common_kld <- scan(file = "results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/lociIDs_signf_snps_kld_common.txt", what = "character")

## EXPORT THE LIST GENES ASSOCIATED WITH OUTLIER LIST
lociids_lists <- list(fst = unique(lociids_signf_WC_beta_window),
                      xtx = lociids_signf_snps_core0_xtx,
                      bfmc = lociids_signf_snps_aux0_bfmc,
                      kld = lociids_signf_snps_common_kld)

venn.obj <- venn(lociids_lists, intersections=TRUE)

venn.obj.intersections <- attr(venn.obj, "intersections")

## TWO METHODS
fst.xtx <- venn.obj.intersections$`fst:xtx`
fst.kld <- venn.obj.intersections$`fst:kld`
xtx.bfmc <- venn.obj.intersections$`xtx:bfmc`
xtx.kld <- venn.obj.intersections$`xtx:kld`
bfmc.kld <- venn.obj.intersections$`bfmc:kld`

## THREE METHODS
fst.xtx.kld <- venn.obj.intersections$`fst:xtx:kld`
xtx.bfmc.kld <- venn.obj.intersections$`xtx:bfmc:kld`

list_sign_methods <- list(fst_xtx = fst.xtx,
                          fst_kld = fst.kld,
                          xtx_bfmc = xtx.bfmc,
                          xtx_kld = xtx.kld,
                          bfcm_kld = bfmc.kld,
                          fst_xtx_kld = fst.xtx.kld,
                          xtx_bfmc_kld = xtx.bfmc.kld)

## EXPORT THE LIST OF COMMON SNPS IN AT LEAST TWO TESTS
workingdir = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/loci/"
for (i in 1:length(list_sign_methods))
{
  tmp <- list_sign_methods[[i]]
  
  write.table(tmp,
              file = paste0(workingdir, "signf_snps_methods", "_", names(list_sign_methods[i]), ".txt"),
              col.names = F, row.names = F, quote = F)
  rm(tmp)
}

## EXPORT THE ANNOTATION OF SNPS PRESENT IN THOSE LISTS
workingdir = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/loci/"
for (i in 1:length(list_sign_methods))
{
  tmp <- snps_ann[paste0(snps_ann$CHROM, "_", snps_ann$POS) %in% list_sign_methods[[i]], ]
  
  write.table(tmp,
              file = paste0(workingdir, "signf_snps_methods_anno", "_", names(list_sign_methods[i]), ".txt"),
              col.names = F, row.names = F, quote = F)
  rm(tmp)
}

## Savepoint_2
#save.image("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")
#load("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")

##### PART 3: FUNCTIONAL ANNOTATION OF SIGNIFICANT GENES
#####---------------------------------------------------
###
###
### --- GENES IN COMMON AMONG GENOME SCAN METHODS ---
###
###

## GET GENE IDS ASSOCIATED WITH SIGNIFICANT SNPs
## THIS MIGHT CONTAIN GENES THAT HAS WITHIN OR OUTSIDE SNPS
## THE "WINDOW" SIZE MAY VARY
workingdir = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/"
for (i in 1:length(lociids_lists))
{
  tmp <- snps_ann[paste0(snps_ann$CHROM, "_", snps_ann$POS) %in% lociids_lists[[i]], "GENE"]
  res <- unique(sort(do.call("c", strsplit(tmp, split="-"))))
  
  write.table(res,
              file = paste0(workingdir, "geneIDs_signf", "_", names(lociids_lists[i]), ".txt"),
              col.names = F, row.names = F, quote = F)
  rm(tmp)
  rm(res)
}

## GET GENES IN COMMON AMONG METHODS
## Upload files
genesids_fst <- scan(file = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/geneIDs_signf_fst.txt",what = "character")
genesids_xtx <- scan(file = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/geneIDs_signf_xtx.txt", what = "character")
genesids_bfmc <- scan(file = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/geneIDs_signf_bfmc.txt", what = "character")
genesids_kld <- scan(file = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/geneIDs_signf_kld.txt", what = "character")

## EXPORT THE LIST GENES ASSOCIATED WITH OUTLIER LIST
genesids_lists <- list(fst = genesids_fst,
                       xtx = genesids_xtx,
                       bfmc = genesids_bfmc,
                       kld = genesids_kld)

venn.obj.genes <- venn(genesids_lists, intersections=TRUE)

venn.obj.genes.intersections <- attr(venn.obj.genes, "intersections")

## TWO METHODS
g_fst.xtx <-  venn.obj.genes.intersections$`fst:xtx`
g_fst.kld <-  venn.obj.genes.intersections$`fst:kld`
g_xtx.bfmc <- venn.obj.genes.intersections$`xtx:bfmc`
g_xtx.kld <-  venn.obj.genes.intersections$`xtx:kld`
g_bfmc.kld <- venn.obj.genes.intersections$`bfmc:kld`

## THREE METHODS
g_fst.xtx.kld <- venn.obj.genes.intersections$`fst:xtx:kld`
g_xtx.bfmc.kld <- venn.obj.genes.intersections$`xtx:bfmc:kld`

## FOUR METHODS
g_fst.xtx.bfmc.kld <- venn.obj.genes.intersections$`fst:xtx:bfmc:kld`

## OUTPUT A TABLE CONTAINING THE GENES IDENTIFIED BY AT LEAST TWO METHODS
gene_table <- data.frame(geneID = c(g_fst.xtx, g_fst.kld, g_xtx.bfmc, g_xtx.kld, 
                                    g_bfmc.kld, g_fst.xtx.kld, g_xtx.bfmc.kld, g_fst.xtx.bfmc.kld),
                         methods = c(rep("FST_XtX",14), rep("FST_KLD",2), rep("XtX_BFmc",135), rep("XtX_KLD",762), 
                                     rep("BFmc_KLD",136), rep("FST_XtX_KLD",8), rep("XtX_BFmc_KLD",220), rep("FST_XtX_BFmc_KLD",7)),
                         signf_level = c(rep(2, 1049), rep(3, 228), rep(4, 7))
                         )

## Add the scaffold information
gene_table_scaffolds <- snps_ann[snps_ann$GENE %in% gene_table$geneID, c("CHROM", "GENE")] 
colnames(gene_table_scaffolds) <- c("Scaffold","geneID")

gene_table_scaffolds <- gene_table_scaffolds[!duplicated(gene_table_scaffolds$geneID), ]

gene_table <- merge(gene_table, gene_table_scaffolds, by.x = 1, by.y=2, all.x=T)
gene_table <- gene_table[, c(1, 4, 2, 3)]

barplot(table(gene_table$signf_level))
dim(gene_table)
length(unique(gene_table$geneID)) #1284

## EXPORT TABLE
workingdir = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/"
write.table(gene_table,
            file = paste0(workingdir, "signf_common_genes_table", ".txt"),
            col.names = T, row.names = F, quote = F)

## Savepoint_3
#save.image("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")
#load("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")

###
###
### --- INCLUDE FUNCTIONAL ANNOTATION OF CANDIDATE GENES ---
###
###

## BIPAA ANNOTATION
bipaa_annot <- read.table(file="ref/raw/func_annot/bipaa/blast2go.tsv", header = T, sep = "\t", na.strings = "")
bipaa_annot["geneID"] <- do.call(rbind, strsplit(bipaa_annot$Sequence.Name, ".p"))[,1]
bipaa_annot <- bipaa_annot[,c(14, 1:13)]

bipaa_annot <- bipaa_annot[!duplicated(bipaa_annot$geneID), ]

# GO TERMS
bipaa_gaf <- read.table(file="ref/raw/func_annot/bipaa/blast2go.gaf", header = T, na.strings = "")
colnames(bipaa_gaf) <- c("Sequence.Name", "go_term")

bipaa_gaf["geneID"] <- do.call(rbind, strsplit(bipaa_gaf$Sequence.Name, ".p"))[,1]
bipaa_gaf <- bipaa_gaf[,c(3, 1:2)]

## ROGER KEGG TERMS 
roger_kegg <- read.table(file = "ref/raw/func_annot/rogerio/Kegg.txt", header = T, sep = "\t", fill = T, na.strings = "")
roger_kegg["geneID"] <- do.call(rbind, strsplit(roger_kegg$gene, ".t"))[,1]
roger_kegg <- roger_kegg[,c(3, 1:2)]

roger_kegg <- roger_kegg[!duplicated(roger_kegg$geneID), ]

## MERGE "gene_table" with annotations
gene_table_annot <- merge(gene_table, bipaa_annot, by = 1)
gene_table_annot <- merge(gene_table_annot, roger_kegg, by = 1)

## Savepoint_4
#save.image("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")
#load("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")

###
###
### --- INCLUDE THE GENE EXPRESSION AND THE ADAPTIVE RESPONSE ---
###
###

## File path of the expression analysis
plastic_biotype1_file <- "/fs/project/PAS1554/aphidrnaseq/results/response_plastic_b1_RS.txt"
plastic_biotype4_file <- "/fs/project/PAS1554/aphidrnaseq/results/response_plastic_b4_RS.txt"
evolved_susceptible_file <- "/fs/project/PAS1554/aphidrnaseq/results/response_evolved_b1b4_S.txt"
evolved_resistance_file <- "/fs/project/PAS1554/aphidrnaseq/results/response_evolved_b1b4_R.txt"

file_names <- c(plastic_biotype1_file, plastic_biotype4_file, 
                evolved_susceptible_file, evolved_resistance_file)

## Import
data_std.response <- lapply(file_names, function(x){read.table(file= x, header=T)})

## Change names
names(data_std.response) <- c("plastic_b1_RS",
                              "plastic_b4_RS",
                              "evolved_b1b4_S",
                              "evolved_b1b4_R")

## Simplify the data
plastic_b1_RS <- data_std.response$plastic_b1_RS[,c("name", "foldChange", "FDR", "adaptation")]
colnames(plastic_b1_RS) <- c("geneID", "plastic_b1_RS_foldChange", "plastic_b1_RS_b1_FDR", "plastic_b1_RS_adaptation")

plastic_b4_RS <- data_std.response$plastic_b4_RS[,c("name", "foldChange", "FDR", "adaptation")]
colnames(plastic_b4_RS) <- c("geneID", "plastic_b4_RS_foldChange", "plastic_b4_RS_FDR", "plastic_b4_RS_adaptation")

evolved_b1b4_S <- data_std.response$evolved_b1b4_S[,c("name", "foldChange", "FDR", "adaptation")]
colnames(evolved_b1b4_S) <- c("geneID", "evolved_b1b4_S_foldChange", "evolved_b1b4_S_FDR", "evolved_b1b4_S_adaptation")

evolved_b1b4_R <- data_std.response$evolved_b1b4_R[,c("name", "foldChange", "FDR", "adaptation")]
colnames(evolved_b1b4_R) <- c("geneID", "evolved_b1b4_R_foldChange", "evolved_b1b4_R_FDR", "evolved_b1b4_R_adaptation")

## Make a list 
data_std.response.simple <- list(plastic_b1_RS, plastic_b4_RS, 
                                 evolved_b1b4_S, evolved_b1b4_R)

## Combine the datasets
data_std.response.merged <- Reduce(function(x,y) {merge(x,y, by=1, all = TRUE)}, data_std.response.simple)


## Now merge with the table with significant genes
gene_table_annot.response <- merge(gene_table_annot, data_std.response.merged, by = 1, all.x = TRUE)

## Sort by the number of significant tests
gene_table_annot.response.sorted <- gene_table_annot.response[order(gene_table_annot.response$signf_level, decreasing = T), ]

## EXPORT TABLE
workingdir = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/gene/"
write.table(gene_table_annot.response.sorted,
            file = paste0(workingdir, "signf_common_genes_table_annotated_reponse", ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")


## MAKE SOME PLOTS
n_genes=length(unique(gene_table_annot.response.sorted$geneID))

## % OF PLASTIC REPONSE IN EACH BIOTYPE
df_1 <- data.frame(biotype = c(rep("B1", 2), rep("B4",2)),
                  response = c(rep(c("Constitutive", "Plastic"),2)),
                  value = c(99.343832, 0.656168, 99.0909091,0.9090909))

p_1 <- ggplot(df_1 , aes(fill=response, y=value, x=biotype)) + 
       geom_bar(position="dodge",stat="identity") +
       labs(fill="Response") +
       xlab("") +
       scale_x_discrete(labels = c("Avirulent", "Virulent")) +
       ylab("%") +
       ylim(c(0,100)) +
       theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                          axis.line = element_line(colour = "black"),
                          axis.text = element_text(size = 14),
                          axis.title.y=element_text(size=16),
                          axis.title.x=element_text(size=16),
                          axis.text.x = element_text(size=12, angle = 0, vjust = 0.5, hjust=0.5))
p_1

## % OF EVOLVED REPONSE IN EACH PLANT
df_2 <- data.frame(plant = c(rep("Susceptible", 2), rep("Resistant",2)),
                   response = c(rep(c("Constitutive", "Evolved"),2)),
                   value = c(68.14965, 31.85035, 72.77433, 27.22567))

p_2 <- ggplot(df_2 , aes(fill=response, y=value, x=plant)) + 
  geom_bar(position="dodge",stat="identity") +
  labs(fill="Response") +
  xlab("") +
  scale_x_discrete(labels = c("Resistant", "Susceptible")) +
  ylab("%") +
  ylim(c(0,100)) +
  theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(size = 14),
                     axis.title.y=element_text(size=16),
                     axis.title.x=element_text(size=16),
                     axis.text.x = element_text(size=12, angle = 0, vjust = 0.5, hjust=0.5))
p_2

## Savepoint_5
#save.image("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")
#load("results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/sign.multilocusfst.workspace23Feb22.RData")








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