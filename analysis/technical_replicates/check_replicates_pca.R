#############################################################
# Access sequencing consistency across technical replicates #
#     Allele frequency calculation from vcf and PCA         #
#############################################################

# Each library have 5 pooled individuals and they were sequenced 4-5 times (4-5 days.)
# They were the same library, with same barcode, but the sequencing was splitted in 4-5 days. 

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
## Remove last features
rm(list=ls())
ls()

## Load libraries
library(poolfstat)
library(ggpubr)
library(missMDA)
library(ade4)
library(factoextra)
library(zoo)
library(FactoMineR)
packageVersion("poolfstat") 
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

###
###
### ---- MAXIMUM LIKELIHOOD ESTIMATION FROM ALLELE COUNTS ----
###
###

### DATASET: VCF WITH REPLICATED RUNS

## ADD CITATION HERE

# First load vcf using poolfsta function vcf2pooldata from poolfstat package

ALPHA=0.75

poolsizes <- rep(10,87)
poolnames <- c("MN_BIO1_S1_140711", "MN_BIO1_S1_140722", "MN_BIO1_S1_140725", "MN_BIO1_S1_140730", "MN_BIO1_S1_AddRun", 
               "MN_BIO4_S1_140711", "MN_BIO4_S1_140722", "MN_BIO4_S1_140725", "MN_BIO4_S1_140730", 
               "ND_BIO1_S1_140711", "ND_BIO1_S1_140722", "ND_BIO1_S1_140725", "ND_BIO1_S1_140730", 
               "ND_BIO4_S1_140711", "ND_BIO4_S1_140722", "ND_BIO4_S1_140725", "ND_BIO4_S1_140730", 
               "NW_BIO1_S1_140711", "NW_BIO1_S1_140722", "NW_BIO1_S1_140725", "NW_BIO1_S1_140730", 
               "NW_BIO1_S2_140711", "NW_BIO1_S2_140722", "NW_BIO1_S2_140725", "NW_BIO1_S2_140730", 
               "NW_BIO4_S1_140711", "NW_BIO4_S1_140722", "NW_BIO4_S1_140725", "NW_BIO4_S1_140730", "NW_BIO4_S1_AddRun", 
               "NW_BIO4_S2_140711", "NW_BIO4_S2_140722", "NW_BIO4_S2_140725", "NW_BIO4_S2_140730", 
               "NW_BIO4_S3_140711", "NW_BIO4_S3_140722", "NW_BIO4_S3_140725", "NW_BIO4_S3_140730", 
               "PA_BIO1_S1_140711", "PA_BIO1_S1_140722", "PA_BIO1_S1_140725", "PA_BIO1_S1_140730", 
               "PA_BIO4_S1_140711", "PA_BIO4_S1_140722", "PA_BIO4_S1_140725", "PA_BIO4_S1_140730", 
               "PA_BIO4_S2_140711", "PA_BIO4_S2_140722", "PA_BIO4_S2_140725", "PA_BIO4_S2_140730",
               "WI_BIO1_S1_140711", "WI_BIO1_S1_140722", "WI_BIO1_S1_140725", "WI_BIO1_S1_140730", 
               "WI_BIO1_S2_140711", "WI_BIO1_S2_140722", "WI_BIO1_S2_140725", "WI_BIO1_S2_140730", "WI_BIO1_S2_AddRun", 
               "WI_BIO4_S1_140711", "WI_BIO4_S1_140722", "WI_BIO4_S1_140725", "WI_BIO4_S1_140730", 
               "WI_BIO4_S2_140711", "WI_BIO4_S2_140722", "WI_BIO4_S2_140725", "WI_BIO4_S2_140730", 
               "WI_BIO4_S3_140711", "WI_BIO4_S3_140722", "WI_BIO4_S3_140725", "WI_BIO4_S3_140730", 
               "WO_BIO1_S1_140711", "WO_BIO1_S1_140722", "WO_BIO1_S1_140725", "WO_BIO1_S1_140730", 
               "WO_BIO4_S1_140711", "WO_BIO4_S1_140722", "WO_BIO4_S1_140725", "WO_BIO4_S1_140730", 
               "WO_BIO4_S2_140711", "WO_BIO4_S2_140722", "WO_BIO4_S2_140725", "WO_BIO4_S2_140730", 
               "WO_BIO4_S3_140711", "WO_BIO4_S3_140722", "WO_BIO4_S3_140725", "WO_BIO4_S3_140730")

dt.repl <- vcf2pooldata(
          vcf.file = "vcf/technical_replicates/minmaxcov_3_99/aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz",
          poolsizes = poolsizes,
          poolnames = poolnames,
          min.cov.per.pool = -1,
          min.rc = 1,
          max.cov.per.pool = 1e+06,
          min.maf = 0.05,
          remove.indels = FALSE,
          nlines.per.readblock = 1e+06
          )
dt.repl
# * * * PoolData Object * * * 
# * Number of SNPs   =  40105 
# * Number of Pools  =  87 

nbr.gene.copies = 10
nbr.loci = dt.repl@nsnp

###
###
### ---- PCA ON ML ALLELE FREQUENCIES OF TECHNICAL REPLICATES----
###
###

dt.repl.imputedRefMLCount <- imputedRefMLCount(dt.repl)
dt.repl.imputedRefMLFreq <- dt.repl.imputedRefMLCount[[1]]/dt.repl.imputedRefMLCount[[3]]

## Price et al 2010 - vcov of matrix individuals x markers

# Prepare the vectors: color, points, names
# for the plot and for the legends
run.colors <- c(rep("#0000FF",5), rep("#FF0000",4), #MN
                rep("#8282FF",4), rep("#FF8282",4), #ND
                rep("#10B717",4), rep("#64C366",4), rep("#FD63CB",5), rep("#F98AD1",4), rep("#F4A8D6",4), #NW
                rep("#B78560",4), rep("#9E5C00",4), rep("#7D3200",4), #PA
                rep("#B2B2FF",4), rep("#D5D5FF",5), rep("#FFB2B2",4), rep("#FFD5D5",4), rep("#FFF2F2",4), #WI
                rep("#90CE91",4), rep("#EFC0DB",4), rep("#E9D2DF",4), rep("#E4DFE2",4) #WO
                )

lib.pch <- c(c(8, 15, 17, 18, 19), c(8, 15, 17, 18),
             c(8, 15, 17, 18), c(8, 15, 17, 18),
             c(8, 15, 17, 18), c(8, 15, 17, 18), c(8, 15, 17, 18, 19), c(8, 15, 17, 18), c(8, 15, 17, 18),
             c(8, 15, 17, 18), c(8, 15, 17, 18), c(8, 15, 17, 18),
             c(8, 15, 17, 18), c(8, 15, 17, 18, 19), c(8, 15, 17, 18), c(8, 15, 17, 18), c(8, 15, 17, 18),
             c(8, 15, 17, 18), c(8, 15, 17, 18), c(8, 15, 17, 18), c(8, 15, 17, 18))

legend.colors <- c("#0000FF", "#FF0000", #MN
                   "#8282FF", "#FF8282", #ND
                   "#10B717", "#64C366", "#FD63CB", "#F98AD1", "#F4A8D6", #NW
                   "#B78560", "#9E5C00", "#7D3200", #PA
                   "#B2B2FF", "#D5D5FF", "#FFB2B2", "#FFD5D5", "#FFF2F2", #WI
                   "#90CE91", "#EFC0DB", "#E9D2DF", "#E4DFE2", #WO
                   rep("black", 5)
                   )

legend.names <- c("MN-B1.1", "MN-B4.1",
                  "ND-B1.1", "ND-B4.1", 
                  "NW-B1.1", "NW-B1.2", "NW-B4.1", "NW-B4.2", "NW-B4.3", 
                  "PA-B1.1", "PA-B4.1", "PA-B4.2", 
                  "WI-B1.1", "WI-B1.2", "WI-B4.1", "WI-B4.2", "WI-B4.3", 
                  "WO-B1.1", "WO-B4.1", "WO-B4.2", "WO-B4.3",
                  paste0("r_",1:5)
                  )


legend.pch <- c(rep(19, 21), c(8, 15, 17, 18, 19))

W<-scale(t(dt.repl.imputedRefMLFreq), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")

l1<- 100*lambda[1]/sum(lambda) 
l2<- 100*lambda[2]/sum(lambda) 

# PCA plot
pdf("results/technical_replicates/minmaxcov_3_99/pca_imputed_allelefreq_check_repl.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=run.colors,
     xlim = c(-0.25, 0.2), ylim = c(-0.35, 0.35),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"), 
     cex=1.5, pch=lib.pch)
legend("topleft", 
       legend = legend.names, col = legend.colors, 
       pch = legend.pch, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Remove the _AddRun_ from the allele frequency matrix (5th technical replication)
W<-scale(t(dt.repl.imputedRefMLFreq[,-c(5,30,59)]), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")

l1<- 100*lambda[1]/sum(lambda) 
l2<- 100*lambda[2]/sum(lambda) 

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=run.colors[-c(5,30,59)],
     xlim = c(-0.25, 0.2), ylim = c(-0.35, 0.35),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"),
     cex=1.5, pch=lib.pch[-c(5,30,59)])
legend("topleft", 
       legend = legend.names[-c(26)], col = legend.colors[-c(26)], 
       pch = legend.pch[-c(26)], bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)

## Savepoint_1
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
############## END THIS PART

###
###
### ---- % OF MISSING DATA AND CORRELATION WITH REPLICATE COVERAGE----
###
###

#####
### Total % of of NA Markers in each replicated pool
#####
colnames(dt.repl.imputedRefMLCount[[1]]) <- poolnames
dt.repl.missing <- 100*(colSums(is.na(dt.repl.imputedRefMLCount[[1]]))/nbr.loci)

# Mean coverage of each replicated pool
dt.repl.coverage <- dt.repl@readcoverage
mean.dt.repl.coverage <- apply(dt.repl.coverage, 2, mean)

# Coverage CI 95%
mean(mean.dt.repl.coverage) +c(-1.96, +1.96)*sd(mean.dt.repl.coverage)/sqrt(87) # 8.293861 10.807788

## MISSING DATA AS A FUNCTION OF MEAN COVERAGE FROM SNP TOTAL READS
# IDENTIFY PROBLEMATIC POOLS/REPLICATES
missing_coverage_data <- data.frame(MISSING=dt.repl.missing, COV=mean.dt.repl.coverage)
missing_coverage_data[,'STATUS_COL'] <- rep(NA, dt.repl@npools)

missing_coverage_data$STATUS_COL[missing_coverage_data$MISSING > 30 & missing_coverage_data$COV < 4] <- "orange"
missing_coverage_data$STATUS_COL[missing_coverage_data$MISSING > 30 & missing_coverage_data$COV > 4] <- "red"
missing_coverage_data$STATUS_COL[missing_coverage_data$MISSING > 20 & missing_coverage_data$COV > 15] <- "purple"
missing_coverage_data$STATUS_COL[is.na(missing_coverage_data$STATUS_COL)] <- "black"

COL <- adjustcolor(missing_coverage_data$STATUS_COL, alpha.f = 0.5)

pdf("results/technical_replicates/minmaxcov_3_99/missing_read_count_poolfstat.pdf",         # File name
    width = 11, height = 8.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
plot(MISSING ~COV,
     xlab="Mean coverage (number of reads)",ylab="% of Missing Data", cex=1.5,
     col=COL, pch=19, cex.lab=1.4, cex.axis=1.2, data=missing_coverage_data)
lines(lowess(missing_coverage_data$MISSING ~ missing_coverage_data$COV, f=0.99), col='blue', lwd=2)
abline(v=4, col="black", lty=3)
abline(h=30, col="black", lty=3)
dev.off()


#####
### MEAN COVERAGE PER SCAFFOLD FROM THE READ ALIGNMENTS
#####
pipeline.folder = "dest_mapped/pipeline_output/"
sulfixe.file.name = "_meandepth_original.bam.txt"
header.names <- c("rname", "startpos", "endpos",  "numreads", 
                  "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")

mapped.files <- paste0(pipeline.folder, "/", poolnames, "/", poolnames, sulfixe.file.name)
mapped.files.list <- lapply(mapped.files, function(x){read.table(file= x, col.names = header.names, na.strings = "NA")})

## MISSING DATA AS A FUNCTION OF MEAN SCAFFOLD COVERAGE %
mapped.files.coverage <- mapped.files.list[[1]][,c(1,2,3,6)]
names(mapped.files.coverage)[4] <- poolnames[1]

c = 5 
for (l in 2:length(mapped.files.list))
{
  mapped.files.coverage[, c] <- mapped.files.list[[l]][, 6]
  names(mapped.files.coverage)[c] <- poolnames[l]
  c = c + 1
}

mean_scaffold_coverage <- colMeans(mapped.files.coverage[-c(831:1518), -c(1:3)])
mean_scaffold_coverage_data <- data.frame(MISSING=dt.repl.missing, SCAFFOLD_COV=mean_scaffold_coverage)

mean_scaffold_coverage_data['STATUS_COL'] <- rep(NA, dt.repl@npools)

mean_scaffold_coverage_data$STATUS_COL[mean_scaffold_coverage_data$MISSING > 30 & mean_scaffold_coverage_data$SCAFFOLD_COV < 25] <- "orange"
mean_scaffold_coverage_data$STATUS_COL[mean_scaffold_coverage_data$MISSING > 30 & mean_scaffold_coverage_data$SCAFFOLD_COV > 25] <- "red"
mean_scaffold_coverage_data$STATUS_COL[mean_scaffold_coverage_data$MISSING > 25 & mean_scaffold_coverage_data$SCAFFOLD_COV > 50] <- "purple"
mean_scaffold_coverage_data$STATUS_COL[is.na(mean_scaffold_coverage_data$STATUS_COL)] <- "black"

COL_SCAFFOLD_COV <- adjustcolor(mean_scaffold_coverage_data$STATUS_COL, alpha.f = 0.5)

plot(MISSING ~SCAFFOLD_COV,
     xlab="% of Scaffold covered",ylab="% of Missing Data", cex=1.5,
     col=COL_SCAFFOLD_COV, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_scaffold_coverage_data)
lines(lowess(mean_scaffold_coverage_data$MISSING ~ mean_scaffold_coverage_data$SCAFFOLD_COV, f=0.99), col='blue', lwd=2)
abline(v=25, col="black", lty=3)
abline(h=30, col="black", lty=3)

#####
### MEAN SITE DEPTH FROM THE READ ALIGNMENTS
#####

## MISSING DATA AS A FUNCTION OF MEAN SITE DEPTH
mapped.files.depth <- mapped.files.list[[1]][,c(1,2,3,7)]
names(mapped.files.depth)[4] <- poolnames[1]

c = 5 
for (l in 2:length(mapped.files.list))
{
  mapped.files.depth[, c] <- mapped.files.list[[l]][, 7]
  names(mapped.files.depth)[c] <- poolnames[l]
  c = c + 1
}

mean_site_depth <- colMeans(mapped.files.depth[-c(831:1518), -c(1:3)])
mean_site_depth_data <- data.frame(MISSING=dt.repl.missing, SITE_DEPTH=mean_site_depth)

mean_site_depth_data['STATUS_COL'] <- rep(NA, dt.repl@npools)

mean_site_depth_data$STATUS_COL[mean_site_depth_data$MISSING > 30 & mean_site_depth_data$SITE_DEPTH < 1] <- "orange"
mean_site_depth_data$STATUS_COL[mean_site_depth_data$MISSING > 30 & mean_site_depth_data$SITE_DEPTH > 1] <- "red"
mean_site_depth_data$STATUS_COL[is.na(mean_site_depth_data$STATUS_COL)] <- "black"

COL_SITE_DEPTH <- adjustcolor(mean_site_depth_data$STATUS_COL, alpha.f = 0.5)

plot(MISSING ~SITE_DEPTH,
     xlab="Mean Site Depth",ylab="% of Missing Data", cex=1.5,
     col=COL_SITE_DEPTH, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_site_depth_data)
lines(lowess(mean_site_depth_data$MISSING ~ mean_site_depth_data$SITE_DEPTH, f=0.99), col='blue', lwd=2)
abline(v=1, col="black", lty=3)
abline(h=30, col="black", lty=3)

#####
### MEAN BASE QUALITY FROM THE READ ALIGNMENTS
#####

## MISSING DATA AS A FUNCTION OF MEAN BASE Q
mapped.files.baseq <- mapped.files.list[[1]][,c(1,2,3,8)]
names(mapped.files.baseq)[4] <- poolnames[1]

c = 5 
for (l in 2:length(mapped.files.list))
{
  mapped.files.baseq[, c] <- mapped.files.list[[l]][, 8]
  names(mapped.files.baseq)[c] <- poolnames[l]
  c = c + 1
}

mean_baseq <- colMeans(mapped.files.baseq[-c(831:1518), -c(1:3)])
mean_baseq_data <- data.frame(MISSING=dt.repl.missing, BASEQ=mean_baseq)

mean_baseq_data['STATUS_COL'] <- rep(NA, dt.repl@npools)

mean_baseq_data$STATUS_COL[mean_baseq_data$MISSING > 30] <- "orange"
mean_baseq_data$STATUS_COL[is.na(mean_baseq_data$STATUS_COL)] <- "black"

COL_BASEQ <- adjustcolor(mean_baseq_data$STATUS_COL, alpha.f = 0.5)

plot(MISSING ~BASEQ,
     xlab="Mean BaseQ",ylab="% of Missing Data", cex=1.5,
     col=COL_BASEQ, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_baseq_data)
lines(lowess(mean_baseq_data$MISSING ~ mean_baseq_data$BASEQ, f=0.99), col='blue', lwd=2)
abline(h=30, col="black", lty=3)

#####
### MEAN BASE MAPQ FROM THE READ ALIGNMENTS
#####

## MISSING DATA AS A FUNCTION OF MEAN BASE MAPQ
mapped.files.mapq <- mapped.files.list[[1]][,c(1,2,3,9)]
names(mapped.files.mapq)[4] <- poolnames[1]

c = 5 
for (l in 2:length(mapped.files.list))
{
  mapped.files.mapq[, c] <- mapped.files.list[[l]][, 8]
  names(mapped.files.mapq)[c] <- poolnames[l]
  c = c + 1
}

mean_mapq <- colMeans(mapped.files.mapq[-c(831:1518), -c(1:3)])
mean_mapq_data <- data.frame(MISSING=dt.repl.missing, MAPQ=mean_mapq)

mean_mapq_data['STATUS_COL'] <- rep(NA, dt.repl@npools)

mean_mapq_data$STATUS_COL[mean_mapq_data$MISSING > 30] <- "orange"
mean_mapq_data$STATUS_COL[is.na(mean_mapq_data$STATUS_COL)] <- "black"

COL_MAPQ <- adjustcolor(mean_mapq_data$STATUS_COL, alpha.f = 0.5)

plot(MISSING ~MAPQ,
     xlab="Mean MAPQ",ylab="% of Missing Data", cex=1.5,
     col=COL_MAPQ, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_mapq_data)
lines(lowess(mean_mapq_data$MISSING ~ mean_mapq_data$MAPQ, f=0.99), col='blue', lwd=2)
abline(h=30, col="black", lty=3)

## COMBINED PLOT
pdf("results/technical_replicates/minmaxcov_3_99/missing_read_mapping_stats.pdf",         # File name
    width = 8.50, height = 11, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(2, 2))
# PLOT 1
plot(MISSING ~SCAFFOLD_COV,
     xlab="% of Scaffold covered",ylab="% of Missing Data", cex=1.5,
     col=COL_SCAFFOLD_COV, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_scaffold_coverage_data)
lines(lowess(mean_scaffold_coverage_data$MISSING ~ mean_scaffold_coverage_data$SCAFFOLD_COV, f=0.99), col='blue', lwd=2)
abline(v=25, col="black", lty=3)
abline(h=30, col="black", lty=3)

# PLOT 2
plot(MISSING ~SITE_DEPTH,
     xlab="Mean Site Depth",ylab="% of Missing Data", cex=1.5,
     col=COL_SITE_DEPTH, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_site_depth_data)
lines(lowess(mean_site_depth_data$MISSING ~ mean_site_depth_data$SITE_DEPTH, f=0.99), col='blue', lwd=2)
abline(v=1, col="black", lty=3)
abline(h=30, col="black", lty=3)

# PLOT 3
plot(MISSING ~BASEQ,
     xlab="Mean BaseQ",ylab="% of Missing Data", cex=1.5,
     col=COL_BASEQ, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_baseq_data)
lines(lowess(mean_baseq_data$MISSING ~ mean_baseq_data$BASEQ, f=0.99), col='blue', lwd=2)
abline(h=30, col="black", lty=3)

# PLOT 4
plot(MISSING ~MAPQ,
     xlab="Mean MAPQ",ylab="% of Missing Data", cex=1.5,
     col=COL_MAPQ, pch=19, cex.lab=1.4, cex.axis=1.2, data=mean_mapq_data)
lines(lowess(mean_mapq_data$MISSING ~ mean_mapq_data$MAPQ, f=0.99), col='blue', lwd=2)
abline(h=30, col="black", lty=3)
dev.off()

## CONCLUSION: WHAT WERE THE CAUSES OF MISSING DATA?
## LOW SEQUENCE INPUT (for NW_BIO1_S1) AND LOW NUMBER q20 READS
## Pools from NW_BIO1_S1 and WO_BIO4_3 had less overall number of reads which caused:
## 1) Lower scaffold coverage: less than ~30% of the scaffold were coverade in these pools.
## 2) Lower depth: less than ~1 read/base were coverade in most of scaffolds of these pools.

## UNEVEN SEQUENCING
## Pools from PA_BIO1_S1 had a decent mean depth/based per scaffolds, decent of q20 mapping, but lower % of scaffold coverage
## which indicated that some  parts were more sequenced than other because probably genome shearing was not good.

## Savepoint_2
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
############## END THIS PART

###
###
### ---- MAXIMUM LIKELIHOOD ESTIMATION FROM AGGREGATED ALLELE COUNTS ----
###
###

### DATASET: AGGREGATED DATA OF EACH POOL

## Aggregate read counts data for each pool
## and repeat the allele frequency estimates
## and PCA.
coverage <- dt.repl@readcoverage
counts <- dt.repl@refallele.readcount

aggr.coverage <- cbind(coverage[,1]+coverage[,2]+coverage[,3]+coverage[,4]+coverage[,5], #MN_B1
                       coverage[,6]+coverage[,7]+coverage[,8]+coverage[,9],              #MN_B4
                       coverage[,10]+coverage[,11]+coverage[,12]+coverage[,13],          #ND_B1
                       coverage[,14]+coverage[,15]+coverage[,16]+coverage[,17],          #ND_B4
                       coverage[,18]+coverage[,19]+coverage[,20]+coverage[,21],          #NW_BIO1
                       coverage[,22]+coverage[,23]+coverage[,24]+coverage[,25],          #NW_BIO1
                       coverage[,26]+coverage[,27]+coverage[,28]+coverage[,29]+coverage[,30], #NW_BIO4
                       coverage[,31]+coverage[,32]+coverage[,33]+coverage[,34],               #NW_BIO4
                       coverage[,35]+coverage[,36]+coverage[,37]+coverage[,38],               #NW_BIO4
                       coverage[,39]+coverage[,40]+coverage[,41]+coverage[,42],               #PA_BIO1
                       coverage[,43]+coverage[,44]+coverage[,45]+coverage[,46],               #PA_BIO4
                       coverage[,47]+coverage[,48]+coverage[,49]+coverage[,50],               #PA_BIO4
                       coverage[,51]+coverage[,52]+coverage[,53]+coverage[,54],               #WI_BIO1
                       coverage[,55]+coverage[,56]+coverage[,57]+coverage[,58]+coverage[,59], #WI_BIO1
                       coverage[,60]+coverage[,61]+coverage[,62]+coverage[,63],               #WI_BIO4
                       coverage[,64]+coverage[,65]+coverage[,66]+coverage[,67],               #WI_BIO4
                       coverage[,68]+coverage[,69]+coverage[,70]+coverage[,71],               #WI_BIO4
                       coverage[,72]+coverage[,73]+coverage[,74]+coverage[,75],               #WO_BIO1
                       coverage[,76]+coverage[,77]+coverage[,78]+coverage[,79],               #WO_BIO4
                       coverage[,80]+coverage[,81]+coverage[,82]+coverage[,83],               #WO_BIO4
                       coverage[,84]+coverage[,85]+coverage[,86]+coverage[,87]               #WO_BIO4
)

aggr.counts <- cbind(counts[,1] +counts[,2] +counts[,3] +counts[,4] +counts[,5], #MN_B1
                     counts[,6] +counts[,7] +counts[,8] +counts[,9],              #MN_B4
                     counts[,10]+counts[,11]+counts[,12]+counts[,13],          #ND_B1
                     counts[,14]+counts[,15]+counts[,16]+counts[,17],          #ND_B4
                     counts[,18]+counts[,19]+counts[,20]+counts[,21],          #NW_BIO1
                     counts[,22]+counts[,23]+counts[,24]+counts[,25],          #NW_BIO1
                     counts[,26]+counts[,27]+counts[,28]+counts[,29]+counts[,30], #NW_BIO4
                     counts[,31]+counts[,32]+counts[,33]+counts[,34],               #NW_BIO4
                     counts[,35]+counts[,36]+counts[,37]+counts[,38],               #NW_BIO4
                     counts[,39]+counts[,40]+counts[,41]+counts[,42],               #PA_BIO1
                     counts[,43]+counts[,44]+counts[,45]+counts[,46],               #PA_BIO4
                     counts[,47]+counts[,48]+counts[,49]+counts[,50],               #PA_BIO4
                     counts[,51]+counts[,52]+counts[,53]+counts[,54],               #WI_BIO1
                     counts[,55]+counts[,56]+counts[,57]+counts[,58]+counts[,59], #WI_BIO1
                     counts[,60]+counts[,61]+counts[,62]+counts[,63],               #WI_BIO4
                     counts[,64]+counts[,65]+counts[,66]+counts[,67],               #WI_BIO4
                     counts[,68]+counts[,69]+counts[,70]+counts[,71],               #WI_BIO4
                     counts[,72]+counts[,73]+counts[,74]+counts[,75],               #WO_BIO1
                     counts[,76]+counts[,77]+counts[,78]+counts[,79],               #WO_BIO4
                     counts[,80]+counts[,81]+counts[,82]+counts[,83],               #WO_BIO4
                     counts[,84]+counts[,85]+counts[,86]+counts[,87]              #WO_BIO4
 )

nbr.gene.copies = dt.repl@poolsizes[1]
nbr.loci = dt.repl@nsnp
nbr.pops = dim(aggr.counts)[2]

# Receive the outputs
refcount.matrix = matrix(nrow=nbr.loci, ncol = nbr.pops)
altcount.matrix = matrix(nrow=nbr.loci, ncol = nbr.pops)
haplocount.matrix = matrix(nrow=nbr.loci, ncol = nbr.pops)

for (k in 1:nbr.pops)
{
  reads.allele.1 <- aggr.counts[,k]
  reads.allele.2 <- aggr.coverage[,k] - aggr.counts[,k]
  total.counts <- aggr.coverage[,k]
  
  allele.counts.1 <- vector(length = nbr.loci, mode = "integer")
  allele.counts.2 <- vector(length = nbr.loci, mode = "integer")
  haploid.counts  <- vector(length = nbr.loci, mode = "integer")
  
  likelihood <- matrix(nrow = nbr.loci, ncol = (nbr.gene.copies + 1))
  
  for (i in 1:nbr.loci) {
    if (total.counts[i] == 0){
      allele.counts.1[i] <- NA
      allele.counts.2[i] <- NA
      haploid.counts[i]  <- NA
      
    } else if ( total.counts[i] <= nbr.gene.copies) {
      allele.counts.1[i] <- reads.allele.1[i]
      allele.counts.2[i] <- reads.allele.2[i]
      haploid.counts[i] <- total.counts[i]  
      
    } else {
      if (reads.allele.1[i] == 0 | reads.allele.1[i] == total.counts[i]){
        allele.counts.1[i] <- reads.allele.1[i]
        allele.counts.2[i] <- reads.allele.2[i]
        haploid.counts[i] <- total.counts[i]  
        
      } else {
        for (j in 1:(nbr.gene.copies + 1)) 
        {
          likelihood[i,j] <- dbinom(x = reads.allele.1[i], size = total.counts[i],
                                    prob = (j - 1) / nbr.gene.copies, log = FALSE)
        } # end of for i,j
        
        allele.counts.1[i] <- (which(likelihood[i,] == max(likelihood[i,])) - 1)
        allele.counts.2[i] <- (nbr.gene.copies - allele.counts.1[i])
        haploid.counts[i] <- nbr.gene.copies
      } # end of inside if
      
    } # end of if
    
  } # end of for i
  
  refcount.matrix[,k] <- allele.counts.1
  altcount.matrix[,k] <- allele.counts.2
  haplocount.matrix[,k] <- haploid.counts
  
} # end of FOR k

out_list <- list(ref_count = refcount.matrix, 
                 alt_count = altcount.matrix,
                 hap_count = haplocount.matrix)


dt.aggregated.imputedRefMLCount <- out_list
dt.aggregated.imputedRefMLFreq <- dt.aggregated.imputedRefMLCount[[1]]/dt.aggregated.imputedRefMLCount[[3]]


###
###
### ---- PCA ON ML ALLELE FREQUENCIES FOR AGGREGATED TECHNICAL REPLICATES READS ----
###
###

legend.colors.aggr <- c("#0000FF", "#FF0000", #MN
                       "#8282FF", "#FF8282", #ND
                       "#10B717", "#64C366", "#FD63CB", "#F98AD1", "#F4A8D6", #NW
                       "#B78560", "#9E5C00", "#7D3200", #PA
                       "#B2B2FF", "#D5D5FF", "#FFB2B2", "#FFD5D5", "#FFF2F2", #WI
                       "#90CE91", "#EFC0DB", "#E9D2DF", "#E4DFE2" #WO
                     )

legend.names.aggr <- c("MN-B1.1", "MN-B4.1",
                       "ND-B1.1", "ND-B4.1", 
                       "NW-B1.1", "NW-B1.2", "NW-B4.1", "NW-B4.2", "NW-B4.3", 
                       "PA-B1.1", "PA-B4.1", "PA-B4.2", 
                       "WI-B1.1", "WI-B1.2", "WI-B4.1", "WI-B4.2", "WI-B4.3", 
                       "WO-B1.1", "WO-B4.1", "WO-B4.2", "WO-B4.3"
                    )

## EXPORT THE IMPUTED ML ALLELE FREQUENCY MATRIX TO A FILE
## WITH THE SNP INFORMATION AND POOL NAMES (FOR THE COMPARATIVE PCA ANALYSIS)
dt.aggregated.imputedRefMLFreq.data <- dt.aggregated.imputedRefMLFreq
colnames(dt.aggregated.imputedRefMLFreq.data) <- legend.names.aggr

dt.aggregated.imputedRefMLFreq.dataFrame <- data.frame(dt.repl@snp.info, dt.aggregated.imputedRefMLFreq.data)

# Save to a file
#write.table(dt.aggregated.imputedRefMLFreq.dataFrame,
#            file="results/technical_replicates/minmaxcov_3_99/imputedML_REFAlleleFrequencies_repl-aggregated_pools.txt", sep = "\t",
#            quote = F, row.names = F, col.names = T)


## PRINCIPAL COMPONENT ANALYSIS
W<-scale(t(dt.aggregated.imputedRefMLFreq), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")

l1 <- 100*lambda[1]/sum(lambda) 
l2 <- 100*lambda[2]/sum(lambda)

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=legend.colors.aggr,
     xlim = c(-0.4, 0.5), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6
     )
legend("topleft", 
       legend = legend.names.aggr, col = legend.colors.aggr, 
       pch = legend.pch, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)


## Without the outlier population
W<-scale(t(dt.aggregated.imputedRefMLFreq[,-c(5,21)]), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")

l1 <- 100*lambda[1]/sum(lambda) 
l2 <- 100*lambda[2]/sum(lambda)

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=legend.colors.aggr[-c(5,21)],
     #xlim = c(-0.75, 0.45), ylim = c(-0.55, 0.65),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
legend("bottomright", 
       legend = legend.names.aggr[-c(5,21)], col = legend.colors.aggr[-c(5,21)], 
       pch = legend.pch[-c(5,21)], bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)


## PCA PLOT BY BIOTYPE
biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53")

pdf("results/technical_replicates/minmaxcov_3_99/pca_imputed_allelefreq_repl_aggregated_biotypecolors.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=biotype.col[-c(5,21)],
     xlim = c(-0.5, 0.75), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec[,1],eig.vec[,2], 
     legend.names.aggr[-c(5,21)], pos=2 , cex = 0.6)
legend("bottomright", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## PCA with FactorMineR - EXPERIMENTAL

biotype.sym <- c(15,17,
                 15,17,
                 15,15,17,17,17,
                 15,17,17,
                 15,15,17,17,17,
                 15,17,17,17)

# Convert NA cells into loci means (without outlier populations)
dt.aggregated.imputedRefMLFreq.NAimputed = na.aggregate(t(dt.aggregated.imputedRefMLFreq[,-c(5,21)]))

# PCA
pca.aggregated <- PCA(dt.aggregated.imputedRefMLFreq.NAimputed, scale.unit = F, graph = F)

barplot(pca.aggregated$eig[,1],main="Eigenvalues",names.arg=1:nrow(pca.aggregated$eig))
summary(pca.aggregated)

## Run cluster analysis

# 1. Loading and preparing data
pca.aggregated.ind_coord <- pca.aggregated$ind$coord

cluster_discovery = fviz_nbclust(pca.aggregated.ind_coord, kmeans, method = "gap_stat")
#ggsave("cluster_discovery.pdf",cluster_discovery,  width =4, height = 4)

# 2. Compute k-means at K=2, K=3 and at K=4
set.seed(123)
kmeans_k2 <- kmeans(pca.aggregated.ind_coord, 2, nstart = 20)
kmeans_k3 <- kmeans(pca.aggregated.ind_coord, 3, nstart = 20)
kmeans_k4 <- kmeans(pca.aggregated.ind_coord, 4, nstart = 20)
kmeans_k5 <- kmeans(pca.aggregated.ind_coord, 5, nstart = 20)

# Combine Pool ID with cluster assigment
data.frame(k2_clusters = kmeans_k2$cluster,
           k3_clusters = kmeans_k3$cluster,
           k4_clusters = kmeans_k4$cluster,
           k5_clusters = kmeans_k5$cluster) %>% mutate(sampleId = legend.names.aggr[-c(5,21)]) -> pca.aggregated_cluster_data


###COMBINED PLOTS
pdf("results/technical_replicates/minmaxcov_3_99/pca_imputed_allelefreq_repl_aggregated_factorminer_clusters.pdf",  
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(2, 2))
## PLOT PCA COLORED BY BIOTYPE
plot(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], col=biotype.col[-c(5,21)],
     xlim = c(-2, 2), ylim = c(-3, 2),
     xlab=paste0("PC1 (", round(pca.aggregated$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.aggregated$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], 
     pca.aggregated_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K2
plot(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], col=pca.aggregated_cluster_data$k2_clusters,
     xlim = c(-2, 2), ylim = c(-3, 2),
     xlab=paste0("PC1 (", round(pca.aggregated$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.aggregated$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], 
     pca.aggregated_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2"), col = c(1,2), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K3
plot(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], col=pca.aggregated_cluster_data$k3_clusters,
     xlim = c(-2, 2), ylim = c(-3, 2),
     xlab=paste0("PC1 (", round(pca.aggregated$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.aggregated$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], 
     pca.aggregated_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2", "K3"), col = c(1,2,3), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K4
plot(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], col=pca.aggregated_cluster_data$k4_clusters,
     xlim = c(-2, 2), ylim = c(-3, 2),
     xlab=paste0("PC1 (", round(pca.aggregated$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.aggregated$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.aggregated$svd$U[,1], pca.aggregated$svd$U[,2], 
     pca.aggregated_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2", "K3", "K4"), col = c(1,2,3,4), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Savepoint_3
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
############## END THIS PART

#####
### PCA with prcomp - EXPERIMENTAL
## First try to impute missing data
#dt.aggregated.imputedRefMLFreq_nafilled <- imputePCA(dt.aggregated.imputedRefMLFreq[,-c(5,21)], method="EM",ncp=1, seed=4401)
#
## PCA
#biotypes.pr <- prcomp(t(dt.aggregated.imputedRefMLFreq_nafilled$fittedX), center = TRUE, scale = TRUE)
#summary(biotypes.pr)
#
#screeplot(biotypes.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
#abline(h = 1, col="red", lty=5)
#legend("topright", legend=c("Eigenvalue = 1"),
#       col=c("red"), lty=5, cex=0.6)
#
#cumpro <- cumsum(biotypes.pr$sdev^2 / sum(biotypes.pr$sdev^2))
#plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
#abline(v = 2, col="blue", lty=5)
#abline(h = 1, col="blue", lty=5)
#legend("topleft", legend=c("Cut-off @ PC6"),
#       col=c("blue"), lty=5, cex=0.6)
#
#plot(biotypes.pr$x[,1], 
#     biotypes.pr$x[,2], 
#     xlab="PC1 (81.6%)", ylab = "PC2 (18.3%)", main = "PC1 / PC2 - plot")
#
#biotype.vec <- c("B1","B4","B1","B4","B1","B1","B4","B4","B4","B1",
#                 "B4","B4","B1","B1","B4","B4","B4","B1","B4","B4",
#                 "B4")
## final plot
#fviz_pca_ind(biotypes.pr, geom.ind = "point", pointshape = 21, 
#             pointsize = 2, 
#             fill.ind = biotype.vec[-c(5,21)], 
#             col.ind = "black", 
#             palette = "jco", 
#             addEllipses = TRUE,
#             label = "var",
#             col.var = "black",
#             repel = TRUE,
#             legend.title = "Biotypes") +
#  ggtitle("2D PCA-plot from 30 feature dataset") +
#  theme(plot.title = element_text(hjust = 0.5))
#####


###
###
### ---- % OF MISSING DATA AND CORRELATION WITH REPLICATE-AGGREGATED COVERAGE----
###
###

# NA's count
colnames(refcount.matrix) <- legend.names.aggr
dt.missing.aggr <- 100*(colSums(is.na(refcount.matrix))/nbr.loci)

barplot(dt.missing.aggr)
abline(h=30, col="red")

# Mean replicate coverage
mean.coverage.aggregated <- apply(aggr.coverage, 2, mean)
mean(mean.coverage.aggregated) +c(-1.96, +1.96)*sd(mean.coverage.aggregated)/sqrt(21)

## MISSING DATA AS A FUNCTION OF COVERAGE
plot(dt.missing.aggr ~mean.coverage.aggregated,
     xlab="Mean coverage (number of reads)",ylab="% of Missing Data", cex=1.5,
     col="black", pch=19, cex.lab=1.4, cex.axis=1.2)
lines(lowess(dt.missing.aggr ~mean.coverage.aggregated, f=0.99), col='blue', lwd=2)

quantile(aggr.coverage, probs = c(0.025, 0.5, 0.975))
hist(aggr.coverage, breaks = 100)
abline(v= 212, col="red")

# FUNCTION TO PRODUCE THE PLOTS - LOG and LOGIT SCALE
plot.mult.densities <- function(table=NULL, par.name=NULL, 
                                col.vect=NULL, name.vect=NULL, y_lim=c(0,0.2), cex_axis = 1.2, cex_lab  = 1.2,
                                plot.legend=TRUE, legend.side="topright")
{
  xlim.min = min(apply(table, 2, min))
  xlim.max = min(apply(table, 2, max))
  
  plot(density(table[,1]), lty=3, col=col.vect[1], xlim=c(xlim.min, xlim.max), ylim=y_lim, 
       main = "", ylab = "Density", xlab = par.name, cex.axis = cex_axis, cex.lab  = cex_lab)
  lines(density(table[,1]), col=col.vect[1])
  
  for (p in 2:dim(table)[2])
  {
    lines(density(table[,p]),lty=3, col=col.vect[p])
    lines(density(table[,p]), col=col.vect[p])
  }
  if (plot.legend) {
    legend(legend.side, col = col.vect, lty = 1, cex = 0.7,
           legend = name.vect, box.lwd = 0, box.col = NULL, bg = NULL)
  }
}

plot.mult.densities(aggr.coverage, col.vect = legend.colors.aggr, name.vect = legend.names.aggr)

## Savepoint_4
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
############## END THIS PART

###
###
### ---- SNP INFORMATION: CHRM, POS, MEAN_COV, MIN/MAX_COV, % MISSING DATA ----
###
###

# copy of read coverage
copy_aggr.coverage <- aggr.coverage
copy_aggr.coverage[copy_aggr.coverage==0] <- NA

# SNP COVERAGE
snp_mean_cov <- apply(copy_aggr.coverage, 1, function(x) mean(x, na.rm = T))
snp_min_cov <- apply(copy_aggr.coverage, 1, function(x) min(x, na.rm = T))
snp_max_cov <- apply(copy_aggr.coverage, 1, function(x) max(x, na.rm = T))

# SNP REFERENCE FREQUENCY
snp_mean_refMLFreq <- apply(dt.aggregated.imputedRefMLFreq, 1, function(x) mean(x, na.rm = T))
snp_min_refMLFreq <- apply(dt.aggregated.imputedRefMLFreq, 1, function(x) min(x, na.rm = T))
snp_max_refMLFreq <- apply(dt.aggregated.imputedRefMLFreq, 1, function(x) max(x, na.rm = T))

# SNP % MISSING DATA
snp_mean_missing <- apply(dt.aggregated.imputedRefMLFreq, 1, function(x) sum(is.na(x))/21)*100

snp_info_techaggregated <- data.frame(dt.repl@snp.info, 
                                      MEAN_COV=snp_mean_cov, MIN_COV=snp_min_cov, MAX_COV=snp_max_cov, 
                                      MEAN_REFFREQ=snp_mean_refMLFreq, MIN_REFFREQ=snp_min_refMLFreq, MAX_REFFREQ=snp_max_refMLFreq,
                                      MISSING=snp_mean_missing)

# Save to a file all the information
write.table(snp_info_techaggregated,
            file="results/technical_replicates/minmaxcov_3_99/snpinfo_replicates_aggregated_data.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

# Save to a file only the combined chromosome and SNP position
write.table(paste0(snp_info_techaggregated$Chromosome, "_", snp_info_techaggregated$Position),
            file="results/technical_replicates/minmaxcov_3_99/snpinfo_replicates_aggregated_names.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)

## Savepoint_5
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/technical_replicates/minmaxcov_3_99/check.replicates.pca.workspace27Jul21.RData")
############## END THIS PART
