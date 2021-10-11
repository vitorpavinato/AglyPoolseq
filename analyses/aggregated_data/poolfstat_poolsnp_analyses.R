##########################################################
#              POOLSNP-POOLFSTAT Analyses                # 
#         AF, HE, Theta, pF_ST, gF_ST and genome scan    #
##########################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

###
###
### ---- SETUP ANALYZES ENVIRONMENT ----
###
###

### Recover R-renv environment
setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()

### Remove last features
rm(list=ls())
ls()

### Load R libraries
library(testthat)
library(poolfstat)
#packageVersion("poolfstat") 
library(tidyverse)
library(missMDA)
library(ade4)
library(factoextra)
library(PopGenReport)
#library(vegan)
#library(ggplot2)
#library(viridis)

### Import auxiliary R functions  
source("DEST-AglyPoolseq/analyses/aggregated_data/aux_func.R")

### SETUP PYTHON INTEGRATION
library(reticulate) # To use python code within R
use_python("~/.conda/envs/local_py3/bin/python")

### Import Python packages
# Python installation: anaconda python3
sns <- import('seaborn')
plt <- import('matplotlib.pyplot')
pd <- import('pandas')

### POOL'S INFORMATION
# Sample's size
poolsizes <- rep(10,21)

# Sample's names
poolnames <- c("MN-Av.1", "MN-V.1",
               "ND-Av.1", "ND-V.1", 
               "NW-Av.1", "NW-Av.2", "NW-V.1", "NW-V.2", "NW-V.3", 
               "PA-Av.1", "PA-V.1", "PA-V.2", 
               "WI-Av.1", "WI-Av.2", "WI-V.1", "WI-V.2", "WI-V.3", 
               "WO-Av.1", "WO-V.1", "WO-V.2", "WO-V.3")

# Sample's colors - for PCA
biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53")


###
###
### ---- POPULATION GENOMICS ANALYZES ----
###
###

### Dataset 2: 24-Jun-2021; Min_cov=4; Max_cov=0.99; MAF=0.05; MAC=5; 21 pools; miss_fraction=0.50
### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz

### LOAD VCF FILE
dt.1 <- vcf2pooldata(
                     vcf.file = "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz",
                     poolsizes = poolsizes,
                     poolnames = poolnames,
                     min.cov.per.pool = -1,
                     min.rc = 1,
                     max.cov.per.pool = 1e+06,
                     min.maf = 0.05,
                     remove.indels = FALSE,
                     nlines.per.readblock = 1e+06
)
dt.1
# * * * PoolData Object * * * 
# * Number of SNPs   = 262866
# * Number of Pools  = 21 
# Number of scaffolds:

length(unique(sort(dt.1@snp.info$Chromosome))) #705

### EXPORT DATASET SNPS INFORMATION AS A SORTED BED FILE -LIKE
### AFTER READ VCF WITH POOLFSTAT. POOLFSTAT FILTER OUT SOME SNPS
### (SNPS WITH MORE THAN 2 ALLELES). THE LIST OF SNPS IS ALSO
### USEFUL TO IDENTIFY THE VARIATIONS CLOSE TO MULTILOCUS
### FST VALUES THAT ARE SIGNIFICANT IN THE GENOME SCAN AS
### THE MULTILOCUS FST SCAN ONLY REPORTS THE WINDOW.

# Get all SNPs information
snps_ordered <- dt.1@snp.info[order(dt.1@snp.info$Chromosome, dt.1@snp.info$Position), -c(3,4)]

# Create a BED file -like table
snps_ordered.bed <- data.frame(CHROM=snps_ordered$Chromosome,
                               Start=as.integer(snps_ordered$Position-1),
                               Ends=as.integer(snps_ordered$Position))
# Save to a file
#write.table(snps_ordered.bed,
#            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/snps_ordered.bed", sep = "\t",
#            quote = F, row.names = F, col.names = F)
############## END THIS PART

###
###
### ---- MAXIMUM LIKELIHOOD ESTIMATION OF THE ALLELE COUNTS ----
###
###

### COMPUTE MAXIMUM LIKELIHOOD IMPUTED SAMPLE ALLELE COUNT
dt.1.imputedRefMLCount <- imputedRefMLCount(dt.1)
dt.1.imputedRefMLFreq <- dt.1.imputedRefMLCount[[1]]/dt.1.imputedRefMLCount[[3]]

### EXPORT THE IMPUTED ML ALLELE FREQUENCY MATRIX TO A FILE
### WITH THE SNP INFORMATION AND POOL NAMES (FOR THE PCA ANALYSIS)
dt.1.imputedRefMLFreq.data <- dt.1.imputedRefMLFreq
colnames(dt.1.imputedRefMLFreq.data) <- dt.1@poolnames

dt.1.imputedRefMLFreq.dataFrame <- data.frame(dt.1@snp.info, dt.1.imputedRefMLFreq.data)

# Save to a file
#write.table(dt.1.imputedRefMLFreq.dataFrame,
#            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_pools.txt", sep = "\t",
#            quote = F, row.names = F, col.names = T)

### SANITY CHECK MAF AND MAC VALUES AFTER POOLFSTAT LOADING

# CHECK RAW MAF
dt.1.overall.ALTallele.freq <- apply(dt.1.imputedRefMLFreq, 1,  function(x) 1-mean(x, na.rm = T))
dt.1.overall.REFallele.freq <- apply(dt.1.imputedRefMLFreq, 1,  function(x) mean(x, na.rm = T))

dt.1.MAF <- pmin(dt.1.overall.REFallele.freq, dt.1.overall.ALTallele.freq)
min(dt.1.MAF); max(dt.1.MAF)
hist(dt.1.MAF)

# CHECK RAW MAC
dt.1.ALTallele.count <- (dt.1@readcoverage - dt.1@refallele.readcount)

dt.1.overall.ALTallele.count <- apply(dt.1.ALTallele.count, 1, sum)
dt.1.overall.REFallele.count <- apply(dt.1@refallele.readcount, 1, sum)

dt.1.MAC <- pmin(dt.1.overall.REFallele.count, dt.1.overall.ALTallele.count)
min(dt.1.MAC); max(dt.1.MAC)
hist(dt.1.MAC)

### SANITY CHECK THE MAF SPECTRUM DISTRIBUTION - LOOK FOR THE U-SHAPE SFS - ADMIXTURE
#hist(pmin(dt.1@refallele.readcount[,1]/dt.1@readcoverage[,1], (1-(dt.1@refallele.readcount[,1]/dt.1@readcoverage[,1]))))
#hist((1-(dt.1@refallele.readcount[,1]/dt.1@readcoverage[,1])))
#hist(1-dt.1.imputedRefMLFreq[,21], breaks = 10)
#
#hist(pmin(dt.1.imputedRefMLCount[[1]][,21], (dt.1.imputedRefMLCount[[2]][,21])), breaks = 7)
#hist(dt.1.imputedRefMLCount[[2]][,21], breaks = 10)

### CHECK HOW MANY SNPS IS FIXED IN 95% OF POOLS
hist(apply(dt.1.imputedRefMLFreq == 1, 1, function(x) sum(x, na.rm = T)))
table(apply(dt.1.imputedRefMLFreq == 1, 1, function(x) sum(x, na.rm = T)))
quasi.fixed <- which(apply(dt.1.imputedRefMLFreq == 1, 1, function(x) sum(x, na.rm = T)) == 20) 

# Check the MAF for almost fixed markers
#dt.1.fixed.ALTallele.freq <- apply(dt.1.imputedRefMLFreq[quasi.fixed ,], 1,  function(x) 1-mean(x, na.rm = T))
#dt.1.fixed.REFallele.freq <- apply(dt.1.imputedRefMLFreq[quasi.fixed ,], 1,  function(x) mean(x, na.rm = T))
#dt.1.MAF.fixed <- pmin(dt.1.fixed.REFallele.freq, dt.1.fixed.ALTallele.freq)
#min(dt.1.MAF.fixed); max(dt.1.MAF.fixed)

### CHECK THE % OF MISSING DATA
# POOLS
dt.1.missing.pools <- (apply(dt.1.imputedRefMLFreq, 2, function(x) sum(is.na(x)))/dt.1@nsnp)*100
names(dt.1.missing.pools) <- poolnames

barplot(dt.1.missing.pools, ylim = c(0,100))
abline(h=30, col="red")

# SNPs
dt.1.missing.snps <- (apply(dt.1.imputedRefMLFreq, 1, function(x) sum(is.na(x)))/dt.1@npools)*100
hist(dt.1.missing.snps, breaks = 10)

#### REMOVE 2 SAMPLES IDENTIFIED WITH > 30% MISSING DATA (DONT REMOVE PA SAMPLES)
#reduced.dt.1.imputedRefMLFreq <- dt.1.imputedRefMLFreq[,-c(5,21)]
## Check Missing %
#dt.1.missing.snps.reduced <- (apply(reduced.dt.1.imputedRefMLFreq, 1, function(x) sum(is.na(x)))/(dt.1@npools-2))*100
#hist(dt.1.missing.snps.reduced, breaks = 10)

## Savepoint_1
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
############## END THIS PART

###
###
### ---- RUN MULTILOCUS-FST GENOME SCAN ----
###
###

### MULTI-LOCUS FST GENOME SCAN
multilocusfst_window_size = 10
dt.1.fst.scan <- computeFST(dt.1, sliding.window.size = multilocusfst_window_size)

### Multilocus FST thresholds
# 99% quantile
dt.1.fst.scan.thresh <-  quantile(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst, probs = 0.99, na.rm = T)

# Arbritrary FST > 0.25
dt.1.fst.scan.025.thres <- 0.25

### MANHATTAN PLOT
# Color scaffolds
chrm_colors_multifst <- colorChromosomes(x=dt.1.fst.scan$sliding.windows.fst$Chr)

# SAVE MANHATTAN PLOT TO A PDF FILE
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/genome_scan_multilocusFST.pdf", # File name
    width = 11, height = 8.50,                                                                # Width and height in inches
    bg = "white",                                                                             # Background color
    colormodel = "cmyk",                                                                      # Color model (cmyk is required for most publications)
    )

par(mar=c(5,5,4,1)+.1)
plot(dt.1.fst.scan$sliding.windows.fst$CumulatedPosition/1e6,
     dt.1.fst.scan$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab=expression(paste("Muli-locus ", italic(F)[ST])),
     col=chrm_colors_multifst,pch=20, cex.lab=1.4, cex.axis=1.2)

# Add thresholds
abline(h=c(dt.1.fst.scan$FST, dt.1.fst.scan.thresh, dt.1.fst.scan.025.thres),
       lty=c(2,2,2), col=c("lightgreen", "blue", "red"), lwd=2)

# Add legends
legend("topright", 
       legend = c(expression(paste("Global ", italic(F)[ST])), 
                              expression(paste(italic(F)[ST], " > 99% quantile")),
                              expression(paste(italic(F)[ST], " > 0.25"))),
       lty = 2, lwd=2, col=c("lightgreen", "blue", "red"), bty="n")

#text(x=15, y=c(dt.1.fst.scan$FST, dt.1.fst.scan.thresh, dt.1.fst.scan.025.thres)-0.02,
#     c(expression(paste("Global ", italic(F)[ST])), 
#       expression(paste(italic(F)[ST], " > 99% quantile")),
#       expression(paste(italic(F)[ST], " > 0.25"))),
#     cex=0.75, pos=3,col=c("lightgreen", "blue", "red"))

dev.off()

### ISOLATE THE MULTILOCUS FST > 0.25
high.multilocus.fst <- data.frame(CHRM=dt.1.fst.scan$sliding.windows.fst$Chr[which(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)], 
                                  POS=dt.1.fst.scan$sliding.windows.fst$Position[which(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)], 
                                  FST= dt.1.fst.scan$sliding.windows.fst$MultiLocusFst[which(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)])

# Order by scaffold and position
high.multilocus.fst.ordered <- high.multilocus.fst[order(high.multilocus.fst$CHRM, high.multilocus.fst$POS), ]

# Save to a file
#write.table(high.multilocus.fst.ordered, file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst_values.txt",
#            quote = F, col.names = T, row.names = F)

# Create a BED file -like table
high.multilocus.fst.ordered.bed <- data.frame(CHROM=high.multilocus.fst.ordered$CHRM,
                                              Start=high.multilocus.fst.ordered$POS-1,
                                              Ends=high.multilocus.fst.ordered$POS)
# Save to a file
#write.table(high.multilocus.fst.ordered.bed,
#            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst.bed", sep = "\t",
#            quote = F, row.names = F, col.names = F)

### MANHATTAN PLOT OF SIGNIFICANT MULTILOCUS FST SCAFFOLDS

# Isolate the significant scaffolds
intralocus_fst <- data.frame(dt.1@snp.info, fst=dt.1.fst.scan$snp.FST)
intralocus_fst_sign <- intralocus_fst[intralocus_fst$Chromosome %in% unique(high.multilocus.fst.ordered$CHRM), ]
#table(intralocus_fst_sign$Chromosome)

# Multi-manhattam plots for all significant scaffolds
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/significant_intralocus_fst_scaffolds.pdf",  # File name
    width = 8.5, height = 11,                                                                             # Width and height in inches
    bg = "white",                                                                                         # Background color
    colormodel = "cmyk",                                                                                  # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(4, 4)) # 16 scaffolds
for (s in seq_along(unique(intralocus_fst_sign$Chromosome)))
{
  chr_name <- unique(intralocus_fst_sign$Chromosome)[s]
  n <- intralocus_fst_sign[intralocus_fst_sign$Chromosome %in% chr_name,]
  
  plot(n$fst ~ n$Position, 
       xlab="Position (in bp)",ylab=expression(paste(italic(F)[ST])), main=chr_name,
       col="black",pch=20, cex.lab=1.4, cex.axis=1.2)
  abline(h=c(dt.1.fst.scan$FST, dt.1.fst.scan.thresh, dt.1.fst.scan.025.thres),
         lty=c(2,2,2), col=c("lightgreen", "blue", "red"), lwd=2)
}
dev.off()

# EXPORT ALL INTRA-LOCUS FST VALUES
#write.table(intralocus_fst, file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/intralocus_fst_data.txt",
#            quote = F, col.names = T, row.names = F)

### SANITY CHECK HOW MISSING DATA AFFECTED THE NUMBER OF CALCULATED INTRA-LOCUS FSTS
# CHECK THE NUMBER OF SNPS WITHOUT NaN RESULTS
sum(!is.na(dt.1.fst.scan$snp.FST)) #19842
sum(!is.na(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst)) #18207

### CHECK THE NUMBER OF SNPS WITH NO MISSING DATA
# Data set with all 21 pools
no.missing.all <- dt.1.imputedRefMLFreq[complete.cases(dt.1.imputedRefMLFreq), ]
dim(no.missing.all) #19842    21

# Data set with 19 pools
#no.missing.reduced <- reduced.dt.1.imputedRefMLFreq[complete.cases(reduced.dt.1.imputedRefMLFreq), ]
#dim(no.missing.reduced) #33313    19

#Data set with 16 pools (no PA samples) # IDEAL SCENARIO
no.missing.reduced.2 <- dt.1.imputedRefMLFreq[complete.cases(dt.1.imputedRefMLFreq[,-c(5,10:12,21)]) , ]
dim(no.missing.reduced.2) #82534    21

## Savepoint_2
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
############## END THIS PART

###
###
### ---- EXPORT DATA TO BAYPASS ----
###
###

## pooldata2genobaypass
#pooldata2genobaypass(pooldata = dt.1,
#                     writing.dir = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/",
#                     prefix = "aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.complete50")
#
############## END THIS PART

###
###
### ---- DIVERSTITY AND COVERAGE STATS FOR EACH SNP ASSOCIATED WITH SIGNIFICANT SCAFFOLDS ----
###
###

### OVERALL HETEROZYGOSITY ACROSS SAMPLES - TOTAL HETEROZYGOSITY H_T

## GENOME-WIDE MEAN HE FOR EACH SAMPLE
# OVERALL
locus_totalhe_pools       <- apply(dt.1.imputedRefMLFreq, 1, totalHE) # 1 indicates rows: locus

# ACROSS SAMPLES/BIOTYPES
locus_totalhe_B1_pools    <- apply(dt.1.imputedRefMLFreq[,-c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, totalHE)
locus_totalhe_B4_pools    <- apply(dt.1.imputedRefMLFreq[, c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, totalHE)
# NOTE: FOR THE MOMENT USE THIS, BUT IDEALLY CALCULATE MEAN LOCUS-HE ACROSS WINDOWS WITH BIOTYPES

### AVERAGE PER-SAMPLE DEPTH OF BASES
dt.1.read_coverage_data <- dt.1@readcoverage
dt.1.read_coverage_data[ dt.1.read_coverage_data==0 ] <- NA

locus_ADP_pools <- apply(dt.1.read_coverage_data, 1, function(x) mean(x, na.rm = T))
locus_ADP_B1_pools <- apply(dt.1.read_coverage_data[,-c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) mean(x, na.rm = T))
locus_ADP_B4_pools <- apply(dt.1.read_coverage_data[, c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) mean(x, na.rm = T))

# Preapare dataFrame with significant FST scan scaffolds, intra-locus HEs and mean DEPTH
intralocus_he_adp     <- data.frame(dt.1@snp.info, 
                                    he_pools=locus_totalhe_pools, he_b1=locus_totalhe_B1_pools, he_b4=locus_totalhe_B4_pools,
                                    adp_pools=locus_ADP_pools, adp_b1=locus_ADP_B1_pools, adp_b4=locus_ADP_B4_pools)

intralocus_he_adp_signfst <- intralocus_he_adp[intralocus_he_adp$Chromosome %in% unique(high.multilocus.fst.ordered$CHRM), ]
#table(intralocus_he_adp_signfst$Chromosome)

### PLOT Multi LOCUS-HE plots for all significant scaffolds
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/significant_intralocus_fst_he_scaffolds.pdf",         # File name
    width = 8.5, height = 11, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(4, 4)) # 16 chromosomes
for (s in seq_along(unique(intralocus_he_adp_signfst$Chromosome)))
{
  chr_name <- unique(intralocus_he_adp_signfst$Chromosome)[s]
  n <- intralocus_he_adp_signfst[intralocus_he_adp_signfst$Chromosome %in% chr_name,]
  
  plot(n$he_b1 ~ n$Position, type="n", 
       xlab="Position (in bp)",ylab=expression(paste(italic(H)[E])), main=chr_name,
       ylim=c(0,0.6),
       col="black",pch=20, cex.lab=1.4, cex.axis=1.2)
  lines(lowess(n$he_pools ~ n$Position, f=0.10), col="darkgray", lwd=2, lty=2)
  lines(lowess(n$he_b1 ~ n$Position, f=0.10), col="#247F00", lwd=2)
  lines(lowess(n$he_b4 ~ n$Position, f=0.10), col="#AB1A53", lwd=2)
  abline(h=mean(locus_totalhe_pools), lty=3, col="black", lwd=1)
}
# Add legend on the last plot
legend("bottomleft", legend = c("Combined", "Avirulent", "Virulent", expression(paste("Average ", italic(H)[E], ))), 
       lty = c(2,1,1,3), lwd=c(2,2,2,1), col=c("darkgray", "#247F00", "#AB1A53", "black"), bty="n")
dev.off()

### PLOT Multi AVERAGE PER-SAMPLE DEPTH OF BASES plots for all significant scaffolds
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/significant_intralocus_fst_adp_scaffolds.pdf", # File name
    width = 8.5, height = 11,                                                                                # Width and height in inches
    bg = "white",                                                                                            # Background color
    colormodel = "cmyk",                                                                                     # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(4, 4)) # 16 chromosomes
for (s in seq_along(unique(intralocus_he_adp_signfst$Chromosome)))
{
  chr_name <- unique(intralocus_he_adp_signfst$Chromosome)[s]
  n <- intralocus_he_adp_signfst[intralocus_he_adp_signfst$Chromosome %in% chr_name,]
  
  plot(n$adp_b1 ~ n$Position, type="n", ylim=c(4,40),
       xlab="Position (in bp)",ylab="ADP", main=chr_name,
       col="black",pch=20, cex.lab=1.4, cex.axis=1.2)
  lines(lowess(n$adp_pools ~ n$Position, f=0.10), col="darkgray", lwd=2, lty=2)
  lines(lowess(n$adp_b1 ~ n$Position, f=0.10), col="#247F00", lwd=2)
  lines(lowess(n$adp_b4 ~ n$Position, f=0.10), col="#AB1A53", lwd=2)
  abline(h=4, lty=3, col="black", lwd=1)
}
# Add legend on the last plot
legend("topleft", legend = c("Combined", "Avirulent", "Virulent", "Min read cov"), 
       lty = c(2,1,1,3), lwd=c(2,2,2,1), col=c("darkgray", "#247F00", "#AB1A53", "black"), bty="n")
dev.off()

### ALLELE FREQUENCY CORRELATION ACROSS THE SIGNIFICAN SCAFFOLDS
## WITH THE LINE PLOT OF ALLELE FREQ CORRELATIONSI WANT TO ANSWER THESE QUESTIONS:
## THEY ARE POSITIVELY CORRELATED ACROSS THE SIGNIFICAN SCAFFOLDS? OR
## THE AF IS NEGATIVELY CORRELATION (AF GOES IN OPPOSITE DIRECTIONS)?

### BIOTYPE-AVERAGES reference allele frequency from the references allele frequencies 
locus_RefFREQ_B1_pools <- apply(dt.1.imputedRefMLFreq.data[,-c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) mean(x, na.rm = T))
locus_RefFREQ_B4_pools <- apply(dt.1.imputedRefMLFreq.data[, c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) mean(x, na.rm = T))

intralocus_RefFREQ     <- data.frame(CHROM = dt.1@snp.info[,1], POS = dt.1@snp.info[,2],
                                     AF_1 = locus_RefFREQ_B1_pools, AF_2 = locus_RefFREQ_B4_pools)

intralocus_RefFREQ_signfst <- intralocus_RefFREQ[intralocus_RefFREQ$CHROM %in% unique(high.multilocus.fst.ordered$CHRM), ]
#table(intralocus_RefFREQ_signfst$CHROM)

# Split the table in scaffolds and produce a list - to parallel the code
listOf_signfstScaffolds <- split(intralocus_RefFREQ_signfst, intralocus_RefFREQ_signfst$CHROM)

# Apply the avgFstSlidingWindow function on each scaffold
intralocus_RefFREQ_signfst_window_cor <- do.call(rbind, lapply(listOf_signfstScaffolds, 
                                                               alleleRefFreqCorSlidingWindow, 
                                                               step=100, windowSize=10000))

### Multi ALLELE FREQUENCY CORRELATIONS ALONG CHROM plots for all significant scaffolds
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/significant_intralocus_fst_RefFREQcor_scaffolds.pdf", # File name
    width = 8.5, height = 11,                                                                                       # Width and height in inches
    bg = "white",                                                                                                   # Background color
    colormodel = "cmyk",                                                                                            # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(4, 4)) # 16 chromosomes
for (s in seq_along(unique(intralocus_RefFREQ_signfst_window_cor$CHROM)))
{
  chr_name <- unique(intralocus_RefFREQ_signfst_window_cor$CHROM)[s]
  n <- intralocus_RefFREQ_signfst_window_cor[intralocus_RefFREQ_signfst_window_cor$CHROM %in% chr_name,]
  
  plot(n$windowCor ~ n$Step, type="n", ylim=c(-1,1), 
       xlab="Window step",ylab="Correlation", main=chr_name,
       col="black",pch=20, cex.lab=1.4, cex.axis=1.2)
  lines(lowess(n$windowCor ~ n$Step, f=0.10), col="orange", lwd=2, lty=1)
}
dev.off()

### BIOTYPE-AVERAGES allele frequencies from the ML reference allele counts pooled by biotype
# BIOTYPE 1
ref_count_b1 <- apply(dt.1.imputedRefMLCount$ref_count[,-c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) sum(x, na.rm = T))
read_count_b1 <- apply(dt.1.imputedRefMLCount$hap_count[,-c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) sum(x, na.rm = T))

ref_allele_freq_b1 <- ref_count_b1/read_count_b1

# BIOTYPE 4
ref_count_b4 <- apply(dt.1.imputedRefMLCount$ref_count[,c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) sum(x, na.rm = T))
read_count_b4 <- apply(dt.1.imputedRefMLCount$hap_count[,c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) sum(x, na.rm = T))

ref_allele_freq_b4 <- ref_count_b4/read_count_b4

biotypes_ref_allele_freq  <- data.frame(CHROM = dt.1@snp.info[,1], POS = dt.1@snp.info[,2],
                                         biotype_1 = ref_allele_freq_b1, biotype_4 = ref_allele_freq_b4)

#max(abs(biotypes_ref_allele_freq$biotype_1 - biotypes_ref_allele_freq$biotype_4))
#biotypes_ref_allele_freq[which(abs(biotypes_ref_allele_freq$biotype_1 - biotypes_ref_allele_freq$biotype_4) > 0.5), ]

### EXPORT ADDITIONAL SNP INFORMATION
# SNP COVERAGE
locus_min_cov_pools <- apply(dt.1.read_coverage_data, 1, function(x) min(x, na.rm = T))
locus_max_cov_pools <- apply(dt.1.read_coverage_data, 1, function(x) max(x, na.rm = T))

# SNP REFERENCE FREQUENCY
locus_min_refMLFreq_pools <- apply(dt.1.imputedRefMLFreq, 1, function(x) min(x, na.rm = T))
locus_max_refMLFreq_pools <- apply(dt.1.imputedRefMLFreq, 1, function(x) max(x, na.rm = T))

snp_info_aggregated_pools <- data.frame(dt.1@snp.info, 
                                        MEAN_COV=locus_ADP_pools, MIN_COV=locus_min_cov_pools, MAX_COV=locus_max_cov_pools, 
                                        MEAN_REFFREQ=dt.1.overall.REFallele.freq, MIN_REFFREQ=locus_min_refMLFreq_pools, MAX_REFFREQ=locus_max_refMLFreq_pools,
                                        MISSING=dt.1.missing.snps)

# Save to a file all the information
write.table(snp_info_aggregated_pools,
            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/snpinfo_aggregated_data_pools.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

# Save to a file only the combined chromosome and SNP position
write.table(paste0(snp_info_aggregated_pools$Chromosome, "_", snp_info_aggregated_pools$Position),
            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/snpinfo_aggregated_names_pools.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)


###
###
### ---- COMPUTE GLOBAL AND PAIRWISE FST WITH ALL POOLS ----
###
###

### ALL 21 POOLS
### Dataset 2: 24-Jun-2021; Min_cov=4; Max_cov=0.99; MAF=0.05; MAC=5; 21 pools; miss_fraction=0.50
### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz

### COMPUTE THE GLOBAL FST
dt.1.fst <- computeFST(dt.1, method = "Anova", nsnp.per.bjack.block = 10, verbose=FALSE)

# genome-wide estimate
dt.1.fst$FST 
# 0.004362626

# mean Block-Jackknife estimation
dt.1.fst$mean.fst
# 0.004459117

# 95% CI
dt.1.fst$mean.fst+c(-1.96,1.96)*dt.1.fst$se.fst
# 0.002462227 0.006456007

### COMPUTE PAIRWISE FST
#dt.1.pfst <- compute.pairwiseFST(dt.1,verbose=FALSE)
#dt.1.pfst <- compute.pairwiseFST(dt.1,verbose=FALSE, output.snp.values = T) 
# to have SNP-specific FST for all pairwise comparisons
# it can be used to average across biotypes and get the fst tracks

# with blockjackniff
dt.1.pfst <- compute.pairwiseFST(dt.1, nsnp.per.bjack.block = 100, verbose=TRUE)

dt.1.pfst.bjack.table <- data.frame(dt.1.pfst@values, 
                                    ci95lw = dt.1.pfst@values$`Fst bjack mean`+c(-1.96)*dt.1.pfst@values$`Fst bjack s.e.`,
                                    ci95up = dt.1.pfst@values$`Fst bjack mean`+c(1.96)*dt.1.pfst@values$`Fst bjack s.e.`)

dt.1.pfst.bjack.table[which.min(dt.1.pfst.bjack.table$Fst.Estimate), ]
dt.1.pfst.bjack.table[which.max(dt.1.pfst.bjack.table$Fst.Estimate), ]

# HEATMAP OF PAIRWISE FST
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pFST_pools_R-heatmap.pdf",   # File name
    width = 8.5, height = 11,                                                              # Width and height in inches
    bg = "white",                                                                          # Background color
    colormodel = "cmyk",                                                                   # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
heatmap(dt.1.pfst)
dev.off()

# PAIRWISE FST IC 95% PLOT
#plot_fstats(dt.1.pfst, cex=0.4)

### PYTHON SEABORN HEATMAP
# Convert Matrix of pFST to Numpy Array
pfst_matrix <- dt.1.pfst@PairwiseFSTmatrix
pfst_matrix[is.na(pfst_matrix)] <- 0
pfst_matrix <- as.matrix(pfst_matrix)
pfst_numpay_array <- r_to_py(pfst_matrix)

# Convert Numpy Array to Pandas DataFrame
pfst_pandas_df <- pd$DataFrame(data=pfst_numpay_array)
colnames(pfst_pandas_df) <- poolnames
rownames(pfst_pandas_df) <- poolnames

# SEABORN HEATMAP OF FSTS
sns$clustermap(pfst_pandas_df, 
              center=0, cmap="vlag",
              row_colors = biotype.col, col_colors = biotype.col,
              dendrogram_ratio=c(.1, .2),
              #cbar_pos=c(.01, .32, .02, .2),
              linewidths=.75, figsize=c(10, 11),
              row_cluster=TRUE)
plt$savefig("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pFST_pools_Seaborn-heatmap.pdf")

  
### ONE SAMPLE OF EACH BIOTYPE / LOCALITY
### SUBSET THE POOLFSTAT OBJECT TO KEEP ONLY ONE BIOTYPE/LOCALITY
### USE THE POOL WITH LOW % OF MISSING DATA IF MORE THAN ONE BIOTYPE POOL / LOCALITY
### Dataset 2: 24-Jun-2021; Min_cov=4; Max_cov=0.99; MAF=0.05; MAC=5; 21 pools; miss_fraction=0.50
### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz
### REMOVE SNPS WITH > 25% OF MISSING DATA AFTER POOL REMOVAL

# KEEP Biotype pools with less missing data in each locality
biotype.pools2keep <- c(1, 2, 3, 4, 6, 7,
                        10,11,14,16,18,20)

biotype.pools2remove <- setdiff(x=seq(1,21,1), y=biotype.pools2keep)

# remove the other pools from the imputed allele frequency matrix
dt.1.imputedRefMLFreq.biotypes  <- dt.1.imputedRefMLFreq[,-c(biotype.pools2remove)]

### REMOVE SNPS WITH > 25% AFTER THE REMOVAL OF THE SAMPLES

# Calculate the % of missing data 
dt.1.missing.snps.biotypes <- (apply(dt.1.imputedRefMLFreq.biotypes, 1, function(x) sum(is.na(x)))/(length(biotype.pools2keep)))*100
hist(dt.1.missing.snps.biotypes, breaks = 10)

length(dt.1.missing.snps.biotypes)
min(dt.1.missing.snps.biotypes) #0
max(dt.1.missing.snps.biotypes) #75

# Filter SNPs
keep.snps.biotypes <- which(dt.1.missing.snps.biotypes <= 25)
length(keep.snps.biotypes) # 201573

# How many were removed?
length(dt.1.missing.snps.biotypes) - length(keep.snps.biotypes) # 61293

# SANITY CHECK MISSING % AFTER SNP REMOVAL
dt.1.imputedRefMLFreq.biotypes.filtered025 <- dt.1.imputedRefMLFreq.biotypes[keep.snps.biotypes, ] 
dt.1.missing.snps.biotypes.filtered025 <- (apply(dt.1.imputedRefMLFreq.biotypes.filtered025, 1, function(x) sum(is.na(x)))/length(biotype.pools2keep))*100

length(dt.1.missing.snps.biotypes.filtered025) # 201573
min(dt.1.missing.snps.biotypes.filtered025)
max(dt.1.missing.snps.biotypes.filtered025)

### SUBSET THE DATA
dt.1.biotypes  <- pooldata.subset(
                                  pooldata=dt.1,
                                  pool.index = biotype.pools2keep,
                                  snp.index = keep.snps.biotypes,
                                  min.cov.per.pool = -1,
                                  max.cov.per.pool = 1e+06,
                                  min.maf = -1,
                                  cov.qthres.per.pool = c(0, 1),
                                  return.snp.idx = TRUE,
                                  verbose = TRUE
)
dt.1.biotypes

### EXPORT THE IMPUTED ML ALLELE FREQUENCY MATRIX TO A FILE
### WITH THE SNP INFORMATION AND POOL NAMES (FOR THE PCA ANALYSIS)
dt.1.imputedRefMLFreq.biotypes.filtered025.data <- dt.1.imputedRefMLFreq.biotypes.filtered025
colnames(dt.1.imputedRefMLFreq.biotypes.filtered025.data) <- dt.1.biotypes@poolnames

dt.1.imputedRefMLFreq.biotypes.filtered025.dataFrame <- data.frame(dt.1.biotypes@snp.info, dt.1.imputedRefMLFreq.biotypes.filtered025.data)

# Save to a file
#write.table(dt.1.imputedRefMLFreq.biotypes.filtered025.dataFrame,
#            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_biotypes.txt", sep = "\t",
#            quote = F, row.names = F, col.names = T)

### COMPUTE THE GLOBAL FST
dt.1.biotypes.fst <- computeFST(dt.1.biotypes, method = "Anova", nsnp.per.bjack.block = 10, verbose=FALSE)

# genome-wide estimate
dt.1.biotypes.fst$FST 
# 0.001619748

# mean Block-Jackknife estimation
dt.1.biotypes.fst$mean.fst
# 0.001700997

# 95% CI
dt.1.biotypes.fst$mean.fst + c(-1.96,1.96) * dt.1.biotypes.fst$se.fst
# 0.0002086295 0.0031933653

### COMPUTE PAIRWISE FST
dt.1.biotypes.pfst <- compute.pairwiseFST(dt.1.biotypes,verbose=FALSE)

# with blockjackniff
dt.1.biotypes.pfst <- compute.pairwiseFST(dt.1.biotypes, nsnp.per.bjack.block = 100, verbose=FALSE)

# HEATMAP OF PAIRWISE FST
pdf("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pFST_biotypes_R-heatmap.pdf",   # File name
    width = 8.5, height = 11,                                                                 # Width and height in inches
    bg = "white",                                                                             # Background color
    colormodel = "cmyk",                                                                      # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
heatmap(dt.1.biotypes.pfst)
dev.off()

# PAIRWISE FST IC 95% PLOT
#plot_fstats(dt.1.biotypes.pfst, cex=0.4)

### PYTHON SEABORN HEATMAP
# Convert Matrix of pFST to Numpy Array
pfst_biotypes_matrix <- dt.1.biotypes.pfst@PairwiseFSTmatrix
pfst_biotypes_matrix[is.na(pfst_biotypes_matrix)] <- 0
pfst_biotypes_matrix <- as.matrix(pfst_biotypes_matrix)
pfst_biotypes_numpay_array <- r_to_py(pfst_biotypes_matrix)

# Convert Numpy Array to Pandas DataFrame
pfst_biotypes_pandas_df <- pd$DataFrame(data=pfst_biotypes_numpay_array)
colnames(pfst_biotypes_pandas_df) <- poolnames[-c(biotype.pools2remove)]
rownames(pfst_biotypes_pandas_df) <- poolnames[-c(biotype.pools2remove)]

# SEABORN HEATMAP OF FSTS
sns$clustermap(pfst_biotypes_pandas_df, 
               center=0, cmap="vlag",
               row_colors = biotype.col[-c(biotype.pools2remove)], col_colors = biotype.col[-c(biotype.pools2remove)],
               dendrogram_ratio=c(.1, .2),
               #cbar_pos=c(.01, .32, .02, .2),
               linewidths=.75, figsize=c(10, 11),
               row_cluster=TRUE)
plt$savefig("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pFST_biotypes_Seaborn-heatmap.pdf")

## Savepoint_3
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
############## END THIS PART

###
###
### ---- PARTIAL MANTEL TEST AND CAUSAL MODELLING ----
###
###

### ALL 21 SAMPLES
### Dataset 2: 24-Jun-2021; Min_cov=4; Max_cov=0.99; MAF=0.05; MAC=5; 21 samples; miss_fraction=0.50
### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz

## CAUSAL MODEL: 
## Partial Mentel Test: pFST vs geographic distance | cost matrix (same vs different biotypes)
## Partial Mentel Test: pFST vs geographic distance | cost matrix (biotype 1 vs biotype 4 vs different biotypes)
## FOR ALL SAMPLES

# Matrix of linearized pFST
pfst_lin_matrix <- dt.1.pfst@PairwiseFSTmatrix
pfst_lin_matrix <- pfst_lin_matrix / (1- pfst_lin_matrix)
pfst_lin_matrix[is.na(pfst_lin_matrix)] <- 0
pfst_lin_matrix <- as.matrix(pfst_lin_matrix)

# Load the NxN Matrix of the geographic distance between Pools in Km
geo_matrix <- read.table("DEST-AglyPoolseq/populationInfo/filedPools_aggregated_geocoord_pools.matrix.txt", header=TRUE)
geo_matrix <- as.matrix(geo_matrix)

# Load the NxN Matrix of the cost distance
# here is a design matrix were pools collected in the same host have the same identification number
cost_matrix_1 <-read.table("DEST-AglyPoolseq/populationInfo/fieldPools_aggregated_cost_pools.matrix1.txt", header=TRUE)
cost_matrix_1 <-as.matrix(cost_matrix_1)

# Partial mantel test with pFST and cost matrix 1
wassermann(eucl.mat = geo_matrix, cost.mats = list(cost=cost_matrix_1), gen.mat = pfst_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.1834 0.0748
#1 Gen ~cost | Euclidean 0.0759 0.1651

# Partial mantel test with linearized pFST and cost matrix 1
wassermann(eucl.mat = geo_matrix, cost.mats = list(cost=cost_matrix_1), gen.mat = pfst_lin_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.1817 0.081
#1 Gen ~cost | Euclidean 0.0755 0.1724

#  Load the NxN Matrix of the cost distance 2- 1-biotype1; 2-biotype4, 3-different biotypes
cost_matrix_2 <-read.table("DEST-AglyPoolseq/populationInfo/fieldPools_aggregated_cost_pools.matrix2.txt", header=TRUE)
cost_matrix_2 <-as.matrix(cost_matrix_2)

# Partial mantel test with pFST and cost matrix 2
wassermann(eucl.mat = geo_matrix, cost.mats = list(cost=cost_matrix_2), gen.mat = pfst_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.1813 0.0775
#1 Gen ~cost | Euclidean 0.0174 0.302

# Partial mantel test with linearized pFST and cost matrix 2
wassermann(eucl.mat = geo_matrix, cost.mats = list(cost=cost_matrix_2), gen.mat = pfst_lin_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.1797 0.0834
#1 Gen ~cost | Euclidean 0.0182 0.2909

## Mantel test genetic dist vs geo dist
pfst_matrix.dist <- as.dist(pfst_matrix)
pfst_lin_matrix.dist <- as.dist(pfst_lin_matrix) # linearized fst

geo_matrix.dist <- as.dist(geo_matrix)
geo_matrix.log.dist <-log(geo_matrix.dist)

cost_matrix_1.dist <- as.dist(cost_matrix_1)
cost_matrix_2.dist <- as.dist(cost_matrix_2)

mantel_test_table_pools <- data.frame(lin_pFST = as.vector(pfst_lin_matrix.dist),
                                      geo_dist = as.vector(geo_matrix.dist),
                                      cost_dist_1 = as.vector(cost_matrix_1.dist),
                                      cost_dist_2 = as.vector(cost_matrix_2.dist))

### PYTHON SEABORN LMPLOT
sns$set_theme()
sns$set_palette("colorblind")

# MATRIX RELATIONSHIP COST MATRIX 1
plt$figure(figsize=c(5,6))
sns$lmplot(x="geo_dist", y="lin_pFST", hue="cost_dist_1", 
           data=r_to_py(mantel_test_table_pools), 
           legend=FALSE, height=5)
plt$xlim(c(-0.5,1725))
plt$xlabel("Geographical distances (Km)")
plt$ylabel("Pairwise FST")
plt$legend(title="Cost matrix", fontsize=9, 
           labels=c("Same Biotype pools", "Avirulent vs Virulent"),
           loc='upper right')
plt$tight_layout()
plt$savefig("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pMantel_test_cost1_pools_Seaborn-lmplot.pdf")

# MATRIX RELATIONSHIP COST MATRIX 2
plt$figure(figsize=c(5,6))
sns$lmplot(x="geo_dist", y="lin_pFST", hue="cost_dist_2", 
           data=r_to_py(mantel_test_table_pools), 
           legend=FALSE, height=5)
plt$xlim(c(-0.5,1725))
plt$xlabel("Geographical distances (Km)")
plt$ylabel("Pairwise FST")
plt$legend(title="Cost matrix", fontsize=9, 
           labels=c("Avirulent", "Virulent", "Avirulent vs Virulent"),
           loc='upper right')
plt$tight_layout()
plt$savefig("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pMantel_test_cost2_pools_Seaborn-lmplot.pdf")

mantel.randtest(pfst_lin_matrix.dist, geo_matrix.dist, nrepet = 1000)
mantel.randtest(pfst_matrix.dist, geo_matrix.dist, nrepet = 1000)
mantel.randtest(pfst_lin_matrix.dist, geo_matrix.log.dist, nrepet = 1000)
mantel.randtest(pfst_matrix.dist, geo_matrix.log.dist, nrepet = 1000)

## LINEAR REGRESSION pFST ~ GEO DISTANCE | COST MATRIX

# COST MATRIX 1
lm1_cost_1 <- lm(mantel_test_table_pools$lin_pFST[which(mantel_test_table_pools$cost_dist_1 == 1)] ~ mantel_test_table_pools$geo_dist[which(mantel_test_table_pools$cost_dist_1 == 1)])
summary(lm1_cost_1)

lm2_cost_1 <- lm(mantel_test_table_pools$lin_pFST[which(mantel_test_table_pools$cost_dist_1 == 2)] ~ mantel_test_table_pools$geo_dist[which(mantel_test_table_pools$cost_dist_1 == 2)])
summary(lm2_cost_1)

# COST MATRIX 2
lm1_cost_2 <- lm(mantel_test_table_pools$lin_pFST[which(mantel_test_table_pools$cost_dist_2 == 1)] ~ mantel_test_table_pools$geo_dist[which(mantel_test_table_pools$cost_dist_2 == 1)])
summary(lm1_cost_2)

lm2_cost_2 <- lm(mantel_test_table_pools$lin_pFST[which(mantel_test_table_pools$cost_dist_2 == 2)] ~ mantel_test_table_pools$geo_dist[which(mantel_test_table_pools$cost_dist_2 == 2)])
summary(lm2_cost_2)

lm3_cost_2 <- lm(mantel_test_table_pools$lin_pFST[which(mantel_test_table_pools$cost_dist_2 == 3)] ~ mantel_test_table_pools$geo_dist[which(mantel_test_table_pools$cost_dist_2 == 3)])
summary(lm3_cost_2)


### ONE SAMPLE OF EACH BIOTYPE / LOCALITY
### SUBSET THE POOLFSTAT OBJECT TO KEEP ONLY ONE BIOTYPE/LOCALITY (ALREADY DONE ABOVE)
### USE THE SAMPLE WITH LOW % OF MISSING DATA IF MORE THAN ONE BIOTYPE SAMPLE / LOCALITY
### Dataset 2: 24-Jun-2021; Min_cov=4; Max_cov=0.99; MAF=0.05; MAC=5; 21 pools; miss_fraction=0.50
### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz
### REMOVE SNPS WITH > 25% OF MISSING DATA AFTER POOL REMOVAL

## CAUSAL MODEL: 
## Partial Mentel Test: pFST vs geographic distance | cost matrix (same vs different biotypes)
## Partial Mentel Test: pFST vs geographic distance | cost matrix (biotype 1 vs biotype 4 vs different biotypes)
## FOR ALL POOLS

# Matrix of linearized pFST
pfst_biotypes_lin_matrix <- dt.1.biotypes.pfst@PairwiseFSTmatrix
pfst_biotypes_lin_matrix <- pfst_biotypes_lin_matrix / (1- pfst_biotypes_lin_matrix)
pfst_biotypes_lin_matrix[is.na(pfst_biotypes_lin_matrix)] <- 0
pfst_biotypes_lin_matrix <- as.matrix(pfst_biotypes_lin_matrix)

# Load the NxN Matrix of the geographic distance between Pools in Km
geo_biotypes_matrix <- read.table("DEST-AglyPoolseq/populationInfo/filedPools_aggregated_geocoord_biotypes.matrix.txt", header=TRUE)
geo_biotypes_matrix <- as.matrix(geo_biotypes_matrix)

# Load the NxN Matrix of the cost distance
# here is a design matrix were pools collected in the same host have the same identification number
cost_biotypes_matrix_1 <-read.table("DEST-AglyPoolseq/populationInfo/fieldPools_aggregated_cost_biotypes.matrix1.txt", header=TRUE)
cost_biotypes_matrix_1 <-as.matrix(cost_biotypes_matrix_1)

# Change to row and column names of the pfst biotypes matrix
rownames(pfst_biotypes_matrix) <- rownames(geo_biotypes_matrix) # From pFST HEATMAP ABOVE
colnames(pfst_biotypes_matrix) <- colnames(geo_biotypes_matrix)

rownames(pfst_biotypes_lin_matrix) <- rownames(geo_biotypes_matrix)
colnames(pfst_biotypes_lin_matrix) <- colnames(geo_biotypes_matrix)

# Partial mantel test with pFST and cost matrix 1
wassermann(eucl.mat = geo_biotypes_matrix, cost.mats = list(cost=cost_biotypes_matrix_1), gen.mat = pfst_biotypes_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.3576 0.0014
#1 Gen ~cost | Euclidean 0.0759 0.1651

# Partial mantel test with linearized pFST and cost matrix 1
wassermann(eucl.mat = geo_biotypes_matrix, cost.mats = list(cost=cost_biotypes_matrix_1), gen.mat = pfst_biotypes_lin_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.3598  0.001
#1 Gen ~cost | Euclidean 0.0587 0.1353

#  Load the NxN Matrix of the cost distance 2- 1-biotype1; 2-biotype4, 3-different biotypes
cost_biotypes_matrix_2 <-read.table("DEST-AglyPoolseq/populationInfo/fieldPools_aggregated_cost_biotypes.matrix2.txt", header=TRUE)
cost_biotypes_matrix_2 <-as.matrix(cost_biotypes_matrix_2)

# Partial mantel test with pFST and cost matrix 2
wassermann(eucl.mat = geo_biotypes_matrix, cost.mats = list(cost=cost_biotypes_matrix_2), gen.mat = pfst_biotypes_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.3536 0.0013
#1 Gen ~cost | Euclidean 0.0233  0.447

# Partial mantel test with linearized pFST and cost matrix 2
wassermann(eucl.mat = geo_biotypes_matrix, cost.mats = list(cost=cost_biotypes_matrix_2), gen.mat = pfst_biotypes_lin_matrix, nperm=10000, plot=F)$mantel.tab
#                   model      r      p
#2 Gen ~Euclidean | cost 0.3556 0.0011
#1 Gen ~cost | Euclidean 0.0238 0.4509

## Mantel test genetic dist vs geo dist
pfst_biotypes_matrix.dist <- as.dist(pfst_biotypes_matrix)
pfst_biotypes_lin_matrix.dist <- as.dist(pfst_biotypes_lin_matrix) # linearized fst

geo_biotypes_matrix.dist <- as.dist(geo_biotypes_matrix)
geo_biotypes_matrix.log.dist <-log(geo_biotypes_matrix.dist)

cost_biotypes_matrix_1.dist <- as.dist(cost_biotypes_matrix_1)
cost_biotypes_matrix_2.dist <- as.dist(cost_biotypes_matrix_2)

mantel_test_table_biotypes <- data.frame(lin_pFST = as.vector(pfst_biotypes_lin_matrix.dist),
                                         geo_dist = as.vector(geo_biotypes_matrix.dist),
                                         cost_dist_1 = as.vector(cost_biotypes_matrix_1.dist),
                                         cost_dist_2 = as.vector(cost_biotypes_matrix_2.dist))

### PYTHON SEABORN LMPLOT
sns$set_theme()
sns$set_palette("colorblind")

# MATRIX RELATIONSHIP COST MATRIX 1
plt$figure(figsize=c(5,6))
sns$lmplot(x="geo_dist", y="lin_pFST", hue="cost_dist_1", 
           data=r_to_py(mantel_test_table_biotypes), 
           legend=FALSE, height=5)
plt$xlim(c(-0.5,1725))
plt$xlabel("Geographical distances (Km)")
plt$ylabel("Pairwise FST")
plt$legend(title="Cost matrix", fontsize=9, 
           labels=c("Same Biotype pools", "Avirulent vs Virulent"),
           loc='upper right')
plt$tight_layout()
plt$savefig("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pMantel_test_cost1_biotypes_Seaborn-lmplot.pdf")

# MATRIX RELATIONSHIP COST MATRIX 2
plt$figure(figsize=c(5,6))
sns$lmplot(x="geo_dist", y="lin_pFST", hue="cost_dist_2", 
           data=r_to_py(mantel_test_table_biotypes), 
           legend=FALSE, height=5)
plt$xlim(c(-0.5,1725))
plt$xlabel("Geographical distances (Km)")
plt$ylabel("Pairwise FST")
plt$legend(title="Cost matrix", fontsize=9, 
           labels=c("Avirulent", "Virulent", "Avirulent vs Virulent"),
           loc='upper right')
plt$tight_layout()
plt$savefig("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pMantel_test_cost2_biotypes_Seaborn-lmplot.pdf")

mantel.randtest(pfst_biotypes_lin_matrix.dist, geo_biotypes_matrix.dist, nrepet = 1000)
mantel.randtest(pfst_biotypes_matrix.dist, geo_biotypes_matrix.dist, nrepet = 1000)
mantel.randtest(pfst_biotypes_lin_matrix.dist, geo_biotypes_matrix.log.dist, nrepet = 1000)
mantel.randtest(pfst_biotypes_matrix.dist, geo_biotypes_matrix.log.dist, nrepet = 1000)

## LINEAR REGRESSION pFST ~ GEO DISTANCE | COST MATRIX

# COST MATRIX 1
lm1_bio_cost_1 <- lm(mantel_test_table_biotypes$lin_pFST[which(mantel_test_table_biotypes$cost_dist_1 == 1)] ~ mantel_test_table_biotypes$geo_dist[which(mantel_test_table_biotypes$cost_dist_1 == 1)])
summary(lm1_bio_cost_1)

lm2_bio_cost_1 <- lm(mantel_test_table_biotypes$lin_pFST[which(mantel_test_table_biotypes$cost_dist_1 == 2)] ~ mantel_test_table_biotypes$geo_dist[which(mantel_test_table_biotypes$cost_dist_1 == 2)])
summary(lm2_bio_cost_1)

# COST MATRIX 2
lm1_bio_cost_2 <- lm(mantel_test_table_biotypes$lin_pFST[which(mantel_test_table_biotypes$cost_dist_2 == 1)] ~ mantel_test_table_biotypes$geo_dist[which(mantel_test_table_biotypes$cost_dist_2 == 1)])
summary(lm1_bio_cost_2)

lm2_bio_cost_2 <- lm(mantel_test_table_biotypes$lin_pFST[which(mantel_test_table_biotypes$cost_dist_2 == 2)] ~ mantel_test_table_biotypes$geo_dist[which(mantel_test_table_biotypes$cost_dist_2 == 2)])
summary(lm2_bio_cost_2)

lm3_bio_cost_2 <- lm(mantel_test_table_biotypes$lin_pFST[which(mantel_test_table_biotypes$cost_dist_2 == 3)] ~ mantel_test_table_biotypes$geo_dist[which(mantel_test_table_biotypes$cost_dist_2 == 3)])
summary(lm3_bio_cost_2)

## Savepoint_4
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace22Jul21.RData")
############## END THIS PART