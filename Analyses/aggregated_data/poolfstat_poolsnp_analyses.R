##########################################################
#                  POOLFSTAT package                     # 
#         AF, HE, Theta, pF_ST, gF_ST and genome scan    #
##########################################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES

setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()
## Remove last features
rm(list=ls())
ls()

library(poolfstat)
packageVersion("poolfstat")   
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

## POOL'S INFORMATION
poolsizes <- rep(10,21)
poolnames <- c("MN-B1.1", "MN-B4.1",
               "ND-B1.1", "ND-B4.1", 
               "NW-B1.1", "NW-B1.2", "NW-B4.1", "NW-B4.2", "NW-B4.3", 
               "PA-B1.1", "PA-B4.1", "PA-B4.2", 
               "WI-B1.1", "WI-B1.2", "WI-B4.1", "WI-B4.2", "WI-B4.3", 
               "WO-B1.1", "WO-B4.1", "WO-B4.2", "WO-B4.3")


biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53")

####
### Data set 1: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=5; 21 pools; miss_fraction=0.50
####
#
### INVESTIGATE HOW MAC DETERMINE THE MAF
### THIS DATASET WAS ONLY USED FOR IT
#testcase.dt.1 <- vcf2pooldata(
#                              vcf.file = "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz",
#                              poolsizes = poolsizes,
#                              poolnames = poolnames,
#                              min.cov.per.pool = -1,
#                              min.rc = 1,
#                              max.cov.per.pool = 1e+06,
#                              min.maf = 0.001,
#                              remove.indels = FALSE,
#                              nlines.per.readblock = 1e+06
#)
#
#testcase.dt.1
## * * * PoolData Object * * * 
## * Number of SNPs   = 344524 
## * Number of Pools  =  21 
## Number of scaffolds:
#
### CHECK THE FREQUENCY OF EACH MAC CLASS FROM THE ORIGINAL READ COUNT DATA
#testcase.ALTallele.count <- (testcase.dt.1@readcoverage - testcase.dt.1@refallele.readcount)
#
#testcase.overall.ALTallele.count <- apply(testcase.ALTallele.count, 1, sum)
#testcase.overall.REFallele.count <- apply(testcase.dt.1@refallele.readcount, 1, sum)
#
### MAC
#overall.MAC <- pmin(testcase.overall.REFallele.count, testcase.overall.ALTallele.count)
#
### MAF
#testcase.imputedMLCount <- imputedRefMLCount(testcase.dt.1)
#testcase.imputedRefMLFreq <- testcase.imputedMLCount[[1]]/testcase.imputedMLCount[[3]]
#
#testcase.overall.ALTallele.freq <- apply(testcase.imputedRefMLFreq, 1,  function(x) 1-mean(x, na.rm = T))
#testcase.overall.REFallele.freq <- apply(testcase.imputedRefMLFreq, 1,  function(x) mean(x, na.rm = T))
#
#overall.MAF <- pmin(testcase.overall.REFallele.freq, testcase.overall.ALTallele.freq)
#
### HOW MAF CHANGE WITH MAC
#macs <- seq(5, 100, 5)
#mafs.macs.table <- NULL
#for (i in seq_along(macs))
#{
#  # REF/ALT Alleles
#  a <- testcase.overall.REFallele.freq[overall.MAC >= macs[i]] 
#  b <- testcase.overall.ALTallele.freq[overall.MAC >= macs[i]] 
#  ref.min <- min(a); ref.max <- max(a)
#  alt.min <- min(b); alt.max <- max(b)
#  
#  # MAF
#  p <- pmin(a,b)
#  maf <- min(p)
#  
#  #OUTPUT
#  res <- data.frame(MAC=macs[i], MAF=maf, REF_MIN=ref.min, REF_MAX=ref.max, ALT_MIN=alt.min, ALT_MAX=alt.max)
#  mafs.macs.table <- rbind(mafs.macs.table,res)
#}
#
### MAF AS A FUNCTION OF MAC
#plot(mafs.macs.table$MAC, mafs.macs.table$MAF,  type = "l", lty = 1, col="black", xlab="MAC", ylab="MAF", ylim=c(0,0.15))
#lines(mafs.macs.table$MAC, mafs.macs.table$REF_MIN,  type = "l", lty = 1, col="blue")
##lines(mafs.macs.table$MAC, mafs.macs.table$REF_MAX,  type = "l", lty = 2, col="blue")
#lines(mafs.macs.table$MAC, mafs.macs.table$ALT_MIN,  type = "l", lty = 1, col="red")
##lines(mafs.macs.table$MAC, mafs.macs.table$ALT_MAX,  type = "l", lty = 2, col="red")
#
## MAC         MAF    REF_MIN
##   5 0.004761905 0.01904762 
##  10 0.014285714 0.03529412 
##  15 0.019047619 0.04210526 
##  20 0.022222222 0.05555556 
##  25 0.023809524 0.08500000 
##  30 0.023809524 0.09000000 
##  35 0.023809524 0.09000000 
##  40 0.023809524 0.09000000 
##  45 0.023809524 0.09005848 
##  50 0.023809524 0.09473684 
## 100 0.061904762 0.10333333 
#

###
## Data set 2: 24-Jun-2021; Min_cov=4; Max_cov=0.99; MAF=0.05; MAC=5; 21 pools; miss_fraction=0.50
###

## POPULATION GENOMICS ANALYSES
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
# * Number of SNPs   = 324356
# * Number of Pools  =  21 
# Number of scaffolds:

length(unique(sort(dt.1@snp.info$Chromosome))) #705

## CHECK RAW MAC
dt.1.ALTallele.count <- (dt.1@readcoverage - dt.1@refallele.readcount)

dt.1.overall.ALTallele.count <- apply(dt.1.ALTallele.count, 1, sum)
dt.1.overall.REFallele.count <- apply(dt.1@refallele.readcount, 1, sum)

dt.1.MAC <- pmin(dt.1.overall.REFallele.count, dt.1.overall.ALTallele.count)
min(dt.1.MAC); max(dt.1.MAC)

## CHECK RAW MAF
# COMPUTE MAXIMUM LIKELIHOOD IMPUTED SAMPLE ALLELE COUNT
dt.1.imputedRefMLCount <- imputedRefMLCount(dt.1)
dt.1.imputedRefMLFreq <- dt.1.imputedRefMLCount[[1]]/dt.1.imputedRefMLCount[[3]]

dt.1.overall.ALTallele.freq <- apply(dt.1.imputedRefMLFreq, 1,  function(x) 1-mean(x, na.rm = T))
dt.1.overall.REFallele.freq <- apply(dt.1.imputedRefMLFreq, 1,  function(x) mean(x, na.rm = T))

# Sanity check if the MAF correspond to the MAF used in filtering
dt.1.MAF <- pmin(dt.1.overall.REFallele.freq, dt.1.overall.ALTallele.freq)
min(dt.1.MAF); max(dt.1.MAF)
hist(dt.1.MAF)

# HOW MANY SNPS IS FIXED IN 95% OF POOLS
hist(apply(dt.1.imputedRefMLFreq == 1, 1, function(x) sum(x, na.rm = T)))
table(apply(dt.1.imputedRefMLFreq == 1, 1, function(x) sum(x, na.rm = T)))
quasi.fixed <- which(apply(dt.1.imputedRefMLFreq == 1, 1, function(x) sum(x, na.rm = T)) == 20) 

# Check the MAF for almost fixed markers
#dt.1.fixed.ALTallele.freq <- apply(dt.1.imputedRefMLFreq[quasi.fixed ,], 1,  function(x) 1-mean(x, na.rm = T))
#dt.1.fixed.REFallele.freq <- apply(dt.1.imputedRefMLFreq[quasi.fixed ,], 1,  function(x) mean(x, na.rm = T))
#dt.1.MAF.fixed <- pmin(dt.1.fixed.REFallele.freq, dt.1.fixed.ALTallele.freq)
#min(dt.1.MAF.fixed); max(dt.1.MAF.fixed)

## CHECK THE % OF MISSING DATA
# POOLS
dt.1.missing.pools <- (apply(dt.1.imputedRefMLFreq, 2, function(x) sum(is.na(x)))/dt.1@nsnp)*100
names(dt.1.missing.pools) <- poolnames

barplot(dt.1.missing.pools, ylim = c(0,100))
abline(h=30, col="red")

# SNPs
dt.1.missing.snps <- (apply(dt.1.imputedRefMLFreq, 1, function(x) sum(is.na(x)))/dt.1@npools)*100
hist(dt.1.missing.snps, breaks = 10)

# REMOVE 2 SAMPLES IDENTIFIED WITH > 30% MISSING DATA (DONT REMOVE PA SAMPLES)
#reduced.dt.1.imputedRefMLFreq <- dt.1.imputedRefMLFreq[,-c(5,21)]
#
#dt.1.missing.snps.reduced <- (apply(reduced.dt.1.imputedRefMLFreq, 1, function(x) sum(is.na(x)))/(dt.1@npools-2))*100
#hist(dt.1.missing.snps.reduced, breaks = 10)
#
### POPULATION GENOMICS ANALYSES
## COMPUTE THE GLOBAL FST
dt.1.fst <- computeFST(dt.1, method = "Anova", nsnp.per.bjack.block = 10, verbose=FALSE)

# genome-wide estimate
dt.1.fst$FST 
# 0.004362626

# mean Block-Jackknife estimation
dt.1.fst$mean.fst
# 0.004459117

# CI
dt.1.fst$mean.fst+c(-1.96,1.96)*dt.1.fst$se.fst
# 0.002462227 0.006456007

## FST GENOME SCAN
dt.1.fst.scan <- computeFST(dt.1,sliding.window.size=10)

# Multilocus FST thresholds
dt.1.fst.scan.thresh <-  quantile(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst, probs = 0.99, na.rm = T)
dt.1.fst.scan.025.thres <- 0.25

# Manhattan plot prep
chrm.splitted <- strsplit(dt.1.fst.scan$sliding.windows.fst$Chr, split = "_")
chrm.splitted <- as.numeric(do.call(rbind.data.frame, chrm.splitted)[,2])

chrm.splitted.counts <- table(chrm.splitted)

chrm.colors=NULL
for (i in seq_along(chrm.splitted.counts))
{
  if((i %% 2) == 0) 
  {
    #print(paste(i,"is Even"))
    c = rep("black",chrm.splitted.counts[i])
  } else {
    #print(paste(i,"is Odd"))
    c = rep("grey",chrm.splitted.counts[i])
  }
  
  chrm.colors <- c(chrm.colors, c)
}

# Manhattan plot
par(mar=c(5,5,4,1)+.1)
plot(dt.1.fst.scan$sliding.windows.fst$CumulatedPosition/1e6,
     dt.1.fst.scan$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst",
     col=chrm.colors,pch=20, cex.lab=1.4, cex.axis=1.2)
abline(h=c(dt.1.fst.scan$FST, dt.1.fst.scan.thresh, dt.1.fst.scan.025.thres),
       lty=c(2,2,2), col=c("lightgreen", "blue", "red"), lwd=2)
text(x=270, y=c(dt.1.fst.scan$FST, dt.1.fst.scan.thresh, dt.1.fst.scan.025.thres)-0.02,
     c(expression(paste("Global ", italic(F)[ST])), 
       expression(paste(italic(F)[ST], " > 99% quantile")),
       expression(paste(italic(F)[ST], " > 0.25"))),
     cex=0.75, pos=3,col=c("lightgreen", "blue", "red")
  )

# TOP Multilocus FST
high.multilocus.fst <- data.frame(CHRM=dt.1.fst.scan$sliding.windows.fst$Chr[which(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)], 
                                  POS=dt.1.fst.scan$sliding.windows.fst$Position[which(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)], 
                                  FST= dt.1.fst.scan$sliding.windows.fst$MultiLocusFst[which(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)])

#write.table(high.multilocus.fst, file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/multilocus.fst.higher025.txt",
#            quote = F, col.names = F, row.names = F)

high.multilocus.fst.ordered <- high.multilocus.fst[order(high.multilocus.fst$CHRM, high.multilocus.fst$POS), ]

high.multilocus.fst.ordered.bed <- data.frame(CHROM=high.multilocus.fst.ordered$CHRM,
                                              Start=high.multilocus.fst.ordered$POS-1,
                                              Ends=high.multilocus.fst.ordered$POS)
write.table(high.multilocus.fst.ordered.bed,
            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst.bed", sep = "\t",
            quote = F, row.names = F, col.names = F)


## CHECK NUMBER OF SNPS WITHOUT NaN RESULTS
sum(!is.na(dt.1.fst.scan$snp.FST)) #19842
sum(!is.na(dt.1.fst.scan$sliding.windows.fst$MultiLocusFst)) #18207

## CHECK NUMBER OF SNPS WITH NO MISSING DATA
#Data set with all 21 pools
no.missing.all <- dt.1.imputedRefMLFreq[complete.cases(dt.1.imputedRefMLFreq), ]
dim(no.missing.all) #19842    21
#
##Data set with 19 pools
#no.missing.reduced <- reduced.dt.1.imputedRefMLFreq[complete.cases(reduced.dt.1.imputedRefMLFreq), ]
#dim(no.missing.reduced) #33313    19
#
##Data set with 16 pools (no PA samples)
#no.missing.reduced.2 <- dt.1.imputedRefMLFreq[complete.cases(dt.1.imputedRefMLFreq[,-c(5,10:12,21)]) , ]
#dim(no.missing.reduced.2) #82534    21
#
## COMPUTE PAIRWISE FST
dt.1.pfst <- compute.pairwiseFST(dt.1,verbose=FALSE)
heatmap(dt.1.pfst)

dt.1.pfst <- compute.pairwiseFST(dt.1,nsnp.per.bjack.block = 100, verbose=FALSE)
#heatmap(dt.1.pfst)
plot(dt.1.pfst, cex=0.4)

## COMPUTE THE F-STATISTICS 
dt.1.fstats<-compute.fstats(dt.1, nsnp.per.bjack.block = 100, verbose=FALSE)
dt.1.fstats@heterozygosities

## COMPUTE THETA FROM F-STATS HE 
barplot(dt.1.fstats@heterozygosities[,1]/(1-dt.1.fstats@heterozygosities[,1]))

## POOL MEAN HE
pool.meanHE <- apply(dt.1.imputedRefMLFreq, 2, meanHE)
names(pool.meanHE) <- poolnames
barplot(pool.meanHE)

biotype.vec <- c("B1","B4","B1","B4","B1","B1","B4","B4","B4","B1",
                 "B4","B4","B1","B1","B4","B4","B4","B1","B4","B4",
                 "B4")

boxplot(pool.meanHE ~ biotype.vec, 
        xlab="",ylab="Mean Heterozygosity",
        col=chrm.colors,pch=20, cex.lab=1.4, cex.axis=1.2)

## POOL THETA
pool.theta <- thetaHE(pool.meanHE)
barplot(pool.theta)

boxplot(pool.theta ~ biotype.vec,
        xlab="",ylab=expression(paste("Mean ", theta)),
        col=chrm.colors,pch=20, cex.lab=1.4, cex.axis=1.2)

## MEAN LOCUS HE
locus.meanHE <- apply(dt.1.imputedRefMLFreq, 1, meanLocusHE)

## PLOT INTRA-LOCUS FST AS A FUNCTION OF LOCUS MEAN HE
plot(locus.meanHE, dt.1.fst.scan$snp.FST,
     xlab=expression(paste("Locus ", H[E])), ylab=expression(paste("Locus ", F[ST])),
     cex.lab=1.4, cex.axis=1.2)
abline(h=c(dt.1.fst.scan$FST, quantile(dt.1.fst.scan$snp.FST, probs = 0.99, na.rm = T), 0.25),
       lty=2, col=c("lightgreen", "blue", "red"), lwd=2)
## NEED SIMULATION TO CALCULATE THE ENVELOPE

## PCA ON INPUTED ALLELE COUNT
#W<-scale(t(dt.1.imputedRefMLFreq), scale=TRUE) # inputed Reference allele frequency
W<-scale(t(dt.1.imputedRefMLCount[[1]]), scale=TRUE) # inputed Reference allele count
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

l1 <- 100*lambda[1]/sum(lambda) 
l2 <- 100*lambda[2]/sum(lambda) 

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=biotype.col,
     xlim = c(-0.35, 0.95), ylim = c(-0.35, 0.75),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"), cex=1.7, pch=19)
text(eig.vec[,1],eig.vec[,2], poolnames, pos=2 , cex = 0.6)
legend("topright", 
       legend = c("Biotype 1", "Biotype 4"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

## pooldata2genobaypass
pooldata2genobaypass(pooldata = dt.1,
                     writing.dir = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/",
                     prefix = "aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.complete50")


## IF need remove high missing samples, re-do below with the dataset from above, like:
#Data set 2: Min_cov=4; Max_cov=0.99; [start MAF=0.05 (from above); final MAF=0.05]; MAC=5; 19 pools; [start miss_fraction=0.50; final miss_fraction=0.30]





#### REDUCED DATASET
### Data set 2: Min_cov=4; Max_cov=0.99; [start MAF=0.01 (from above); final MAF=0.01]; MAC=5; 19 pools; [start miss_fraction=0.50; final miss_fraction=0.30]
####
#
## REMOVE 2 SAMPLES IDENTIFIED WITH > 30% MISSING DATA - KEEP PA POOLS; REMOVE 1 NW-BIO1, 1 WO-BIO4
## REMOVE SNPS WITH > 30% AFTER THE REMOVAL OF THE SAMPLES
## SUBSET READ COUNT DATA
### PERFORM ALL THE ANALYSES
#
#keep.pools <- c(1:4,6:20)
#nbr.pools.reduced <- dim(reduced.dt.1.imputedRefMLFreq)[2]
#
### REMOVE SNPS WITH > 30% AFTER THE REMOVAL OF THE SAMPLES
#length(dt.1.missing.snps.reduced)
#min(dt.1.missing.snps.reduced)
#max(dt.1.missing.snps.reduced)
#
#keep.snps <- which(dt.1.missing.snps.reduced <= 30)
#length(keep.snps) #254107
#
## How many were removed?
#length(dt.1.missing.snps.reduced) - length(keep.snps)
#
## CHECK MISSING % AFTER SNP REMOVAL
#reduced.dt.1.imputedRefMLFreq.filtered03 <- reduced.dt.1.imputedRefMLFreq[keep.snps, ] 
#dt.1.missing.snps.reduced.filtered03 <- (apply(reduced.dt.1.imputedRefMLFreq.filtered03, 1, function(x) sum(is.na(x)))/nbr.pools.reduced)*100
#
#length(dt.1.missing.snps.reduced.filtered03)
#min(dt.1.missing.snps.reduced.filtered03)
#max(dt.1.missing.snps.reduced.filtered03)
#
## SAVE TO COMPARE WITH BAYPASS ESTIMATES
##save(reduced.dt.1.imputedRefMLFreq.filtered03, file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/reduced.dt.1.imputedRefMLFreq.filtered03.RData")
#
### SUBSET READ COUNT DATA
#s.dt.1 <- pooldata.subset(
#                          pooldata=dt.1,
#                          pool.index = keep.pools,
#                          snp.index = keep.snps,
#                          min.cov.per.pool = -1,
#                          max.cov.per.pool = 1e+06,
#                          min.maf = -1,
#                          cov.qthres.per.pool = c(0, 1),
#                          return.snp.idx = FALSE,
#                          verbose = TRUE
#)
#
#s.dt.1
## * * * PoolData Object * * * 
## * Number of SNPs   =  254107 
## * Number of Pools  =  19
## Number of scaffolds:
#length(unique(sort(s.dt.1@snp.info$Chromosome))) #685
#
#### POPULATION GENOMICS ANALYSES
### COMPUTE THE GLOBAL FST
#s.dt.1.fst <- computeFST(s.dt.1, method = "Anova", nsnp.per.bjack.block = 10, verbose=FALSE)
#
## genome-wide estimate
#s.dt.1.fst$FST 
## 0.005548407
#
## mean Block-Jackknife estimation
#s.dt.1.fst$mean.fst 
## 0.005536966
#
## CI
#s.dt.1.fst$mean.fst+c(-1.96,1.96)*s.dt.1.fst$se.fst
## 0.004055994 0.007017939
#
### FST GENOME SCAN
#s.dt.1.fst.scan <- computeFST(s.dt.1,sliding.window.size=10)
#
#chrm.splitted.s <- strsplit(s.dt.1.fst.scan$sliding.windows.fst$Chr, split = "_")
#chrm.splitted.s <- as.numeric(do.call(rbind.data.frame, chrm.splitted.s)[,2])
#
#chrm.splitted.counts.s <- table(chrm.splitted.s)
#
#chrm.colors.s=NULL
#for (i in seq_along(chrm.splitted.counts.s))
#{
#  if((i %% 2) == 0) 
#  {
#    #print(paste(i,"is Even"))
#    c = rep("black",chrm.splitted.counts.s[i])
#  } else {
#    #print(paste(i,"is Odd"))
#    c = rep("grey",chrm.splitted.counts.s[i])
#  }
#  
#  chrm.colors.s <- c(chrm.colors.s, c)
#}
#
## Manhattan plot
#par(mar=c(5,5,4,1)+.1)
#plot(s.dt.1.fst.scan$sliding.windows.fst$CumulatedPosition/1e6,
#     s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst,
#     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst",
#     col=chrm.colors.s,pch=20, cex.lab=1.4, cex.axis=1.2)
#abline(h=s.dt.1.fst.scan$FST,lty=2, col="lightgreen", lwd=3)
#
## TOP Multilocus FST
#high.multi.fst.reduced <- data.frame(CHRM=s.dt.1.fst.scan$sliding.windows.fst$Chr[which(s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)], 
#                                     POS=s.dt.1.fst.scan$sliding.windows.fst$Position[which(s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)], 
#                                     FST= s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst[which(s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst > 0.25)])
#
#write.table(high.multi.fst, file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/higher025.multilocus.fst.reduced.txt",
#            quote = F, col.names = F, row.names = F)
#
#high.multi.fst.reduced.ordered <- high.multi.fst.reduced[order(high.multi.fst.reduced$CHRM, 
#                                                               high.multi.fst.reduced$POS), ]
#
#high.multi.fst.reduced.ordered.bed <- data.frame(CHROM=high.multi.fst.reduced.ordered$CHRM,
#                                                 Start=high.multi.fst.reduced.ordered$POS-1,
#                                                 Ends=high.multi.fst.reduced.ordered$POS)
#write.table(high.high.multi.fst.reduced.ordered.bed,
#            file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/higher025.multilocus.fst.reduced.bed", sep = "\t",
#            quote = F, row.names = F, col.names = F)
#
### CHECK NUMBER OF SNPS WITH NaN RESULTS
#sum(!is.na(s.dt.1.fst.scan$snp.FST)) #42566
#sum(!is.na(s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst)) #25832
#
##Data set with 19 pools
#dim(no.missing.reduced) #42608    19
#
### COMPUTE PAIRWISE FST
#s.dt.1.pfst <- compute.pairwiseFST(s.dt.1,verbose=FALSE)
#heatmap(s.dt.1.pfst)
## SAVE TO COMPARE WITH BAYPASS OMEGA MATRIX
##save(s.dt.1.pfst, file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/s.dt.1.pfst.RData")
#
#s.dt.1.pfst <- compute.pairwiseFST(s.dt.1,nsnp.per.bjack.block = 100, verbose=FALSE)
#plot(s.dt.1.pfst, cex=0.4)
#
### COMPUTE THE F-STATISTICS 
#s.dt.1.fstats<-compute.fstats(s.dt.1, nsnp.per.bjack.block = 100, verbose=FALSE)
#s.dt.1.fstats@heterozygosities
#
### COMPUTE THETA FROM F-STATS HE 
#barplot(s.dt.1.fstats@heterozygosities[,1]/(1-s.dt.1.fstats@heterozygosities[,1]))
#
### POOL MEAN HE
#s.pool.meanHE <- apply(reduced.dt.1.imputedRefMLFreq.filtered03, 2, meanHE)
#names(s.pool.meanHE) <- s.dt.1@poolnames
#barplot(s.pool.meanHE)
#
#biotype.vec.s <- c("B1","B4","B1","B4","B1","B4","B4","B4","B1",
#                  "B4","B4","B1","B1","B4","B4","B4","B1","B4","B4")
#
#boxplot(s.pool.meanHE ~ biotype.vec.s, 
#       xlab="",ylab="Mean Heterozygosity",
#       col=chrm.colors,pch=20, cex.lab=1.4, cex.axis=1.2)
#
### POOL THETA
#s.pool.theta <- thetaHE(s.pool.meanHE)
#barplot(s.pool.theta)
#
#boxplot(s.pool.theta ~ biotype.vec.s,
#        xlab="",ylab=expression(paste("Mean ", theta)),
#        col=chrm.colors,pch=20, cex.lab=1.4, cex.axis=1.2)
#
### MEAN LOCUS HE
#s.locus.meanHE <- apply(reduced.dt.1.imputedRefMLFreq.filtered03, 1, meanLocusHE)
#
### PLOT INTRA-LOCUS FST AS A FUNCTION OF LOCUS MEAN HE
#plot(s.locus.meanHE, s.dt.1.fst.scan$snp.FST,
#     xlab=expression(paste("Locus ", H[E])), ylab=expression(paste("Locus ", F[ST])),
#     cex.lab=1.4, cex.axis=1.2)
#abline(h=dt.1.fst.scan$FST,lty=2, col="lightgreen", lwd=3)
### NEED SIMULATION TO CALCULATE THE ENVELOPE
#
### PCA ON INPUTED AF
#s.W<-scale(t(reduced.dt.1.imputedRefMLFreq.filtered03), scale=TRUE) #centering
#s.W[is.na(s.W)]<-0
#cov.s.W<-cov(t(s.W))
#s.eig.result<-eigen(cov.s.W)
#s.eig.vec<-s.eig.result$vectors
#s.lambda<-s.eig.result$values
#
## PVE
#par(mar=c(5,5,4,1)+.1)
#plot(s.lambda/sum(s.lambda),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.2,
#     cex.lab=1.6, pch=19, col="black")
#lines(s.lambda/sum(s.lambda), col="red")
#
#s.l1 <- 100*s.lambda[1]/sum(s.lambda) 
#s.l2 <- 100*s.lambda[2]/sum(s.lambda) 
#
## PCA plot
#par(mar=c(5,5,4,1)+.1)
#plot(s.eig.vec[,1],s.eig.vec[,2], col=biotype.col[keep.pools],
#     xlim = c(-0.75, 0.35), ylim = c(-0.75, 0.55),
#     xlab=paste0("PC1 (", round(s.l1,2), "%)"), ylab=paste0("PC1 (", round(s.l2,2), "%)"), cex=1.7, pch=19)
#text(s.eig.vec[,1],s.eig.vec[,2], s.dt.1@poolnames, pos=2 , cex = 0.6)
#legend("bottomleft", 
#       legend = c("Biotype 1", "Biotype 4"), col = c("#247F00","#AB1A53"), 
#       pch = 19, bty = "n", cex = 1.2)
#abline(v=0,h=0,col="grey",lty=3)
#
### pooldata2genobaypass FROM REDUCED DATASET
#pooldata2genobaypass(pooldata = s.dt.1,
#                     writing.dir = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/",
#                     prefix = "aphidpool.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30")
#
##SAVE WORKSPACE
##load("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.workspace_28Jun21.RData")
