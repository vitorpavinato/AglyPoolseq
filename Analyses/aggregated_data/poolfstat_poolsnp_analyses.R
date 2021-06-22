#############################################################
#         AF, HE, Theta_pi, pF_ST, gF_ST and genome scan    #
#                  poolfstat package                        #
#############################################################

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

### COMPLETE DATASET
## Data set 1: Min_cov=4; Max_cov=0.99; [start MAF=0.001; final MAF=0.01]; MAC=5; 21 pools; miss_fraction=0.50
###

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

dt.1 <- vcf2pooldata(
                     vcf.file = "vcf/aggregated_data/minmaxcov_4_99/pools_all/miss_frac05/aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz",
                     poolsizes = poolsizes,
                     poolnames = poolnames,
                     min.cov.per.pool = -1,
                     min.rc = 1,
                     max.cov.per.pool = 1e+06,
                     min.maf = 0.01,
                     remove.indels = FALSE,
                     nlines.per.readblock = 1e+06
)

dt.1
# * * * PoolData Object * * * 
# * Number of SNPs   = 321056 
# * Number of Pools  =  21 
# Number of scaffolds:
length(unique(sort(dt.1@snp.info$Chromosome))) #708

## COMPUTE MAXIMUM LIKELIHOOD IMPUTED SAMPLE ALLELE COUNT
dt.1.refallelecount <- imputedRefMLCount(dt.1)

dt.1.refalllefreq <- dt.1.refallelecount/dt.1@poolsizes[1]

# Sanity check if the MAF correspond to the MAF used in filtering
hist(apply(dt.1.refalllefreq == 1, 1, function(x) sum(x, na.rm = T)))
table(apply(dt.1.refalllefreq == 1, 1, function(x) sum(x, na.rm = T)))
#quasi.fixed <- which(apply(dt.1.refalllefreq == 1, 1, function(x) sum(x, na.rm = T)) == 20) 

high.than01.af <- which(apply(dt.1.refalllefreq, 1, function(x) 1-mean(x, na.rm = T)) >= 0.01) 
length(high.than01.af)

min(apply(dt.1.refalllefreq[high.than01.af,], 1,  function(x) 1-mean(x, na.rm = T)))
max(apply(dt.1.refalllefreq[high.than01.af,], 1,  function(x) 1-mean(x, na.rm = T)))

## KEEP LOCI WITH AF >= 0.01
# SUBSET READ COUNT DATA
dt.1 <- pooldata.subset(
                        pooldata=dt.1,
                        snp.index = high.than01.af,
                        min.cov.per.pool = -1,
                        max.cov.per.pool = 1e+06,
                        min.maf = -1,
                        cov.qthres.per.pool = c(0, 1),
                        return.snp.idx = FALSE,
                        verbose = TRUE
)

dt.1.refallelecount <- dt.1.refallelecount[high.than01.af, ]
dt.1.refalllefreq <- dt.1.refalllefreq[high.than01.af,]

## COMPUTE THE GLOBAL FST
dt.1.fst <- computeFST(dt.1, method = "Anova", nsnp.per.bjack.block = 10, verbose=FALSE)

# genome-wide estimate
dt.1.fst$FST 
# 0.002508262

# mean Block-Jackknife estimation
dt.1.fst$mean.fst
# 0.002532572

# CI
dt.1.fst$mean.fst+c(-1.96,1.96)*dt.1.fst$se.fst
# 0.0007563193 0.0043088252

## FST GENOME SCAN
dt.1.fst.scan <- computeFST(dt.1,sliding.window.size=10)

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
abline(h=dt.1.fst.scan$FST,lty=2, col="lightgreen", lwd=3)


## COMPUTE PAIRWISE FST
dt.1.pfst <- compute.pairwiseFST(dt.1,verbose=FALSE)
heatmap(dt.1.pfst)

dt.1.pfst <- compute.pairwiseFST(dt.1,nsnp.per.bjack.block = 100, verbose=FALSE)
plot(dt.1.pfst, cex=0.5)

## COMPUTE THE F-STATISTICS 
dt.1.fstats<-compute.fstats(dt.1, nsnp.per.bjack.block = 100, verbose=FALSE)
dt.1.fstats@heterozygosities

## COMPUTE THETA FROM F-STATS HE 
barplot(dt.1.fstats@heterozygosities[,1]/(1-dt.1.fstats@heterozygosities[,1]))

## POOL MEAN HE
pool.meanHE <- apply(dt.1.refalllefreq, 2, meanHE)
names(pool.meanHE) <- poolnames
barplot(pool.meanHE)

## POOL THETA
pool.theta <- thetaHE(pool.meanHE)
barplot(pool.theta)

## MEAN LOCUS HE
locus.meanHE <- apply(dt.1.refalllefreq, 1, meanLocusHE)

## PLOT INTRA-LOCUS FST AS A FUNCTION OF LOCUS MEAN HE
plot(locus.meanHE, dt.1.fst.scan$snp.FST,
     xlab="Locus HE",ylab="Locus Fst")
abline(h=dt.1.fst.scan$FST,lty=2, col="lightgreen")
## NEED SIMULATION TO CALCULATE THE ENVELOPE

## PCA ON INPUTED AF
W<-scale(t(dt.1.refalllefreq), scale=TRUE) #centering
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
     xlim = c(-0.65, 0.4), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1,2), "%)"), ylab=paste0("PC1 (", round(l2,2), "%)"), cex=1.7, pch=19)
text(eig.vec[,1],eig.vec[,2], poolnames, pos=2 , cex = 0.6)
legend("topright", 
       legend = c("Biotype 1", "Biotype 4"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)


### pooldata2genobaypass
# Data set 1: Min_cov=4; Max_cov=0.99; [start MAF=0.001; final MAF=0.01]; MAC=5; 21 pools; miss_fraction=0.50

pooldata2genobaypass(pooldata = dt.1,
                     writing.dir = "data/aggregated_data/genobaypass/",
                     prefix = "aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.missfrac05.allpools")


pooldata2genoselestim(pooldata = dt.1,
                      writing.dir = "data/aggregated_data/genoselestim//",
                      prefix = "aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.missfrac05.allpools")

save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pop_genomics/poolfstat.poolsnp.workspace.RData")
load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/pop_genomics/poolfstat.poolsnp.workspace.RData")

### REDUCED DATASET
## Data set 2: Min_cov=4; Max_cov=0.99; [start MAF=0.01 (from above); final MAF=0.01]; MAC=10; 21 pools; miss_fraction=0.30
###

# REMOVE 5 SAMPLES IDENTIFIED WITH > 30% MISSING DATA
# REMOVE SNPS WITH > 30% AFTER THE REMOVAL OF THE SAMPLES
# RAW SNPS WERE ALREADY FILTERED FOR MAF=0.01
# SUBSET READ COUNT DATA
# PERFORM ALL THE ANALYSES

missing.p.pools <- (apply(refcount.matrix, 2, function(x) sum(is.na(x)))/nbr.loci)*100
names(missing.p.pools) <- poolnames

# REMOVE 5 SAMPLES IDENTIFIED WITH > 30% MISSING DATA
reduced.refcount.matrix <- refcount.matrix[,!(missing.p.pools > 30)]
keep.pools <- which(missing.p.pools <= 30)

nbr.pools.reduced <- dim(reduced.refcount.matrix)[2]

# REMOVE SNPS WITH > 50% AFTER THE REMOVAL OF THE SAMPLES
missing.p.loci <- (apply(reduced.refcount.matrix, 1, function(x) sum(is.na(x)))/nbr.pools.reduced)*100
length(missing.p.loci)
min(missing.p.loci)
max(missing.p.loci)

keep.snps <- which(missing.p.loci <= 50)
length(keep.snps)

# How many were removed?
length(missing.p.loci) - length(keep.snps)

reduced.refcount.matrix.03.filtered <- reduced.refcount.matrix[keep.snps, ] 

missing.p.loci.filtered <- (apply(reduced.refcount.matrix.03.filtered, 1, function(x) sum(is.na(x)))/nbr.pools.reduced)*100
length(missing.p.loci.filtered)
min(missing.p.loci.filtered)
max(missing.p.loci.filtered)

# SUBSET READ COUNT DATA
s.dt.1 <- pooldata.subset(
                          pooldata=dt.1,
                          pool.index = keep.pools,
                          snp.index = keep.snps,
                          min.cov.per.pool = -1,
                          max.cov.per.pool = 1e+06,
                          min.maf = -1,
                          cov.qthres.per.pool = c(0, 1),
                          return.snp.idx = FALSE,
                          verbose = TRUE
)

s.dt.1
# * * * PoolData Object * * * 
# * Number of SNPs   =  253361 
# * Number of Pools  =  16
# Number of scaffolds:
length(unique(sort(s.dt.1@snp.info$Chromosome))) #690

## PERFORM ALL THE ANALYSES 

## COMPUTE THE GLOBAL FST
s.dt.1.fst <- computeFST(s.dt.1, method = "Anova", nsnp.per.bjack.block = 100, verbose=FALSE)

# genome-wide estimate
s.dt.1.fst$FST 
# -0.02009295

# mean Block-Jackknife estimation
s.dt.1.fst$mean.fst 
# -0.01920418

# CI
s.dt.1.fst$mean.fst+c(-1.96,1.96)*s.dt.1.fst$se.fst
# -0.02227349 -0.01613487

## FST GENOME SCAN
s.dt.1.fst.scan <- computeFST(s.dt.1,sliding.window.size=10)

chrm.splitted.s <- strsplit(s.dt.1.fst.scan$sliding.windows.fst$Chr, split = "_")
chrm.splitted.s <- as.numeric(do.call(rbind.data.frame, chrm.splitted.s)[,2])

chrm.splitted.counts.s <- table(chrm.splitted.s)

chrm.colors.s=NULL
for (i in seq_along(chrm.splitted.counts.s))
{
  if((i %% 2) == 0) 
  {
    #print(paste(i,"is Even"))
    c = rep("black",chrm.splitted.counts.s[i])
  } else {
    #print(paste(i,"is Odd"))
    c = rep("grey",chrm.splitted.counts.s[i])
  }
  
  chrm.colors.s <- c(chrm.colors.s, c)
}

# Manhattan plot
par(mar=c(5,5,4,1)+.1)
plot(s.dt.1.fst.scan$sliding.windows.fst$CumulatedPosition/1e6,
     s.dt.1.fst.scan$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst",
     col=chrm.colors.s,pch=20, cex.lab=1.4, cex.axis=1.2)
abline(h=s.dt.1.fst.scan$FST,lty=2, col="lightgreen", lwd=3)

## COMPUTE PAIRWISE FST
s.dt.1.pfst <- compute.pairwiseFST(s.dt.1,verbose=FALSE)
heatmap(s.dt.1.pfst)

s.dt.1.pfst <- compute.pairwiseFST(s.dt.1,nsnp.per.bjack.block = 100, verbose=FALSE)
plot(s.dt.1.pfst, cex=0.4)

## COMPUTE THE F-STATISTICS 
s.dt.1.fstats<-compute.fstats(s.dt.1, nsnp.per.bjack.block = 100, verbose=FALSE)
s.dt.1.fstats@heterozygosities

## COMPUTE THETA FROM F-STATS HE 
barplot(s.dt.1.fstats@heterozygosities[,1]/(1-s.dt.1.fstats@heterozygosities[,1]))

## REDUCED DATASET: reduced.refcount.matrix.03.filtered

# POOL MEAN HE
s.pool.meanHE <- apply(reduced.refcount.matrix.03.filtered, 2, meanHE)
names(s.pool.meanHE) <- s.dt.1@poolnames
barplot(s.pool.meanHE)

# POOL THETA
s.pool.theta <- thetaHE(s.pool.meanHE)
barplot(s.pool.theta)

# MEAN LOCUS HE
s.locus.meanHE <- apply(reduced.refcount.matrix.03.filtered, 1, meanLocusHE)

## PLOT INTRA-LOCUS FST AS A FUNCTION OF LOCUS MEAN HE
plot(s.locus.meanHE, s.dt.1.fst.scan$snp.FST,
     xlab="Locus HE",ylab="Locus Fst")
abline(h=dt.1.fst.scan$FST,lty=2, col="lightgreen")


## PCA ON INPUTED AF
s.W<-scale(t(reduced.refcount.matrix.03.filtered), scale=TRUE) #centering
s.W[is.na(s.W)]<-0
cov.s.W<-cov(t(s.W))
s.eig.result<-eigen(cov.s.W)
s.eig.vec<-s.eig.result$vectors
s.lambda<-s.eig.result$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(s.lambda/sum(s.lambda),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(s.lambda/sum(s.lambda), col="red")

s.l1 <- 100*s.lambda[1]/sum(s.lambda) # 9.704683
s.l2 <- 100*s.lambda[2]/sum(s.lambda) # 9.187907

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(s.eig.vec[,1],s.eig.vec[,2], col=biotype.col[keep.pools],
     #xlim = c(-0.35, 0.95), ylim = c(-0.95, 0.35),
     xlab=paste0("PC1 (", round(s.l1,2), "%)"), ylab=paste0("PC1 (", round(s.l2,2), "%)"), cex=1.7, pch=19)
text(s.eig.vec[,1],s.eig.vec[,2], s.dt.1@poolnames, pos=2 , cex = 0.6)
legend("topright", 
       legend = c("Biotype 1", "Biotype 4"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

### pooldata2genobaypass FROM REDUCED DATASET

pooldata2genobaypass(pooldata = s.dt.1,
                     writing.dir = "data/aggregated_data/genobaypass/",
                     prefix = "aphidpool.PoolSeq.PoolSNP.001.10.21Jun2021.vcf.misspools.missfrac03")


pooldata2genoselestim(pooldata = s.dt.1,
                     writing.dir = "data/aggregated_data/genoselestim//",
                     prefix = "aphidpool.PoolSeq.PoolSNP.001.10.21Jun2021.vcf.misspools.missfrac03")

