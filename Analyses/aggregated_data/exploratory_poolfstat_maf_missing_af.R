#############################################################
#             Poolfsta: FST and Pairwise FST                #
#                                                           #
#############################################################

# Each library have 5 pooled individuals and they were sequenced 4-5 times (4-5 days.)
# They were the same library, with same barcode, but the sequencing was splitted in 4-5 days. 

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES

setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
## Remove last features
rm(list=ls())
ls()

## Load libraries
library(poolfstat)

# First load vcf using poolfsta function vcf2pooldata from poolfstat package

poolsizes <- rep(5,21)
poolnames <- c("MN_BIO1_S1", "MN_BIO4_S1",
               "ND_BIO1_S1", "ND_BIO4_S1", 
               "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3", 
               "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2", 
               "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3", 
               "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3")

dt <- vcf2pooldata(
  vcf.file = "vcf/aphidpool.PoolSeq.PoolSNP.001.5.10May2021.vcf",
  poolsizes = poolsizes,
  poolnames = poolnames,
  min.cov.per.pool = -1,
  min.rc = 1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  remove.indels = FALSE,
  nlines.per.readblock = 1e+06
)

# Check if different MAF produce:
# Different number of SNPs;
# Different proportion of missing data: within pool and overall;
# Which pool has more missing data and if it varies with MAF;

# MAF 0.001 (same used for PoolSNPs calling);
# Number of SNPs: 265644 out of 283912

dt@readcoverage[dt@readcoverage==0] <- NA
missing.within.001 <- 100*(colSums(is.na(dt@readcoverage))/dim(dt@readcoverage)[1])
names(missing.within.001)<- poolnames

par(mar=c(5,5,4,1)+.1)
barplot(missing.within.001,
        beside=T, ylim=c(0,100), ylab="Missing data (%)", las=2,
        col = adjustcolor(c(rep("grey",4),"red",rep("grey",4), rep("red",3),rep("grey",8), "red"), alpha.f = 0.6),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 0.7)
abline(h=30,col="black", lty="dashed")

100*(sum(is.na(dt@readcoverage))/(dim(dt@readcoverage)[1]*dim(dt@readcoverage)[2])) #24.89381


# MAF 0.01 (same used for PoolSNPs calling);
# Number of SNPs: 247806 out of 283912

dt@readcoverage[dt@readcoverage==0] <- NA
missing.within.01 <- 100*(colSums(is.na(dt@readcoverage))/dim(dt@readcoverage)[1])
names(missing.within.01)<- poolnames

par(mar=c(5,5,4,1)+.1)
barplot(missing.within.01,
        beside=T, ylim=c(0,100), ylab="Missing data (%)", las=2,
        col = adjustcolor(c(rep("grey",4),"red",rep("grey",4), rep("red",3),rep("grey",8), "red"), alpha.f = 0.6),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 0.7)
abline(h=30,col="black", lty="dashed")

100*(sum(is.na(dt@readcoverage))/(dim(dt@readcoverage)[1]*dim(dt@readcoverage)[2])) #26.04862

# MAF 0.05 (same used for PoolSNPs calling);
# Number of SNPs: 196660 out of 283912

dt@readcoverage[dt@readcoverage==0] <- NA
missing.within.05 <- 100*(colSums(is.na(dt@readcoverage))/dim(dt@readcoverage)[1])
names(missing.within.05)<- poolnames

par(mar=c(5,5,4,1)+.1)
barplot(missing.within.05,
        beside=T, ylim=c(0,100), ylab="Missing data (%)", las=2,
        col = adjustcolor(c(rep("grey",4),"red",rep("grey",4), rep("red",3),rep("grey",8), "red"), alpha.f = 0.6),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 0.7)
abline(h=30,col="black", lty="dashed")

100*(sum(is.na(dt@readcoverage))/(dim(dt@readcoverage)[1]*dim(dt@readcoverage)[2])) #28.4615


# Check the Pairwise FST between replicates
pfst <- computePairwiseFSTmatrix(dt,
                                 method = "Anova",
                                 min.cov.per.pool = -1,
                                 max.cov.per.pool = 1e+06,
                                 min.maf = -1,
                                 output.snp.values = FALSE)

heatmap(pfst$PairwiseFSTmatrix)

fst <- computeFST(dt, method = "Anova", snp.index = NA)

plot(fst$snp.FST)

dim(dt@snp.info)
100*(colSums(dt@readcoverage==0)/132251)
(dt@poolnames)

# Pooling all the information for the allele frequency calculation
coverage = dt@readcoverage
counts = dt@refallele.readcount

nbr.gene.copies = 10
nbr.loci = dim(dt@snp.info)[1]
npools = 21

# for all technical replicates
refcount.matrix = matrix(nrow=nbr.loci, ncol = npools)
for (k in 1:npools){
  reads.ALF.allele.1 <- counts[,k]
  reads.ALF.allele.2 <- coverage[,k] - counts[,k]
  total.counts.ALF <- coverage[,k]
  
  ALF.allele.counts.1 <- vector(length = nbr.loci, mode = "integer")
  ALF.allele.counts.2 <- vector(length = nbr.loci, mode = "integer")
  
  likelihood.ALF <- matrix(nrow = nbr.loci, ncol = (nbr.gene.copies + 1))
  
  for (i in 1:nbr.loci) {
    if (total.counts.ALF[i] == 0){
      ALF.allele.counts.1[i] <- NA
      ALF.allele.counts.2[i] <- NA
    } else {
      for (j in 1:(nbr.gene.copies + 1)) {
        likelihood.ALF[i,j] <- dbinom(x = reads.ALF.allele.1[i], size = total.counts.ALF[i],
                                      prob = (j - 1) / nbr.gene.copies, log = FALSE)
      } # end of for i,j
      ALF.allele.counts.1[i] <- (which(likelihood.ALF[i,] == max(likelihood.ALF[i,])) - 1)
      ALF.allele.counts.2[i] <- (nbr.gene.copies - ALF.allele.counts.1[i])
    }# end o if
  } # end of for i
  
  refcount.matrix[,k] <- ALF.allele.counts.1/10
}

###
###
### ---- PCA ON ML ALLELE FREQUENCIES ----
###
###


## Price et al 2010 - vcov of matrix individuals x markers

pool.colors <- c("#0000FF", "#FF0000", #MN
                     "#8282FF", "#FF8282", #ND
                     "#10B717", "#64C366", "#FD63CB", "#F98AD1", "#F4A8D6", #NW
                     "#B78560", "#9E5C00", "#7D3200", #PA
                     "#B2B2FF", "#D5D5FF", "#FFB2B2", "#FFD5D5", "#FFF2F2", #WI
                     "#90CE91", "#EFC0DB", "#E9D2DF", "#E4DFE2" #WO
)

pool.names <- c("MN_BIO1_S1", "MN_BIO4_S1",
                    "ND_BIO1_S1", "ND_BIO4_S1",
                    "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3",
                    "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2",
                    "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3",
                    "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3"
)

W<-scale(t(refcount.matrix), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
pdf(file = "results/analysis_technical_reps/pve.technical.combined.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")
dev.off()

100*lambda[1]/sum(lambda) # 7.43222
100*lambda[2]/sum(lambda) # 7.053729

# PCA plot
pdf(file = "results/analysis_technical_reps/pca.technical.combined.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=pool.colors,
     #xlim = c(-1.1, 0.2), ylim = c(-0.8, 0.3),
     xlab="PC1 (11.16%)", ylab="PC2 (10.09%)", cex=1.5, pch=19, cex.lab=1.6)
legend("topleft", 
       legend = pool.names, col = pool.colors, 
       pch = 19, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

##mtDNA

# First load vcf using poolfsta function vcf2pooldata from poolfstat package

poolsizes <- rep(5,21)
poolnames <- c("MN_BIO1_S1", "MN_BIO4_S1",
               "ND_BIO1_S1", "ND_BIO4_S1", 
               "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3", 
               "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2", 
               "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3", 
               "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3")

dt <- vcf2pooldata(
  vcf.file = "vcf/aphidpool.mitochondrion_genome.PoolSeq.PoolSNP.01.12.11May2021.vcf.gz",
  poolsizes = poolsizes,
  poolnames = poolnames,
  min.cov.per.pool = -1,
  min.rc = 2,
  max.cov.per.pool = 1e+06,
  min.maf = 0.01,
  remove.indels = FALSE,
  nlines.per.readblock = 1e+06
)

# Pooling all the information for the allele frequency calculation
coverage = dt@readcoverage
counts = dt@refallele.readcount

nbr.gene.copies = 5
nbr.loci = dim(dt@snp.info)[1]
npools = 21

# for all technical replicates
refcount.matrix = matrix(nrow=nbr.loci, ncol = npools)
for (k in 1:npools){
  reads.ALF.allele.1 <- counts[,k]
  reads.ALF.allele.2 <- coverage[,k] - counts[,k]
  total.counts.ALF <- coverage[,k]
  
  ALF.allele.counts.1 <- vector(length = nbr.loci, mode = "integer")
  ALF.allele.counts.2 <- vector(length = nbr.loci, mode = "integer")
  
  likelihood.ALF <- matrix(nrow = nbr.loci, ncol = (nbr.gene.copies + 1))
  
  for (i in 1:nbr.loci) {
    if (total.counts.ALF[i] == 0){
      ALF.allele.counts.1[i] <- NA
      ALF.allele.counts.2[i] <- NA
    } else {
      for (j in 1:(nbr.gene.copies + 1)) {
        likelihood.ALF[i,j] <- dbinom(x = reads.ALF.allele.1[i], size = total.counts.ALF[i],
                                      prob = (j - 1) / nbr.gene.copies, log = FALSE)
      } # end of for i,j
      ALF.allele.counts.1[i] <- (which(likelihood.ALF[i,] == max(likelihood.ALF[i,])) - 1)
      ALF.allele.counts.2[i] <- (nbr.gene.copies - ALF.allele.counts.1[i])
    }# end o if
  } # end of for i
  
  refcount.matrix[,k] <- ALF.allele.counts.1/5
}

###
###
### ---- PCA ON ML ALLELE FREQUENCIES ----
###
###


## Price et al 2010 - vcov of matrix individuals x markers

pool.colors <- c("#0000FF", "#FF0000", #MN
                 "#8282FF", "#FF8282", #ND
                 "#10B717", "#64C366", "#FD63CB", "#F98AD1", "#F4A8D6", #NW
                 "#B78560", "#9E5C00", "#7D3200", #PA
                 "#B2B2FF", "#D5D5FF", "#FFB2B2", "#FFD5D5", "#FFF2F2", #WI
                 "#90CE91", "#EFC0DB", "#E9D2DF", "#E4DFE2" #WO
)

pool.names <- c("MN_BIO1_S1", "MN_BIO4_S1",
                "ND_BIO1_S1", "ND_BIO4_S1",
                "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3",
                "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2",
                "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3",
                "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3"
)

W<-scale(t(refcount.matrix), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
pdf(file = "results/analysis_technical_reps/pve.technical.combined.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")
dev.off()

100*lambda[1]/sum(lambda) # 7.43222
100*lambda[2]/sum(lambda) # 7.053729

# PCA plot
pdf(file = "results/analysis_technical_reps/pca.technical.combined.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=pool.colors,
     #xlim = c(-1.1, 0.2), ylim = c(-0.8, 0.3),
     xlab="PC1 (11.16%)", ylab="PC2 (10.09%)", cex=1.5, pch=19, cex.lab=1.6)
legend("topleft", 
       legend = pool.names, col = pool.colors, 
       pch = 19, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()