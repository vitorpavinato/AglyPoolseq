#############################################################
# Access sequencing consistency across technical replicates #
#     Allele frequency calculation from vcf and PCA         #
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

###
###
### ---- MAXIMUM LIKELIHOOD ESTIMATION FROM ALLELE COUNTS ----
###
###

## ADD CITATION HERE

# First load vcf using poolfsta function vcf2pooldata from poolfstat package

poolsizes <- rep(5,87)
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

dt <- vcf2pooldata(
          vcf.file = "vcf/aphidpool.all.PoolSNP.001.50.15Apr2021.vcf.gz",
          poolsizes = poolsizes,
          poolnames = poolnames,
          min.cov.per.pool = -1,
          min.rc = 1,
          max.cov.per.pool = 1e+06,
          min.maf = 0.01,
          remove.indels = FALSE,
          nlines.per.readblock = 1e+06
          )


# Check the Pairwise FST between replicates
#pfst <- computePairwiseFSTmatrix(
#                        dt,
#                        method = "Anova",
#                        min.cov.per.pool = -1,
#                        max.cov.per.pool = 1e+06,
#                       min.maf = -1,
#                        output.snp.values = FALSE
#                        )
#
#heatmap(pfst$PairwiseFSTmatrix)

# sample 20,000 random SNPs
selectedsnps <- sample(1:dim(dt@snp.info)[1], 20000, replace = F)

# Pooling all the information for the allele frequency calculation
coverage = dt@readcoverage[selectedsnps, ]
counts = dt@refallele.readcount[selectedsnps, ]
snpInfo = dt@snp.info[selectedsnps, ]

write.table(snpInfo, file = "results/analysis_technical_reps/snps.technical.rep.20000.snps.txt", sep = "\t",
            quote = F, col.names = F, row.names = F)

nbr.gene.copies = 10
nbr.loci = length(selectedsnps)

# For the first pool
#reads.ALF.allele.1 <- counts[,1]
#reads.ALF.allele.2 <- coverage[,1] - counts[,1]
#total.counts.ALF <- coverage[,1]

#ALF.allele.counts.1 <- vector(length = nbr.loci, mode = "integer")
#ALF.allele.counts.2 <- vector(length = nbr.loci, mode = "integer")

#likelihood.ALF <- matrix(nrow = nbr.loci, ncol = (nbr.gene.copies + 1))

#for (i in 1:nbr.loci) {
#  if (total.counts.ALF[i] == 0){
#    ALF.allele.counts.1[i] <- NA
#    ALF.allele.counts.2[i] <- NA
#  } else {
#    for (j in 1:(nbr.gene.copies + 1)) {
#      likelihood.ALF[i,j] <- dbinom(x = reads.ALF.allele.1[i], size = total.counts.ALF[i],
#                                    prob = (j - 1) / nbr.gene.copies, log = FALSE)
#    } # end of for i,j
#    ALF.allele.counts.1[i] <- (which(likelihood.ALF[i,] == max(likelihood.ALF[i,])) - 1)
#    ALF.allele.counts.2[i] <- (nbr.gene.copies - ALF.allele.counts.1[i])
#  }# end o if
#} # end of for i
#
#write.table(file="ML_count_ALF", ALF.allele.counts.1, row.names=F, quote=F, col.names=F)

# for all technical replicates
refcount.matrix = matrix(nrow=nbr.loci, ncol = 87)
for (k in 1:87){
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

legend.names <- c("MN_BIO1_S1", "MN_BIO4_S1",
                  "ND_BIO1_S1", "ND_BIO4_S1",
                  "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3",
                  "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2",
                  "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3",
                  "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3",
                  paste0("r_",1:5)
                  )

legend.pch <- c(rep(19, 21), c(8, 15, 17, 18, 19))


W<-scale(t(refcount.matrix), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
pdf(file = "results/analysis_technical_reps/pve.technical.rep.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")
dev.off()

100*lambda[1]/sum(lambda) # 3.159421
100*lambda[2]/sum(lambda) # 2.674815

# PCA plot
pdf(file = "results/analysis_technical_reps/pca.technical.rep.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=run.colors,
     xlim = c(-0.40, 0.25), ylim = c(-0.25, 0.25),
     xlab="PC1 (3.15%)", ylab="PC2 (2.67%)", cex=1.5, pch=lib.pch)
legend("topleft", 
       legend = legend.names, col = legend.colors, 
       pch = legend.pch, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Remove the _AddRun_ from the allele frequency matrix (5th technical replication)
W<-scale(t(refcount.matrix[,-c(5,30,59)]), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
pdf(file = "results/analysis_technical_reps/pve.technical.rep.20000.snps.wtAddRun.pdf")
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")
dev.off()

100*lambda[1]/sum(lambda) # 3.174587
100*lambda[2]/sum(lambda) # 2.676656

# PCA plot
pdf(file = "results/analysis_technical_reps/pca.technical.rep.20000.snps.wtAddRun.pdf")
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=run.colors[-c(5,30,59)],
     xlim = c(-0.40, 0.40), ylim = c(-0.25, 0.25),
     xlab="PC1 (3.17%)", ylab="PC2 (2.67%)", cex=1.5, pch=lib.pch[-c(5,30,59)])
legend("topleft", 
       legend = legend.names[-c(26)], col = legend.colors[-c(26)], 
       pch = legend.pch[-c(26)], bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Aggregate read counts data for each pool
## and repeat the allele frequency estimates
## and PCA.

comb.coverage <- cbind(coverage[,1]+coverage[,2]+coverage[,3]+coverage[,4]+coverage[,5], #MN_B1
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

comb.counts <- cbind(counts[,1] +counts[,2] +counts[,3] +counts[,4] +counts[,5], #MN_B1
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


# for all technical replicates
refcount.matrix = matrix(nrow=nbr.loci, ncol = 21)
for (k in 1:21){
  reads.ALF.allele.1 <- comb.counts[,k]
  reads.ALF.allele.2 <- comb.coverage[,k] - comb.counts[,k]
  total.counts.ALF <- comb.coverage[,k]
  
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

legend.colors.c <- c("#0000FF", "#FF0000", #MN
                     "#8282FF", "#FF8282", #ND
                     "#10B717", "#64C366", "#FD63CB", "#F98AD1", "#F4A8D6", #NW
                     "#B78560", "#9E5C00", "#7D3200", #PA
                     "#B2B2FF", "#D5D5FF", "#FFB2B2", "#FFD5D5", "#FFF2F2", #WI
                     "#90CE91", "#EFC0DB", "#E9D2DF", "#E4DFE2" #WO
                     )

legend.names.c <- c("MN_BIO1_S1", "MN_BIO4_S1",
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

100*lambda[1]/sum(lambda) # 11.0807
100*lambda[2]/sum(lambda) # 10.02013

# PCA plot
pdf(file = "results/analysis_technical_reps/pca.technical.combined.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=legend.colors.c,
     xlim = c(-1.1, 0.2), ylim = c(-0.8, 0.3),
     xlab="PC1 (11.08%)", ylab="PC2 (10.02%)", cex=1.5, pch=19, cex.lab=1.6)
legend("topleft", 
       legend = legend.names.c, col = legend.colors.c, 
       pch = legend.pch, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

# Without the outlier population
W<-scale(t(refcount.matrix[,-c(5,21)]), scale=TRUE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
pdf(file = "results/analysis_technical_reps/pve.technical.combined.popoutlier.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda/sum(lambda), col="red")
dev.off()

100*lambda[1]/sum(lambda) # 11.60427
100*lambda[2]/sum(lambda) # 8.669536

# PCA plot
pdf(file = "results/analysis_technical_reps/pca.technical.combined.popoutlier.20000.snps.pdf")
par(mar=c(5,5,4,1)+.1)
plot(eig.vec[,1],eig.vec[,2], col=legend.colors.c[-c(5,21)],
     xlim = c(-1.1, 0.3), ylim = c(-0.7, 0.4),
     xlab="PC1 (11.60%)", ylab="PC2 (8.66%)", cex=1.5, pch=19, cex.lab=1.6)
legend("topleft", 
       legend = legend.names.c[-c(5,21)], col = legend.colors.c[-c(5,21)], 
       pch = legend.pch, bty = "n", cex = 0.6)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

# NA's count
dt.af.20000 <- refcount.matrix
colnames(dt.af.20000) <- legend.names.c
colSums(is.na(dt.af.20000))

100*(colSums(is.na(dt.af.20000))/20000)

hist(rowSums(is.na(dt.af.20000)))
100*(max(rowSums(is.na(dt.af.20000)))/21)
100*(min(rowSums(is.na(dt.af.20000)))/21)
100*(median(rowSums(is.na(dt.af.20000)))/21)

# Coverage ranking
cov.20000 <- comb.coverage
colnames(cov.20000) <- legend.names.c

hist(100*(colSums(cov.20000)/sum(cov.20000)))
mean(100*(colSums(cov.20000)/sum(cov.20000)))
quantile(100*(colSums(cov.20000)/sum(cov.20000)), c(0.025, 0.5, 0.975))

# SFS
hist(dt.af.20000[,1])
hist(dt.af.20000[,5])
hist(dt.af.20000[,10])
hist(dt.af.20000[,11])
hist(dt.af.20000[,21])

overall.pi <- function(x)
{
  p <- x
  q <- (1 - x)
  pq <- 2*p*q
  pi <- sum(pq, na.rm = T)/sum(!is.na(x))
  return(pi)
}

pop.pi <- apply(dt.af.20000, MARGIN = 2, overall.pi)
barplot(pop.pi)

theta <- pop.pi/(1-pop.pi)
