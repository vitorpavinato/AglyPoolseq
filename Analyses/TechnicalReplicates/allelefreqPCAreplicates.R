#############################################################
# Access sequencing consistency across technical replicates #
#     Allele frequency calculation from vcf and PCA         #
#############################################################

# Each library, that corresponds to pooled sequecing of 5 individuals that share the same
# barcode, were sequenced 4-5x in different days. The libraries were the same, but they were 
# loaded and sequenced in different runs.

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES

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
          vcf.file = "aphidpool.all.PoolSNP.001.88.08Apr2021.vcf",
          poolsizes = poolsizes,
          poolnames = poolnames,
          min.cov.per.pool = -1,
          min.rc = 1,
          max.cov.per.pool = 1e+06,
          min.maf = 0.01,
          remove.indels = FALSE,
          nlines.per.readblock = 1e+06
          )


# Pairwise FST between replicates
pfst <- computePairwiseFSTmatrix(
                        dt,
                        method = "Anova",
                        min.cov.per.pool = -1,
                        max.cov.per.pool = 1e+06,
                        min.maf = -1,
                        output.snp.values = FALSE
                        )

heatmap(pfst$PairwiseFSTmatrix)

#coverage = read.table("SNP_coverage.tab")
coverage = dt@readcoverage

#counts = read.table("SNP_refcounts.tab")
counts = dt@refallele.readcount

# Sum counts per race 
sum_coverage <- cbind(coverage[,1] + coverage[,2], coverage[,3] + coverage[,4], coverage[,5] + coverage[,6])
sum_coverage <- cbind(apply(coverage, MARGIN = 1, sum))

#sum_counts <- cbind(counts[,1] + counts[,2], counts[,3] + counts[,4], counts[,5] + counts[,6])
sum_counts <- cbind(apply(counts, MARGIN = 1, sum))

nbr.gene.copies = 5 # 60 diploid individuals per pool
nbr.loci = nrow(sum_coverage)

# For the first pool
reads.1.allele.1 <- sum_counts[,1]
reads.ALF.allele.2 <- sum_coverage[,1] - sum_counts[,1]
total.counts.ALF <- sum_counts[,1]

ALF.allele.counts.1 <- vector(length = nbr.loci, mode = "integer")
ALF.allele.counts.2 <- vector(length = nbr.loci, mode = "integer")

likelihood.ALF <- matrix(nrow = nbr.loci, ncol = (nbr.gene.copies + 1))

for (i in 1:nbr.loci) {
  for (j in 1:(nbr.gene.copies + 1)) {
    likelihood.ALF[i,j] <- dbinom(x = reads.ALF.allele.1[i], size = total.counts.ALF[i],
                                  prob = (j - 1) / nbr.gene.copies, log = FALSE)
  }
  ALF.allele.counts.1[i] <- (which(likelihood.ALF[i,] == max(likelihood.ALF[i,])) - 1)
  ALF.allele.counts.2[i] <- (nbr.gene.copies - ALF.allele.counts.1[i])
}

write.table(file="ML_count_ALF", ALF.allele.counts.1, row.names=F, quote=F, col.names=F)