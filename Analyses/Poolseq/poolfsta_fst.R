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
  vcf.file = "vcf/aphidpool.PoolSeq.PoolSNP.001.5.10May2021.vcf.gz",
  poolsizes = poolsizes,
  poolnames = poolnames,
  min.cov.per.pool = 5,
  min.rc = 5,
  max.cov.per.pool = 1e+06,
  min.maf = 0.01,
  remove.indels = FALSE,
  nlines.per.readblock = 1e+06
)

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
100*(colSums(dt@readcoverage==0)/4584)
(dt@poolnames)
