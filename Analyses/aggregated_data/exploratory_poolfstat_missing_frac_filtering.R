#############################################################
#            pF_ST, gF_ST and genome scan                   #
#           exploratory missing filtering                   #
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

## Check if samples with >30% missing data affect the final number of SNPs

####
## 21 Pools
####

poolsizes.21 <- rep(5,21)
poolnames.21 <- c("MN_BIO1_S1", "MN_BIO4_S1",
               "ND_BIO1_S1", "ND_BIO4_S1", 
               "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3", 
               "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2", 
               "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3", 
               "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3")

## Data set 1: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=10; 21 pools; miss_fraction=0.50
dt1 <- vcf2pooldata(
                    vcf.file = "vcf/aggregated_data/minmaxcov_4_99/pools_all/miss_frac05/aphidpool.PoolSeq.PoolSNP.001.10.14Jun2021.vcf.gz",
                    poolsizes = poolsizes.21,
                    poolnames = poolnames.21,
                    min.cov.per.pool = -1,
                    min.rc = 1,
                    max.cov.per.pool = 1e+06,
                    min.maf = 0.001,
                    remove.indels = FALSE,
                    nlines.per.readblock = 1e+06
)

## 293830 out of 303635 SNPs
(293830/303635)*100 #96.77079 of SNPs processed

dt1@readcoverage[dt1@readcoverage==0] <- NA
missing.pop.1 <- 100*(colSums(is.na(dt1@readcoverage))/dim(dt1@readcoverage)[1])
names(missing.pop.1)<-poolnames.21
#MN_BIO1_S1 MN_BIO4_S1 ND_BIO1_S1 ND_BIO4_S1 NW_BIO1_S1 NW_BIO1_S2 NW_BIO4_S1 NW_BIO4_S2 NW_BIO4_S3 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 
#3.041214   7.529524   7.294354  16.095021  64.251438  16.615730   3.320968  16.177381  12.218970  67.313753  30.930130  39.562672 
#WI_BIO1_S1 WI_BIO1_S2 WI_BIO4_S1 WI_BIO4_S2 WI_BIO4_S3 WO_BIO1_S1 WO_BIO4_S1 WO_BIO4_S2 WO_BIO4_S3 
#8.080182   2.873771   7.378076   6.248171  12.749549  22.561685  13.559201  10.258313  53.213082 

barplot(missing.pop.1, ylim = c(0,100))
abline(h=30, col="red")

missing.overall.1 <- 100*(sum(is.na(dt1@readcoverage))/(dim(dt1@readcoverage)[1] * dim(dt1@readcoverage)[2]))
# 20.06063

pfst.1 <- computePairwiseFSTmatrix(dt1,
                                 method = "Anova",
                                 min.cov.per.pool = -1,
                                 max.cov.per.pool = 1e+06,
                                 min.maf = -1,
                                 output.snp.values = FALSE)

heatmap(pfst.1$PairwiseFSTmatrix)

fst.1 <- computeFST(dt1, method = "Anova", snp.index = NA)
plot(fst.1$snp.FST, ylim = c(-1.5, 1))

# Relationship between fst and missing %
dt1@readcoverage[dt1@readcoverage==0] <- NA
missing.snp.1 <- 100*(rowSums(is.na(dt1@readcoverage))/dim(dt1@readcoverage)[2])
max(missing.snp.1) #47.61905

plot(fst.1$snp.FST~missing.snp.1, ylim = c(-1.5, 1))

## Data set 2: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=10; 21 pools; miss_fraction=0.30
dt2 <- vcf2pooldata(
                    vcf.file = "vcf/aggregated_data/minmaxcov_4_99/pools_all/miss_frac03/aphidpool.PoolSeq.PoolSNP.001.10.14Jun2021.vcf.gz",
                    poolsizes = poolsizes.21,
                    poolnames = poolnames.21,
                    min.cov.per.pool = -1,
                    min.rc = 1,
                    max.cov.per.pool = 1e+06,
                    min.maf = 0.001,
                    remove.indels = FALSE,
                    nlines.per.readblock = 1e+06
)

## 216819 out of 226248 SNPs
(216819/226248)*100 #95.83245 of SNPs processed

dt2@readcoverage[dt2@readcoverage==0] <- NA
missing.pop.2 <- 100*(colSums(is.na(dt2@readcoverage))/dim(dt2@readcoverage)[1])
names(missing.pop.2)<-poolnames.21
#MN_BIO1_S1 MN_BIO4_S1 ND_BIO1_S1 ND_BIO4_S1 NW_BIO1_S1 NW_BIO1_S2 NW_BIO4_S1 NW_BIO4_S2 NW_BIO4_S3 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 
#2.078692   2.874748   2.667202   7.867853  53.821851   6.966640   1.785821   5.971340   4.658725  57.944184  19.020012  23.533915 
#WI_BIO1_S1 WI_BIO1_S2 WI_BIO4_S1 WI_BIO4_S2 WI_BIO4_S3 WO_BIO1_S1 WO_BIO4_S1 WO_BIO4_S2 WO_BIO4_S3 
#5.935827   1.717100   3.048626   2.555127   5.951508  10.498157   8.795355   3.754284  40.508443

barplot(missing.pop.2, ylim = c(0,100))
abline(h=30, col="red")

missing.overall.2 <- 100*(sum(is.na(dt2@readcoverage))/(dim(dt2@readcoverage)[1] * dim(dt2@readcoverage)[2]))
# 12.95026

pfst.2 <- computePairwiseFSTmatrix(dt2,
                                   method = "Anova",
                                   min.cov.per.pool = -1,
                                   max.cov.per.pool = 1e+06,
                                   min.maf = -1,
                                   output.snp.values = FALSE)

heatmap(pfst.2$PairwiseFSTmatrix)

fst.2 <- computeFST(dt2, method = "Anova", snp.index = NA)
plot(fst.2$snp.FST, ylim = c(-1.5, 1))

# Relationship between fst and missing %
dt2@readcoverage[dt2@readcoverage==0] <- NA
missing.snp.2 <- 100*(rowSums(is.na(dt2@readcoverage))/dim(dt2@readcoverage)[2])
max(missing.snp.2) #28.57143

plot(fst.2$snp.FST~missing.snp.2, ylim = c(-1.5, 1))


####
## 16 Pools
####

## Data set 3: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=10; 16 pools; miss_fraction=0.50
poolsizes.16 <- rep(5,16)
poolnames.16 <- c("MN_BIO1_S1", "MN_BIO4_S1",
               "ND_BIO1_S1", "ND_BIO4_S1", 
               "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3", 
               "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3", 
               "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2")

dt3 <- vcf2pooldata(
                    vcf.file = "vcf/aggregated_data/minmaxcov_4_99/pools_filtered/miss_frac05/aphidpool.PoolSeq.PoolSNP.001.10.14Jun2021.vcf.gz",
                    poolsizes = poolsizes.16,
                    poolnames = poolnames.16,
                    min.cov.per.pool = -1,
                    min.rc = 1,
                    max.cov.per.pool = 1e+06,
                    min.maf = 0.001,
                    remove.indels = FALSE,
                    nlines.per.readblock = 1e+06
)

## 323114 out of 330502 SNPs
(323114/330502)*100 #97.76461 of SNPs processed

dt3@readcoverage[dt3@readcoverage==0] <- NA
missing.pop.3 <- 100*(colSums(is.na(dt3@readcoverage))/dim(dt3@readcoverage)[1])
names(missing.pop.3)<-poolnames.16
#MN_BIO1_S1 MN_BIO4_S1 ND_BIO1_S1 ND_BIO4_S1 NW_BIO1_S2 NW_BIO4_S1 NW_BIO4_S2 NW_BIO4_S3 WI_BIO1_S1 WI_BIO1_S2 WI_BIO4_S1 WI_BIO4_S2 
#3.289551  11.672041  12.025477  22.249113  23.818838   4.996688  23.658214  18.534944   9.261128   3.690029  11.373695   9.650464 
#WI_BIO4_S3 WO_BIO1_S1 WO_BIO4_S1 WO_BIO4_S2 
#17.965176  30.421152  16.691323  16.361099 

barplot(missing.pop.3, ylim = c(0,100))
abline(h=30, col="red")

missing.overall.3 <- 100*(sum(is.na(dt3@readcoverage))/(dim(dt3@readcoverage)[1] * dim(dt3@readcoverage)[2]))
# 14.72868

pfst.3 <- computePairwiseFSTmatrix(dt3,
                                   method = "Anova",
                                   min.cov.per.pool = -1,
                                   max.cov.per.pool = 1e+06,
                                   min.maf = -1,
                                   output.snp.values = FALSE)

heatmap(pfst.3$PairwiseFSTmatrix)

fst.3 <- computeFST(dt3, method = "Anova", snp.index = NA)
plot(fst.3$snp.FST,ylim = c(-1.5, 1))

# Relationship between fst and missing %
dt3@readcoverage[dt3@readcoverage==0] <- NA
missing.snp.3 <- 100*(rowSums(is.na(dt3@readcoverage))/dim(dt3@readcoverage)[2])
max(missing.snp.3) # 50


plot(fst.3$snp.FST~missing.snp.3, ylim = c(-1.5, 1))

## Data set 4: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=10; 16 pools; miss_fraction=0.30

dt4 <- vcf2pooldata(
                    vcf.file = "vcf/aggregated_data/minmaxcov_4_99/pools_filtered/miss_frac03/aphidpool.PoolSeq.PoolSNP.001.10.14Jun2021.vcf.gz",
                    poolsizes = poolsizes.16,
                    poolnames = poolnames.16,
                    min.cov.per.pool = -1,
                    min.rc = 1,
                    max.cov.per.pool = 1e+06,
                    min.maf = 0.001,
                    remove.indels = FALSE,
                    nlines.per.readblock = 1e+06
)

## 251193 out of 258337 SNPs
(251193/258337)*100 #97.23462 of SNPs processed

dt4@readcoverage[dt4@readcoverage==0] <- NA
missing.pop.4 <- 100*(colSums(is.na(dt4@readcoverage))/dim(dt4@readcoverage)[1])
names(missing.pop.4)<-poolnames.16
#MN_BIO1_S1 MN_BIO4_S1 ND_BIO1_S1 ND_BIO4_S1 NW_BIO1_S2 NW_BIO4_S1 NW_BIO4_S2 NW_BIO4_S3 WI_BIO1_S1 WI_BIO1_S2 WI_BIO4_S1 WI_BIO4_S2 
#1.774333   4.921315   4.418515  12.483230  12.507514   1.852759  11.943804   8.677790   6.076602   1.536667   4.839705   3.797478 
#WI_BIO4_S3 WO_BIO1_S1 WO_BIO4_S1 WO_BIO4_S2 
#9.492701  18.344062  10.912724   6.802339 

barplot(missing.pop.4, ylim = c(0,100))
abline(h=30, col="red")

missing.overall.4 <- 100*(sum(is.na(dt4@readcoverage))/(dim(dt4@readcoverage)[1] * dim(dt4@readcoverage)[2]))
# 7.523846

pfst.4 <- computePairwiseFSTmatrix(dt4,
                                   method = "Anova",
                                   min.cov.per.pool = -1,
                                   max.cov.per.pool = 1e+06,
                                   min.maf = -1,
                                   output.snp.values = FALSE)

heatmap(pfst.4$PairwiseFSTmatrix)

fst.4 <- computeFST(dt4, method = "Anova", snp.index = NA)
plot(fst.4$snp.FST, ylim = c(-1.5, 1))

# Relationship between fst and missing %
dt4@readcoverage[dt4@readcoverage==0] <- NA
missing.snp.4 <- 100*(rowSums(is.na(dt4@readcoverage))/dim(dt4@readcoverage)[2])
max(missing.snp.4) # 25

plot(fst.4$snp.FST~missing.snp.4, ylim = c(-1.5, 1))


####
## Subset data: take the dataset with 21 pools and filter out pools with missing data
####

## SubData set 1: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=10; 21 pools; miss_fraction=0.50; 16 pools
sub.dt.1 <- pooldata.subset(
                            dt1,
                            pool.index = c(1,2,3,4,6,7,8,9,13:20),
                            snp.index = 1:dt1@nsnp,
                            min.cov.per.pool = -1,
                            max.cov.per.pool = 1e+06,
                            min.maf = -1,
                            cov.qthres.per.pool = c(0, 1)
                            )


sub.dt.1@readcoverage[sub.dt.1@readcoverage==0] <- NA
missing.pop.sub.dt.1 <- 100*(colSums(is.na(sub.dt.1@readcoverage))/dim(sub.dt.1@readcoverage)[1])
names(missing.pop.sub.dt.1)<-poolnames.16
#MN_BIO1_S1 MN_BIO4_S1 ND_BIO1_S1 ND_BIO4_S1 NW_BIO1_S2 NW_BIO4_S1 NW_BIO4_S2 NW_BIO4_S3 WI_BIO1_S1 WI_BIO1_S2 WI_BIO4_S1 WI_BIO4_S2 
#3.041214   7.529524   7.294354  16.095021  16.615730   3.320968  16.177381  12.218970   8.080182   2.873771   7.378076   6.248171 
#WI_BIO4_S3 WO_BIO1_S1 WO_BIO4_S1 WO_BIO4_S2 
#12.749549  22.561685  13.559201  10.258313

barplot(missing.pop.sub.dt.1, ylim = c(0,100))
abline(h=30, col="red")

missing.overall.sub.dt.1 <- 100*(sum(is.na(sub.dt.1@readcoverage))/(dim(sub.dt.1@readcoverage)[1] * dim(sub.dt.1@readcoverage)[2]))
# 10.37513

pfst.sub.dt.1 <- computePairwiseFSTmatrix(sub.dt.1,
                                          method = "Anova",
                                          min.cov.per.pool = -1,
                                          max.cov.per.pool = 1e+06,
                                          min.maf = -1,
                                          output.snp.values = FALSE)
        
heatmap(pfst.sub.dt.1$PairwiseFSTmatrix)


fst.sub.dt.1 <- computeFST(sub.dt.1, method = "Anova", snp.index = NA)

plot(fst.sub.dt.1$snp.FST, ylim = c(-1.5, 1))

# Relationship between fst and missing %
sub.dt.1@readcoverage[sub.dt.1@readcoverage==0] <- NA
missing.snp.sub.dt.1 <- 100*(rowSums(is.na(sub.dt.1@readcoverage))/dim(sub.dt.1@readcoverage)[2])
max(missing.snp.sub.dt.1) # 62.5

length(which(missing.snp.sub.dt.1 > 50))

plot(fst.sub.dt.1$snp.FST~missing.snp.sub.dt.1, ylim = c(-1.5, 1))

## SubData set 2: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=10; 21 pools; miss_fraction=0.30; 16 pools
sub.dt.2 <- pooldata.subset(
                            dt2,
                            pool.index = c(1,2,3,4,6,7,8,9,13:20),
                            snp.index = 1:dt2@nsnp,
                            min.cov.per.pool = -1,
                            max.cov.per.pool = 1e+06,
                            min.maf = -1,
                            cov.qthres.per.pool = c(0, 1)
)


sub.dt.2@readcoverage[sub.dt.2@readcoverage==0] <- NA
missing.pop.sub.dt.2 <- 100*(colSums(is.na(sub.dt.2@readcoverage))/dim(sub.dt.2@readcoverage)[1])
names(missing.pop.sub.dt.2)<-poolnames.16
#MN_BIO1_S1 MN_BIO4_S1 ND_BIO1_S1 ND_BIO4_S1 NW_BIO1_S2 NW_BIO4_S1 NW_BIO4_S2 NW_BIO4_S3 WI_BIO1_S1 WI_BIO1_S2 
#2.078692   2.874748   2.667202   7.867853   6.966640   1.785821   5.971340   4.658725   5.935827   1.717100 
#WI_BIO4_S1 WI_BIO4_S2 WI_BIO4_S3 WO_BIO1_S1 WO_BIO4_S1 WO_BIO4_S2 
#3.048626   2.555127   5.951508  10.498157   8.795355   3.754284

barplot(missing.pop.sub.dt.2, ylim = c(0,100))
abline(h=30, col="red")

missing.overall.sub.dt.2 <- 100*(sum(is.na(sub.dt.2@readcoverage))/(dim(sub.dt.2@readcoverage)[1] * dim(sub.dt.2@readcoverage)[2]))
# 4.820438

pfst.sub.dt.2 <- computePairwiseFSTmatrix(sub.dt.2,
                                          method = "Anova",
                                          min.cov.per.pool = -1,
                                          max.cov.per.pool = 1e+06,
                                          min.maf = -1,
                                          output.snp.values = FALSE)

heatmap(pfst.sub.dt.2$PairwiseFSTmatrix)


fst.sub.dt.2 <- computeFST(sub.dt.2, method = "Anova", snp.index = NA)

plot(fst.sub.dt.2$snp.FST, ylim = c(-1.5, 1))


# Relationship between fst and missing %
sub.dt.2@readcoverage[sub.dt.2@readcoverage==0] <- NA
missing.snp.sub.dt.2 <- 100*(rowSums(is.na(sub.dt.2@readcoverage))/dim(sub.dt.2@readcoverage)[2])
max(missing.snp.sub.dt.2) # 37.5

length(which(missing.snp.sub.dt.2 > 30))

plot(fst.sub.dt.2$snp.FST~missing.snp.sub.dt.2, ylim = c(-1.5, 1))

## Overall missing
overall.missing <- c(missing.overall.1,missing.overall.2,
                         missing.overall.3,missing.overall.4,
                         missing.overall.sub.dt.1, missing.overall.sub.dt.2)
names(overall.missing) <- c("p21m05", "p21m03", "p16m05", "p16m03", "p21m05s16", "p21m03s16")

barplot(overall.missing, ylim = c(0,100))

