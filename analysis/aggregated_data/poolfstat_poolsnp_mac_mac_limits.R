##########################################################
#          RELATIONSHIP BETWEEN MAF AND MAC AND          # 
#   THE ACTUAL MAF THRESHOLD AFTER POOLFSTAT LOADING     #
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
setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()

# Remove last features
rm(list=ls())
ls()

# Load R libraries
library(poolfstat)
packageVersion("poolfstat") 

# Import auxiliary R functions  
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

### POOL'S INFORMATION
# Pool's size
poolsizes <- rep(10,21)

# Pool's names
poolnames <- c("MN-B1.1", "MN-B4.1",
               "ND-B1.1", "ND-B4.1", 
               "NW-B1.1", "NW-B1.2", "NW-B4.1", "NW-B4.2", "NW-B4.3", 
               "PA-B1.1", "PA-B4.1", "PA-B4.2", 
               "WI-B1.1", "WI-B1.2", "WI-B4.1", "WI-B4.2", "WI-B4.3", 
               "WO-B1.1", "WO-B4.1", "WO-B4.2", "WO-B4.3")

###
###
### ---- INVESTIGATE THE MINIMUM MAF AND MAC VALUES AFTER POOLFSTAT DATASET LOADING  ----
###
###

### Dataset 1: Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=5; 21 pools; miss_fraction=0.50
### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz

testcase.dt.1 <- vcf2pooldata(
  vcf.file = "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz",
  poolsizes = poolsizes,
  poolnames = poolnames,
  min.cov.per.pool = -1,
  min.rc = 1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.001,
  remove.indels = FALSE,
  nlines.per.readblock = 1e+06
)

testcase.dt.1
# * * * PoolData Object * * * 
# * Number of SNPs   = 344524 
# * Number of Pools  =  21 
# Number of scaffolds:

### CHECK THE FREQUENCY OF EACH MAC CLASS FROM THE ORIGINAL READ COUNT DATA
testcase.ALTallele.count <- (testcase.dt.1@readcoverage - testcase.dt.1@refallele.readcount)

testcase.overall.ALTallele.count <- apply(testcase.ALTallele.count, 1, sum)
testcase.overall.REFallele.count <- apply(testcase.dt.1@refallele.readcount, 1, sum)

# MAC
overall.MAC <- pmin(testcase.overall.REFallele.count, testcase.overall.ALTallele.count)

# MAF
testcase.imputedMLCount <- imputedRefMLCount(testcase.dt.1)
testcase.imputedRefMLFreq <- testcase.imputedMLCount[[1]]/testcase.imputedMLCount[[3]]

testcase.overall.ALTallele.freq <- apply(testcase.imputedRefMLFreq, 1,  function(x) 1-mean(x, na.rm = T))
testcase.overall.REFallele.freq <- apply(testcase.imputedRefMLFreq, 1,  function(x) mean(x, na.rm = T))

overall.MAF <- pmin(testcase.overall.REFallele.freq, testcase.overall.ALTallele.freq)

### HOW MAF CHANGE WITH MAC

macs <- seq(5, 100, 5)
mafs.macs.table <- NULL
for (i in seq_along(macs))
{
  # REF/ALT Alleles
  a <- testcase.overall.REFallele.freq[overall.MAC >= macs[i]] 
  b <- testcase.overall.ALTallele.freq[overall.MAC >= macs[i]] 
  ref.min <- min(a); ref.max <- max(a)
  alt.min <- min(b); alt.max <- max(b)
  
  # MAF
  p <- pmin(a,b)
  maf <- min(p)
  
  #OUTPUT
  res <- data.frame(MAC=macs[i], MAF=maf, REF_MIN=ref.min, REF_MAX=ref.max, ALT_MIN=alt.min, ALT_MAX=alt.max)
  mafs.macs.table <- rbind(mafs.macs.table,res)
}

### MAF AS A FUNCTION OF MAC
plot(mafs.macs.table$MAC, mafs.macs.table$MAF,  type = "l", lty = 1, col="black", xlab="MAC", ylab="MAF", ylim=c(0,0.15))
lines(mafs.macs.table$MAC, mafs.macs.table$REF_MIN,  type = "l", lty = 1, col="blue")
#lines(mafs.macs.table$MAC, mafs.macs.table$REF_MAX,  type = "l", lty = 2, col="blue")
lines(mafs.macs.table$MAC, mafs.macs.table$ALT_MIN,  type = "l", lty = 1, col="red")
#lines(mafs.macs.table$MAC, mafs.macs.table$ALT_MAX,  type = "l", lty = 2, col="red")

# MAC         MAF    REF_MIN
#   5 0.004761905 0.01904762 
#  10 0.014285714 0.03529412 
#  15 0.019047619 0.04210526 
#  20 0.022222222 0.05555556 
#  25 0.023809524 0.08500000 
#  30 0.023809524 0.09000000 
#  35 0.023809524 0.09000000 
#  40 0.023809524 0.09000000 
#  45 0.023809524 0.09005848 
#  50 0.023809524 0.09473684 
# 100 0.061904762 0.10333333 

## Savepoint
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.macmaflimits.workspace22Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/poolfstat.poolsnp.macmaflimits.workspace22Jul21.RData")

############## END SCRIPT