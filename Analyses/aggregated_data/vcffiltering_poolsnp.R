##########################################################
#         PREPARE POOLSNP VCF FILE FOR GRENEDALF         # 
#              REMOVE NON-BIALLELIC SNPS                 #
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
library(data.table)
library(vcfR)
library(SeqArray)

###
###
### ---- PART 1: UPLOAD AND FILTER SNP FILE ----
###
###

### open GDS file
genofile <- seqOpen("vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.gds")

### get target populations
samps <- fread("DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv")

### get subsample of data to work on
seqResetFilter(genofile)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

### choose number of alleles
snps.dt <- snps.dt[nAlleles==2]


seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id[6])

### import CHRM and POS information of SNPs retained by Poolfstat
poofstat.retained.snps <- readLines("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/snpinfo_aggregated_names_pools.txt")

### Take only Poolfstat retained SNPs
snps.dt.poolfstat <- snps.dt[paste0(snps.dt$chr, "_", snps.dt$pos) %in% poofstat.retained.snps, ]

### select sites
seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt.poolfstat$variant.id)

###
###
### ---- PART 2: EXPORT VCF FILE ----
###
###

seqGDS2VCF(genofile, 
           vcf.fn= "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.filtered.vcf", 
           info.var=NULL, fmt.var=NULL, verbose=TRUE)



