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
library(SeqArray)



### open GDS file
genofile <- seqOpen("vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.gds")

seqResetFilter(genofile)

### get subsample of data to work on
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

### choose number of alleles
snps.dt <- snps.dt[nAlleles==2]

### Apply filter
seqSetFilter(genofile, variant.id=snps.dt$variant.id)

### import CHRM and POS information of SNPs retained by Poolfstat
poofstat.retained.snps <- readLines("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/snpinfo_aggregated_names_pools.txt")

### Take only Poolfstat retained SNPs
snps.dt.poolfstat <- snps.dt[paste0(snps.dt$chr, "_", snps.dt$pos) %in% poofstat.retained.snps, ]

### select sites
seqSetFilter(genofile, variant.id=snps.dt.poolfstat$variant.id)

###
###
### ---- PART 2: EXPORT VCF FILE ----
###
###

seqGDS2VCF(genofile, 
           vcf.fn= "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.filtered.vcf", 
           info.var=NULL, fmt.var=NULL, verbose=TRUE)

#### PoolSNPs - Popoolation Summary Statistics Dataset:
#### Dataset 2: 03-Ago-2021; Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=1; 21 pools; miss_fraction=0.50
#### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.1.03Ago2021.vcf.gz

###
###
### ---- PART 1: UPLOAD AND FILTER SNP FILE ----
###
###

### open GDS file
genofile <- seqOpen("vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.01.1.04Ago2021.ann.gds")

seqResetFilter(genofile)

### get the data to work on
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile, .progress=T))

dim(snps.dt)
# R=1,503,409; C=5

###
###
### ---- PART 2: EXPORT VCF AND BED FILE ----
###
###

# No VCF since the file wasn't filtered
#seqGDS2VCF(genofile, 
#           vcf.fn= "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.01.1.04Ago2021.ann.filtered.vcf", 
#           info.var=NULL, fmt.var=NULL, verbose=TRUE)

# Get all SNPs information
snps_ordered4bed <- snps.dt[order(snps.dt$chr, snps.dt$pos), -c(3,4,5)]

# Create a BED file -like table
snps_ordered.bed <- data.frame(CHROM=snps_ordered4bed$chr,
                               Start=as.integer(snps_ordered4bed$pos-1),
                               Ends=as.integer(snps_ordered4bed$pos))
# Save to a file
#write.table(snps_ordered.bed,
#            file="results/aggregated_data/minmaxcov_4_99/popoolation_sumstats/snps_ordered4mpileup_filtering.01.1.04Ago2021.bed", sep = "\t",
#            quote = F, row.names = F, col.names = F)

