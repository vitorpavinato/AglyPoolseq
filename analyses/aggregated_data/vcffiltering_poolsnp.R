##########################################################
#         PREPARE POOLSNP VCF FILE FOR POPOOLATION       # 
#                                                        #
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


###
###
### ---- VCF 1: MAF=0.001 and MAC=2 ----
###
###

#### PoolSNPs - Popoolation Summary Statistics Dataset:
#### Dataset diversity 1: 13-Ago-2021; Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=2; 21 pools; miss_fraction=0.50
#### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.2.13Ago2021.vcf.gz


### open GDS file
genofile <- seqOpen("vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.2.13Ago2021.ann.gds")

seqResetFilter(genofile)
# of selected samples: 21
# of selected variants: 627,825

### get subsample of data to work on
snps.dt.vcf1 <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                           variant.id=seqGetData(genofile, "variant.id"),
                           nAlleles=seqNumAllele(genofile),
                           missing=seqMissing(genofile, .progress=T))


### Calculate sequencing error rate (approximate with % > 2-allelic SNPs)
total_snps.vcf1 <- length(snps.dt.vcf1$variant.id)
trialle_snps.vcf1 <- length(snps.dt.vcf1$variant.id[snps.dt.vcf1$nAlleles > 2])

(trialle_snps.vcf1/total_snps.vcf1)*100
# 9.079441 % It is a lot

### choose number of alleles
snps.dt.vcf1 <- snps.dt.vcf1[nAlleles==2]

### Apply filter
seqSetFilter(genofile, variant.id=snps.dt.vcf1$variant.id)

### Get SNPs chrom and pos and export .bed file
snps_ordered4bed.vcf1 <- snps.dt.vcf1[order(snps.dt.vcf1$chr, snps.dt.vcf1$pos), -c(3,4,5)]

## Create a BED file -like table
snps_ordered.bed.vcf1 <- data.frame(CHROM=snps_ordered4bed.vcf1$chr,
                                    Start=as.integer(snps_ordered4bed.vcf1$pos-1),
                                    Ends=as.integer(snps_ordered4bed.vcf1$pos))
# Save to a file
write.table(snps_ordered.bed.vcf1,
            file="results/aggregated_data/minmaxcov_4_99/popoolation_sumstats/snps_ordered4mpileup_filtering.001.2.13Ago2021.bed", sep = "\t",
            quote = F, row.names = F, col.names = F)

###
###
### ---- VCF 2: MAF=0.001 and MAC=1 ----
###
###

#### PoolSNPs - Popoolation Summary Statistics Dataset:
#### Dataset diversity 2: 13-Ago-2021; Min_cov=4; Max_cov=0.99; MAF=0.001; MAC=1; 21 pools; miss_fraction=0.50
#### vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.1.13Ago2021.vcf.gz


### open GDS file
genofile <- seqOpen("vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.001.1.13Ago2021.ann.gds")

seqResetFilter(genofile)
# of selected samples: 21
# of selected variants: 6,595,229

### get subsample of data to work on
snps.dt.vcf2 <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                           variant.id=seqGetData(genofile, "variant.id"),
                           nAlleles=seqNumAllele(genofile),
                           missing=seqMissing(genofile, .progress=T))


### Calculate sequencing error rate (approximate with % > 2-allelic SNPs)
total_snps.vcf2 <- length(snps.dt.vcf2$variant.id)
trialle_snps.vcf2 <- length(snps.dt.vcf2$variant.id[snps.dt.vcf2$nAlleles > 2])

(trialle_snps.vcf2/total_snps.vcf2)*100
# 4.028382 % It is a lot

### choose number of alleles
snps.dt.vcf2 <- snps.dt.vcf2[nAlleles==2]

### Apply filter
seqSetFilter(genofile, variant.id=snps.dt.vcf2$variant.id)

### Get SNPs chrom and pos and export .bed file
snps_ordered4bed.vcf2 <- snps.dt.vcf2[order(snps.dt.vcf2$chr, snps.dt.vcf2$pos), -c(3,4,5)]

## Create a BED file -like table
snps_ordered.bed.vcf2 <- data.frame(CHROM=snps_ordered4bed.vcf2$chr,
                                    Start=as.integer(snps_ordered4bed.vcf2$pos-1),
                                    Ends=as.integer(snps_ordered4bed.vcf2$pos))
# Save to a file
write.table(snps_ordered.bed.vcf2,
            file="results/aggregated_data/minmaxcov_4_99/popoolation_sumstats/snps_ordered4mpileup_filtering.001.1.13Ago2021.bed", sep = "\t",
            quote = F, row.names = F, col.names = F)
