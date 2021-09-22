##########################################################
#            BASIC SEQUENCING SUMMARY STATISTICS         # 
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
library(poolfstat)
library(foreach)
library(data.table)

# Import auxiliary R functions  
#source("DEST-AglyPoolseq/analyses/aggregated_data/aux_func.R")

#ALPHA=0.75


###
###
### --- TECHNICAL REPLICATES ---
###
###

poolsizes <- rep(10,87)
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

dt.repl <- vcf2pooldata(
  vcf.file = "vcf/technical_replicates/minmaxcov_3_99/aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz",
  poolsizes = poolsizes,
  poolnames = poolnames,
  min.cov.per.pool = -1,
  min.rc = 1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  remove.indels = FALSE,
  nlines.per.readblock = 1e+06
)

dt.repl
# * * * PoolData Object * * * 
# * Number of SNPs   =  40105 
# * Number of Pools  =  87 

nbr.gene.copies = 10
nbr.loci = dt.repl@nsnp

### --- Average Read Coverage and % of missing data ---

dp <- dt.repl@readcoverage

dp[dp==0] <- NA

seqStats <- data.frame(sampleId = dt.repl@poolnames, 
                       AveReadDepth = colMeans(dp, na.rm = T), 
                       propMissing=colSums(is.na(dp))/nbr.loci)

### --- Proportion of PCR duplicates and Number of Mapped reads ---

### PCR Duplicates
fns <- system("ls /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/*_B*_S*_*/*duplicates_report.txt", intern=T)

pcr <- foreach(fn=fns)%do%{
  #fn <- fns[1]
  data.table(pcrDup=fread(fn, skip="LIBRARY", nrows=1)$PERCENT_DUPLICATION, sampleId=tstrsplit(fn, "/")[[8]])
}
pcr <- rbindlist(pcr)

seqStats <- merge(seqStats, pcr, by="sampleId")

### Mapped Reads
fns <- system("ls /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/*_B*_S*_*/*flagstat_agly.txt", intern=T)

mpr <- foreach(fn=fns)%do%{
  #fn <- fns[1]
  data.table(mappedReads=fread(fn, nrow=1)$V1/1e+06, sampleId=tstrsplit(fn, "/")[[8]])
}
mpr <- rbindlist(mpr)

seqStats <- merge(seqStats, mpr, by="sampleId")

### --- Effective Number of Chromosomes NE ---

seqStats["effRD"] <- (seqStats$AveReadDepth * nbr.gene.copies) / (seqStats$AveReadDepth + nbr.gene.copies)


### Fixed sample ids
sample_ids <- c("MN-Av.1.140711", "MN-Av.1.140722", "MN-Av.1.140725", "MN-Av.1.140730", "MN-Av.1.AddRun", 
                "MN-V.1.140711",  "MN-V.1.140722",  "MN-V.1.140725",  "MN-V.1.140730", 
                "ND-Av.1.140711", "ND-Av.1.140722", "ND-Av.1.140725", "ND-Av.1.140730", 
                "ND-V.1.140711",  "ND-V.1.140722",  "ND-V.1.140725",  "ND-V.1.140730", 
                "NW-Av.1.140711", "NW-Av.1.140722", "NW-Av.1.140725", "NW-Av.1.140730", 
                "NW-Av.2.140711", "NW-Av.2.140722", "NW-Av.2.140725", "NW-Av.2.140730", 
                "NW-V.1.140711",  "NW-V.1.140722",  "NW-V.1.140725",  "NW-V.1.140730",  "NW-V.1.AddRun", 
                "NW-V.2.140711",  "NW-V.2.140722",  "NW-V.2.140725",  "NW-V.2.140730", 
                "NW-V.3.140711",  "NW-V.3.140722",  "NW-V.3.140725",  "NW-V.3.140730", 
                "PA-Av.1.140711", "PA-Av.1.140722", "PA-Av.1.140725", "PA-Av.1.140730", 
                "PA-V.1.140711",  "PA-V.1.140722",  "PA-V.1.140725",  "PA-V.1.140730", 
                "PA-V.2.140711",  "PA-V.2.140722",  "PA-V.2.140725",  "PA-V.2.140730",
                "WI-Av.1.140711", "WI-Av.1.140722", "WI-Av.1.140725", "WI-Av.1.140730", 
                "WI-Av.2.140711", "WI-Av.2.140722", "WI-Av.2.140725", "WI-Av.2.140730",  "WI-Av.2.AddRun", 
                "WI-V.1.140711",  "WI-V.1.140722",  "WI-V.1.140725",  "WI-V.1.140730", 
                "WI-V.2.140711",  "WI-V.2.140722",  "WI-V.2.140725",  "WI-V.2.140730", 
                "WI-V.3.140711",  "WI-V.3.140722",  "WI-V.3.140725",  "WI-V.3.140730", 
                "WO-Av.1.140711", "WO-Av.1.140722", "WO-Av.1.140725", "WO-Av.1.140730", 
                "WO-V.1.140711",  "WO-V.1.140722",  "WO-V.1.140725",  "WO-V.1.140730", 
                "WO-V.2.140711",  "WO-V.2.140722",  "WO-V.2.140725",  "WO-V.2.140730", 
                "WO-V.3.140711",  "WO-V.3.140722",  "WO-V.3.140725",  "WO-V.3.140730")

seqStats["fixSampleId"] <- sample_ids

# Save to a file
write.table(seqStats,
            file="results/technical_replicates/minmaxcov_3_99/sequencing_summaryStats.txt", sep = "\t",
            quote = F, row.names = F, col.names = TRUE)

###
###
### --- AGGREGATED DATA ---
###
###

poolsizes.aggre <- rep(10,21)

poolnames.aggre <- c("MN_BIO1_S1", "MN_BIO4_S1",
                     "ND_BIO1_S1", "ND_BIO4_S1", 
                     "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3", 
                     "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2", 
                     "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3", 
                     "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3")

dt.aggre <- vcf2pooldata(
  vcf.file = "vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz",
  poolsizes = poolsizes.aggre,
  poolnames = poolnames.aggre,
  min.cov.per.pool = -1,
  min.rc = 1,
  max.cov.per.pool = 1e+06,
  min.maf = 0.05,
  remove.indels = FALSE,
  nlines.per.readblock = 1e+06
)

dt.aggre
# * * * PoolData Object * * * 
# * Number of SNPs   =  262866 
# * Number of Pools  =  21

nbr.gene.copies = 10
nbr.loci = dt.aggre@nsnp

### --- Average Read Coverage and % of missing data ---

dp <- dt.aggre@readcoverage

dp[dp==0] <- NA

seqStats <- data.frame(sampleId = dt.aggre@poolnames, 
                       AveReadDepth = colMeans(dp, na.rm = T), 
                       propMissing=colSums(is.na(dp))/nbr.loci)

### --- Proportion of PCR duplicates and Number of Mapped reads ---

### PCR Duplicates
fns <- system("ls /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/aggregated/*_B*_S*/*duplicates_report.txt", intern=T)

pcr <- foreach(fn=fns)%do%{
  #fn <- fns[1]
  data.table(pcrDup=fread(fn, skip="LIBRARY", nrows=1)$PERCENT_DUPLICATION, sampleId=tstrsplit(fn, "/")[[9]])
}
pcr <- rbindlist(pcr)

seqStats <- merge(seqStats, pcr, by="sampleId")

### Mapped Reads
fns <- system("ls /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/aggregated/*_B*_S*/*.flagstat.txt", intern=T)

mpr <- foreach(fn=fns)%do%{
  #fn <- fns[1]
  data.table(mappedReads=fread(fn, nrow=1)$V1/1e+06, sampleId=tstrsplit(fn, "/")[[9]])
}
mpr <- rbindlist(mpr)

seqStats <- merge(seqStats, mpr, by="sampleId")

### --- Effective Number of Chromosomes NE ---

seqStats["effRD"] <- (seqStats$AveReadDepth * nbr.gene.copies) / (seqStats$AveReadDepth + nbr.gene.copies)


### Fixed sample ids
sample_ids <- c("MN-Av.1", "MN-V.1",
                "ND-Av.1", "ND-V.1", 
                "NW-Av.1", "NW-Av.2", "NW-V.1", "NW-V.2", "NW-V.3", 
                "PA-Av.1", "PA-V.1", "PA-V.2", 
                "WI-Av.1", "WI-Av.2", "WI-V.1", "WI-V.2", "WI-V.3", 
                "WO-Av.1", "WO-V.1", "WO-V.2", "WO-V.3")

seqStats["fixSampleId"] <- sample_ids

# Save to a file
write.table(seqStats,
            file="results/aggregated_data/minmaxcov_4_99/sequencing_summaryStats.txt", sep = "\t",
            quote = F, row.names = F, col.names = TRUE)
