#############################################################
#                   Exploratory analysis                    #
#           MAC, MAF, pN/pS and % missing data              #
#############################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()
## Remove last features
rm(list=ls())
ls()

library(poolfstat)
library(tidyverse)
library(gridExtra)
library(scales)
library(ggpubr)
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

## MIN_COV = 4; MAX_COV=0.99; SAMPLE OF SNPS

## Set global Alpha value
ALPHA=0.75

MAC = c(5, 10, 15, 20, 50, 100)
MAF = c(0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2)

###
# MISSING DATA AS A FUNCTION OF MAF & MAC
###

## REMEMBER THIS VCFs CONTAINS ONLY A SAMPLE OF SNPs

poolsizes <- rep(10,21)

poolnames <- c("MN_BIO1_S1", "MN_BIO4_S1",
               "ND_BIO1_S1", "ND_BIO4_S1", 
               "NW_BIO1_S1", "NW_BIO1_S2", "NW_BIO4_S1", "NW_BIO4_S2", "NW_BIO4_S3", 
               "PA_BIO1_S1", "PA_BIO4_S1", "PA_BIO4_S2", 
               "WI_BIO1_S1", "WI_BIO1_S2", "WI_BIO4_S1", "WI_BIO4_S2", "WI_BIO4_S3", 
               "WO_BIO1_S1", "WO_BIO4_S1", "WO_BIO4_S2", "WO_BIO4_S3")

#paste0("paramTest/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.", 0, ".", MAC[i], ".paramTest.vcf.gz")
##paste0("paramTest/aggregated_data/minmaxcov_4_95/aphidpool.PoolSeq.PoolSNP.", 0, ".", MAC[i], ".paramTest.vcf.gz")
##paste0("paramTest/aggregated_data/minmaxcov_5_95/snps_sample/aphidpool.PoolSeq.PoolSNP.", 0, ".", MAC[i], ".paramTest.vcf.gz")
##paste0("paramTest/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.", 0, ".", MAC[i], ".paramTest.vcf.gz")

missing_data=NULL
#t=NULL
for (i in seq_along(MAC))
{
  for (j in seq_along(MAF))
  {
   print(paste0("MAC:", MAC[i], " MAF:", MAF[j]))
   dt <- vcf2pooldata(vcf.file = paste0("paramTest/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.",0, ".", MAC[i],".paramTest.vcf"),
                      poolsizes = poolsizes,
                      poolnames = poolnames,
                      min.cov.per.pool = -1,
                      min.rc = 1,
                      max.cov.per.pool = 1e+20,
                      min.maf = MAF[j],
                      remove.indels = FALSE,
                      nlines.per.readblock = 1e+20)
   
   dt@readcoverage[dt@readcoverage==0] <- NA
   missing.pop <- 100*(colSums(is.na(dt@readcoverage))/dim(dt@readcoverage)[1])
   missing.overall <- 100*(sum(is.na(dt@readcoverage))/(dim(dt@readcoverage)[1] * dim(dt@readcoverage)[2]))
   res <- data.frame(MAC=MAC[i], MAF=MAF[j],  POP=c(poolnames, "Overall"), Missing=c(missing.pop, missing.overall))
   
   missing_data <- rbind(missing_data, res)
   #t_ <- paste0("MAC:", MAC[i], " MAF:", MAF[j])
   #t <- c(t, t_)
  }
}

write.table(missing_data, file = "paramTest/aggregated_data/minmaxcov_4_99/maf_mac_missing.mincov4.maxcov99.txt",
            row.names = F, col.names = T, quote = F)


missing_data <- read.table("paramTest/aggregated_data/minmaxcov_4_99/maf_mac_missing.mincov4.maxcov99.txt", header=TRUE, na.strings = NA)

## FOR EACH POOL
# Remove Overall data
DATA.missing.pools <- missing_data %>%
  filter(POP !="Overall") 

# Plot all pools
pools.missing.mac <- ggplot(DATA.missing.pools, aes(x=MAC, y=Missing, color=factor(MAF))) + 
                     geom_point() +
                     geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                     labs(color="MAF") + 
                     ylim(0,75) +
                     xlab("Minor allele count (MAC)") +
                     ylab("% of missing data") +
                     facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                     theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                        legend.position = "right", axis.text = element_text(size = 14),
                                        axis.title=element_text(size=16))
                   
pools.missing.maf <- ggplot(DATA.missing.pools, aes(x=MAF, y=Missing, color=factor(MAC))) + 
                     geom_point() +
                     geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                     labs(color="MAC") + 
                     ylim(0,75) +
                     xlab("Minor allele frequency (MAF)") +
                     ylab("% of missing data") +
                     facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                     theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                        legend.position = "right", axis.text = element_text(size = 14),
                                        axis.title=element_text(size=16))



table(DATA.missing.pools[which(DATA.missing.pools$Missing >= 30), 3])

# Mincov=3; Maxcov=99
#NW_BIO1_S1 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 WO_BIO4_S3 
#65         65         48         52         52 

# Mincov=4; Maxcov=95
#NW_BIO1_S1 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 WO_BIO1_S1 WO_BIO4_S3 
#78         78         52         52         17         65

# Mincov=4; Maxcov=99
#NW_BIO1_S1 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 WO_BIO4_S3 
#65         65         52         52         65

# Mincov=4; Maxcov=99 (more mafs)
#NW_BIO1_S1 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 WO_BIO4_S3 
#75         75         60         60         75 

# Mincov=5; Maxcov=95
#NW_BIO1_S1 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 WO_BIO1_S1 WO_BIO4_S3 
#78         78         60         52         36         66 

# Filter pools with > 30% missing data
DATA.missing.pools.filtered <- DATA.missing.pools %>%
  filter(POP !="NW_BIO1_S1" & POP !="PA_BIO1_S1" & POP !="PA_BIO4_S1" & POP !="PA_BIO4_S2" & POP !="WO_BIO4_S3")

# Plot only kept pools
pools.missing.mac.filtered <- ggplot(DATA.missing.pools.filtered, aes(x=MAC, y=Missing, color=factor(MAF))) + 
                              geom_point() +
                              geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                              labs(color="MAF") + 
                              ylim(0,75) +
                              xlab("Minor allele count (MAC)") +
                              ylab("% of missing data") +
                              facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                              theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                 legend.position = "right", axis.text = element_text(size = 14),
                                                 axis.title=element_text(size=16))
                             
pools.missing.maf.filtered <- ggplot(DATA.missing.pools.filtered, aes(x=MAF, y=Missing, color=factor(MAC))) + 
                              geom_point() +
                              geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                              labs(color="MAC") + 
                              ylim(0,75) +
                              xlab("Minor allele frequency (MAF)") +
                              ylab("% of missing data") +
                              facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                              theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                 legend.position = "right", axis.text = element_text(size = 14),
                                                 axis.title=element_text(size=16))
                             

par(mfrow=c(3, 1))
missingP<-ggarrange(ggarrange(pools.missing.mac,
                              pools.missing.mac.filtered,
                              labels = "A",
                              common.legend = T,
                              legend = "bottom",
                              font.label = list(size = 24, face = "bold"),
                              nrow=2),
                    ggarrange(pools.missing.maf,
                              pools.missing.maf.filtered,
                              labels = "B",
                              common.legend = T,
                              legend = "bottom",
                              font.label = list(size = 24, face = "bold"),
                              nrow=2))

missingP

ggsave("results/aggregated_data/minmaxcov_4_99/Figure_missing_maf_mac_pools.pdf",
       missingP,
       device="pdf",
       width=15,
       height=9)


## FOR THE COMBINED POOLS
# Keep only Overall data
DATA.missing.overall <- missing_data %>%
  filter(POP =="Overall")

# Plot with all pools
overall.missing.mac <- ggplot(DATA.missing.overall, aes(x=MAC, y=Missing, color=factor(MAF))) + 
                              geom_point() +
                              geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                              labs(color="MAF") + 
                              ylim(0,75) +
                              xlab("Minor allele count (MAC)") +
                              ylab("% of missing data") +
                              facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                              theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                 legend.position = "right", axis.text = element_text(size = 14),
                                                 axis.title=element_text(size=16))

overall.missing.maf <- ggplot(DATA.missing.overall, aes(x=MAF, y=Missing, color=factor(MAC))) + 
                       geom_point() +
                       geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                       labs(color="MAC") + 
                       ylim(0,75) +
                       xlab("Minor allele frequency (MAF)") +
                       ylab("% of missing data") +
                       facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                       theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                          legend.position = "right", axis.text = element_text(size = 14),
                                          axis.title=element_text(size=16))
                                                
## FILTER OUT POOL > 30% DATA AND CALCULATE OVERALL
missing_data_filtered=NULL
for (i in seq_along(MAC))
{
  for (j in seq_along(MAF))
  {
   print(paste0("MAC:", MAC[i], " MAF:", MAF[j]))
   dt <- vcf2pooldata(vcf.file = paste0("paramTest/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.",0, ".", MAC[i],".paramTest.vcf"),
                      poolsizes = poolsizes,
                      poolnames = poolnames,
                      min.cov.per.pool = -1,
                      min.rc = 1,
                      max.cov.per.pool = 1e+20,
                      min.maf = MAF[j],
                      remove.indels = FALSE,
                      nlines.per.readblock = 1e+20)
   
   dt@readcoverage[dt@readcoverage==0] <- NA
   dt.rc <- dt@readcoverage
   dt.rc.filtered <- dt.rc[,-c(5,10,11,12,21)]
   
   missing.pop <- 100*(colSums(is.na(dt.rc.filtered))/dim(dt.rc.filtered)[1])
   missing.overall <- 100*(sum(is.na(dt.rc.filtered))/(dim(dt.rc.filtered)[1] * dim(dt.rc.filtered)[2]))
   res <- data.frame(MAC=MAC[i], MAF=MAF[j],  POP=c(poolnames[-c(5,10,11,12,21)], "Overall"), 
                     Missing=c(missing.pop, missing.overall))
   
   missing_data_filtered <- rbind(missing_data_filtered, res)
  }
}

write.table(missing_data_filtered, file = "paramTest/aggregated_data/minmaxcov_4_99/maf_mac_missing_filtered.mincov4.maxcov99.txt",
            row.names = F, col.names = T, quote = F)


missing_data_filtered <- read.table("paramTest/aggregated_data/minmaxcov_4_99/maf_mac_missing_filtered.mincov4.maxcov99.txt", header=TRUE, na.strings = NA)

# Keep only Overall data after filtering pools
DATA.missing.overall.filtered <- missing_data_filtered %>%
  filter(POP =="Overall")

# Plot with all pools
overall.missing.mac.filtered <- ggplot(DATA.missing.overall.filtered, aes(x=MAC, y=Missing, color=factor(MAF))) + 
                                geom_point() +
                                geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                                labs(color="MAF") + 
                                ylim(0,75) +
                                xlab("Minor allele count (MAC)") +
                                ylab("% of missing data") +
                                facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                                theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                   legend.position = "right", axis.text = element_text(size = 14),
                                                   axis.title=element_text(size=16))
                      
overall.missing.maf.filtered <- ggplot(DATA.missing.overall.filtered, aes(x=MAF, y=Missing, color=factor(MAC))) + 
                                geom_point() +
                                geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                                labs(color="MAC") + 
                                ylim(0,75) +
                                xlab("Minor allele frequency (MAF)") +
                                ylab("% of missing data") +
                                facet_wrap(vars("Pool's % of missing data ('min cov' = 4; 'max cov' = 99%)")) +
                                theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                   legend.position = "right", axis.text = element_text(size = 14),
                                                   axis.title=element_text(size=16))

par(mfrow=c(3, 1))
missingP.overall <- ggarrange(ggarrange(overall.missing.mac,
                                        overall.missing.mac.filtered,
                                        labels = "A",
                                        common.legend = T,
                                        legend = "bottom",
                                        font.label = list(size = 24, face = "bold"),
                                        nrow=2),
                              ggarrange(overall.missing.maf,
                                        overall.missing.maf.filtered,
                                        labels = "B",
                                        common.legend = T,
                                        legend = "bottom",
                                        font.label = list(size = 24, face = "bold"),
                                        nrow=2))

missingP.overall

ggsave("results/aggregated_data/minmaxcov_4_99/Figure_missing_maf_mac_global.pdf",
       missingP.overall,
       device="pdf",
       width=15,
       height=9)


###
# FST, SAMPLE SIZE AND REALIZED MAF
###

poolfstat_data=NULL
for (i in seq_along(MAC))
{
  for (j in seq_along(MAF))
  {
    print(paste0("MAC:", MAC[i], " MAF:", MAF[j]))
    dt <- vcf2pooldata(vcf.file = paste0("paramTest/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.",0, ".", MAC[i],".paramTest.vcf"),
                       poolsizes = poolsizes,
                       poolnames = poolnames,
                       min.cov.per.pool = -1,
                       min.rc = 1,
                       max.cov.per.pool = 1e+20,
                       min.maf = MAF[j],
                       remove.indels = FALSE,
                       nlines.per.readblock = 1e+20)
    
    
    # Realized MAF calculations
    imputedMLCount <- imputedRefMLCount(dt)
    imputedRefMLFreq <- imputedMLCount[[1]]/imputedMLCount[[3]]
    
    a <- apply(imputedRefMLFreq, 1,  function(x) mean(x, na.rm = T))
    b <- apply(imputedRefMLFreq, 1,  function(x) 1-mean(x, na.rm = T))
    
    p <- pmin(a, b)
    maf <- min(p)
    
    # Fst
    dt.fst <- computeFST(dt, method = "Anova", nsnp.per.bjack.block = 10, verbose=FALSE)
    
    gfst = dt.fst$FST 
    mean.gfst = dt.fst$mean.fst
    se.gfst = dt.fst$mean.fst+c(-1.96,1.96)*dt.fst$se.fst
    
    
    res <- data.frame(MAC=MAC[i], MAF=MAF[j], NSNPs= dt@nsnp, rMAF=maf, 
                      gFST=gfst, mean_gFST=mean.gfst, l_se_gFST=se.gfst[1], u_se_gFST=se.gfst[2])
    
    poolfstat_data <- rbind(poolfstat_data, res)
  }
}

write.table(poolfstat_data, file = "paramTest/aggregated_data/minmaxcov_4_99/nsnps_rmaf_fst.mincov4.maxcov99.txt",
            row.names = F, col.names = T, quote = F)


poolfstat_data <- read.table("paramTest/aggregated_data/minmaxcov_4_99/nsnps_rmaf_fst.mincov4.maxcov99.txt", header=TRUE, na.strings = NA)

# Number of SNPs
number.snps.mac <- ggplot(poolfstat_data, aes(x=MAC, y=NSNPs, color=factor(MAF))) + 
                          geom_point() +
                          geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                          labs(color="MAF") + 
                          #ylim(0,75) +
                          xlab("Minor allele count (MAC)") +
                          ylab("Number of SNPs") +
                          facet_wrap(vars("Pool's number of SNPs ('min cov' = 4; 'max cov' = 99%)")) +
                          theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                             legend.position = "right", axis.text = element_text(size = 14),
                                             axis.title=element_text(size=16))
                      
number.snps.maf <- ggplot(poolfstat_data, aes(x=MAF, y=NSNPs, color=factor(MAC))) + 
                   geom_point() +
                   geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                   labs(color="MAC") + 
                   #ylim(0,75) +
                   xlab("Minor allele frequency (MAF)") +
                   ylab("Number of SNPs") +
                   facet_wrap(vars("Pool's number of SNPs ('min cov' = 4; 'max cov' = 99%)")) +
                   theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                      legend.position = "right", axis.text = element_text(size = 14),
                                      axis.title=element_text(size=16))

## Realized MAF
#rMAF.mac <- ggplot(poolfstat_data, aes(x=MAC, y=rMAF, color=factor(MAF))) + 
#            geom_point() +
#            geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
#            labs(color="MAF") + 
#            #ylim(0,75) +
#            xlab("Minor allele count (MAC)") +
#            ylab("Realized MAF") +
#            facet_wrap(vars("Pool's realized MAF ('min cov' = 4; 'max cov' = 99%)")) +
#            theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                               panel.grid.major = element_blank(),
#                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                               legend.position = "right", axis.text = element_text(size = 14),
#                               axis.title=element_text(size=16))
#          
#rMAF.maf <- ggplot(poolfstat_data, aes(x=MAF, y=rMAF, color=factor(MAC))) + 
#                   geom_point() +
#                   geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
#                   labs(color="MAC") + 
#                   #ylim(0,75) +
#                   xlab("Minor allele frequency (MAF)") +
#                   ylab("Realized MAF") +
#                   facet_wrap(vars("Pool's realized MAF ('min cov' = 4; 'max cov' = 99%)")) +
#                   theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                                      panel.grid.major = element_blank(),
#                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                                      legend.position = "right", axis.text = element_text(size = 14),
#                                      axis.title=element_text(size=16))
#
poolfstat_data[poolfstat_data$MAF == 0.001, ]
#MAC   MAF NSNPs        rMAF  
#2    5 0.001 52350 0.004761905
#17  10 0.001 45098 0.014285714
#32  15 0.001 41126 0.022222222
#47  20 0.001 37473 0.022222222
#62  50 0.001 17847 0.055555556
#77 100 0.001  5709 0.077777778

# global FST
gfst.mac <- ggplot(poolfstat_data, aes(x=MAC, y=gFST, color=factor(MAF))) + 
            geom_point() +
            geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
            labs(color="MAF") + 
            #ylim(0,75) +
            xlab("Minor allele count (MAC)") +
            ylab("gFST") +
            facet_wrap(vars("Pool's gFST ('min cov' = 4; 'max cov' = 99%)")) +
            theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                               legend.position = "right", axis.text = element_text(size = 14),
                               axis.title=element_text(size=16))
          
gfst.maf <- ggplot(poolfstat_data, aes(x=MAF, y=gFST, color=factor(MAC))) + 
                   geom_point() +
                   geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                   labs(color="MAC") + 
                   #ylim(0,75) +
                   xlab("Minor allele frequency (MAF)") +
                   ylab("gFST") +
                   facet_wrap(vars("Pool's gFST ('min cov' = 4; 'max cov' = 99%)")) +
                   theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                      legend.position = "right", axis.text = element_text(size = 14),
                                      axis.title=element_text(size=16))
                                                    
par(mfrow=c(3, 1))
nsnps.fst.rmaf <- ggarrange(ggarrange(number.snps.mac,
                                      gfst.mac,
                                      labels = "A",
                                      common.legend = T,
                                      legend = "bottom",
                                      font.label = list(size = 24, face = "bold"),
                                      nrow=2),
                              ggarrange(number.snps.maf,
                                        gfst.maf,
                                        labels = "B",
                                        common.legend = T,
                                        legend = "bottom",
                                        font.label = list(size = 24, face = "bold"),
                                        nrow=2))

nsnps.fst.rmaf

ggsave("results/aggregated_data/minmaxcov_4_99/Figure_nsnps_gfst_maf_mac.pdf",
       nsnps.fst.rmaf,
       device="pdf",
       width=15,
       height=9)

###
# pN/pS DATA AS A FUNCTION OF MAF & MAC
###

## Six MAC many MAFs

# Combine multiples files
macs.mafs = NULL
for (i in seq_along(MAC)){
  mac.file <- read.table(file = paste0("paramTest/aggregated_data/minmaxcov_4_99/PoolSNP.pnps.mafs.", MAC[i],".mincov4.maxcov99.gz"), 
                         header = T, na.strings = NA)
  mac.file$MAC <- MAC[i]
  
  macs.mafs <- rbind(macs.mafs,mac.file)
  
}

write.table(macs.mafs, file = "paramTest/aggregated_data/minmaxcov_4_99/maf_mac_pnps.mincov4.maxcov99.txt",
            row.names = F, col.names = T, quote = F)

macs.mafs <- read.table("paramTest/aggregated_data/minmaxcov_4_99/maf_mac_pnps.mincov4.maxcov99.txt", header=TRUE, na.strings = NA)


# Keep only genomewid data
DATA.pnps.pools.mafs <- macs.mafs %>%
  filter(Chrom =="genomewide")

# Plot 
genomewide.pspn.mac <- ggplot(DATA.pnps.pools.mafs, aes(x=MAC, y=pNpS, color=factor(MAF))) + 
                       geom_point() +
                       geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                       labs(color="MAF") + 
                       ylim(0.6,2.2) +
                       xlab("Minor allele count (MAC)") +
                       ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                       facet_wrap(vars("Pool's pNpS ('min cov' = 4; 'max cov' = 99%)")) +
                       theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                          legend.position = "right", axis.text = element_text(size = 14),
                                          axis.title=element_text(size=16))
                      
genomewide.pspn.maf <- ggplot(DATA.pnps.pools.mafs, aes(x=MAF, y=pNpS, color=factor(MAC))) + 
                       geom_point() +
                       geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                       labs(color="MAC") + 
                       ylim(0.6,2.2) +
                       xlab("Minor allele frequency (MAF)") +
                       ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                       facet_wrap(vars("Pool's pNpS ('min cov' = 4; 'max cov' = 99%)")) +
                       theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                          legend.position = "right", axis.text = element_text(size = 14),
                                          axis.title=element_text(size=16))


# Keep only genomewide data for pools that passed missing data filtering
DATA.pnps.pools.mafs.filtered <- DATA.pnps.pools.mafs %>%
  filter(POP !="NW_BIO1_S1" & POP !="PA_BIO1_S1" & POP !="PA_BIO4_S1" & POP !="PA_BIO4_S2" & POP !="WO_BIO4_S3")
  

genomewide.pspn.mac.filtered <- ggplot(DATA.pnps.pools.mafs.filtered, aes(x=MAC, y=pNpS, color=factor(MAF))) + 
                                geom_point() +
                                geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                                labs(color="MAF") + 
                                ylim(0.6,2.2) +
                                xlab("Minor allele count (MAC)") +
                                ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                                facet_wrap(vars("Pool's pNpS ('min cov' = 4; 'max cov' = 99%)")) +
                                theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                   legend.position = "right", axis.text = element_text(size = 14),
                                                   axis.title=element_text(size=16))
                               
genomewide.pspn.maf.filtered <- ggplot(DATA.pnps.pools.mafs.filtered, aes(x=MAF, y=pNpS, color=factor(MAC))) + 
                                geom_point() +
                                geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
                                labs(color="MAC") + 
                                ylim(0.6,2.2) +
                                xlab("Minor allele frequency (MAF)") +
                                ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                                facet_wrap(vars("Pool's pNpS ('min cov' = 4; 'max cov' = 99%)")) +
                                theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                   legend.position = "right", axis.text = element_text(size = 14),
                                                   axis.title=element_text(size=16))
                               


par(mfrow=c(3, 1))
pnpsP <- ggarrange(ggarrange(genomewide.pspn.mac,
                             genomewide.pspn.mac.filtered,
                             labels = "A",
                             common.legend = T,
                             legend = "bottom",
                             font.label = list(size = 24, face = "bold"),
                             nrow=2),
                   ggarrange(genomewide.pspn.maf,
                             genomewide.pspn.maf.filtered,
                             labels = "B",
                             common.legend = T,
                             legend = "bottom",
                             font.label = list(size = 24, face = "bold"),
                             nrow=2))

pnpsP

ggsave("results/aggregated_data/minmaxcov_4_99/Figure_pnps_maf_mac.pdf",
       pnpsP,
       device="pdf",
       width=15,
       height=9)


# Calculate the mean, sd and se for each MAF and MAC
DATA.pnps.mean.maf.mac <- DATA.pnps.pools.mafs %>%
  group_by(MAF,MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS, na.rm = T),
    pnps.sd=sd(pNpS, na.rm = T),
    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
  )

DATA.pnps.mean.maf.mac.selected.maf <- DATA.pnps.mean.maf.mac %>%
  filter(MAF == 0.005 | MAF == 0.01 | MAF == 0.05 | MAF == 0.1)

DATA.pnps.mean.maf.mac.selected.mac <- DATA.pnps.mean.maf.mac %>%
  filter(MAC == 5 | MAC == 20 | MAC == 50 | MAC == 100)

# More fancy plot

# All MAC as a function of MAF
genomewide.mean.pspn.mac <- ggplot(DATA.pnps.mean.maf.mac.selected.maf, aes(x=MAC, y=pnps.m, color=factor(MAF))) + 
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.005), aes(x=MAC,y=pNpS)) +
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.01), aes(x=MAC,y=pNpS)) +
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.05), aes(x=MAC,y=pNpS)) +
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.1), aes(x=MAC,y=pNpS)) +
                            geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAF)), colour = NA, alpha=0.2) +
                            geom_line(lwd=1.5) + 
                            xlab("Minor allele count (MAC)") +
                            ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                            facet_wrap(vars("PoolSeq all pools")) +
                            theme_bw() +
                            theme(axis.title.y = element_text(size = 20, angle = 90),
                                  axis.title.x = element_text(size = 20, angle = 0),
                                  axis.text=element_text(size=18),
                                  strip.text =element_text(size=20),
                                  legend.position = "right")+
                            scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
                          scale_fill_manual(values=c(alpha(c("#999999"),
                                                           alpha=ALPHA),
                                                     alpha(c("#de77ae"),
                                                           alpha=ALPHA),
                                                     alpha(c("#abd9e9"),
                                                           alpha=ALPHA),
                                                     alpha(c("#fee08b"),
                                                           alpha=ALPHA)), name="MAF")+
                          scale_colour_manual(values=c(alpha(c("#999999"),
                                                             alpha=ALPHA),
                                                        alpha(c("#de77ae"),
                                                             alpha=ALPHA),
                                                       alpha(c("#abd9e9"),
                                                             alpha=ALPHA),
                                                       alpha(c("#fee08b"),
                                                             alpha=ALPHA)), name="MAF")
genomewide.mean.pspn.mac





# All MAF as a function of MAC
genomewide.mean.pspn.maf <- ggplot(DATA.pnps.mean.maf.mac.selected.mac, aes(x=MAF, y=pnps.m, color=factor(MAC))) + 
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==5), aes(x=MAF,y=pNpS)) +
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==20), aes(x=MAF,y=pNpS)) +
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==50), aes(x=MAF,y=pNpS)) +
                            geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==100), aes(x=MAF,y=pNpS)) +
                            geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAC)), colour = NA, alpha=0.2) +
                            geom_line(lwd=1.5) + 
                            xlab("Minor allele frequency (MAF)") +
                            ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                            facet_wrap(vars("PoolSeq all pools")) +
                            theme_bw() +
                            theme(axis.title.y = element_text(size = 20, angle = 90),
                                  axis.title.x = element_text(size = 20, angle = 0),
                                  axis.text=element_text(size=18),
                                  strip.text =element_text(size=20),
                                  legend.position = "right") +
                            scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
                            scale_fill_manual(values=c(alpha(c("#999999"),
                                                             alpha=ALPHA),
                                                       alpha(c("#E69F00"),
                                                             alpha=ALPHA),
                                                       alpha(c("#56B4E9"),
                                                             alpha=ALPHA),
                                                       alpha(c("#009E73"),
                                                             alpha=ALPHA)), name="MAC", breaks = c("5", "20", "50", "100"))+
                            scale_colour_manual(values=c(alpha(c("#999999"),
                                                               alpha=ALPHA),
                                                         alpha(c("#E69F00"),
                                                               alpha=ALPHA),
                                                         alpha(c("#56B4E9"),
                                                               alpha=ALPHA),
                                                         alpha(c("#009E73"),
                                                               alpha=ALPHA)), name="MAC", breaks = c("5", "20", "50", "100"))
genomewide.mean.pspn.maf


# Filtered data
# Calculate the mean, sd and se for each MAF and MAC
DATA.pnps.mean.maf.mac.filtered <- DATA.pnps.pools.mafs.filtered %>%
  group_by(MAF,MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS, na.rm = T),
    pnps.sd=sd(pNpS, na.rm = T),
    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
  )

DATA.pnps.mean.maf.mac.selected.maf.filtered <- DATA.pnps.mean.maf.mac.filtered %>%
  filter(MAF == 0.005 | MAF == 0.01 | MAF == 0.05 | MAF == 0.1)

DATA.pnps.mean.maf.mac.selected.mac.filtered <- DATA.pnps.mean.maf.mac.filtered %>%
  filter(MAC == 5 | MAC == 20 | MAC == 50 | MAC == 100)

# More fancy plot
# All MAC as a function of MAF
genomewide.mean.pspn.mac.filtered <- ggplot(DATA.pnps.mean.maf.mac.selected.maf.filtered, aes(x=MAC, y=pnps.m, color=factor(MAF))) + 
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.005), aes(x=MAC,y=pNpS)) +   
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.01), aes(x=MAC,y=pNpS)) +
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.05), aes(x=MAC,y=pNpS)) +
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.1), aes(x=MAC,y=pNpS)) +
                                     geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAF)), colour = NA, alpha=0.2) +
                                     geom_line(lwd=1.5) + 
                                     xlab("Minor allele count (MAC)") +
                                     ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                                     facet_wrap(vars("PoolSeq pools < 30% missing data")) +
                                     theme_bw() +
                                     theme(axis.title.y = element_text(size = 20, angle = 90),
                                           axis.title.x = element_text(size = 20, angle = 0),
                                           axis.text=element_text(size=18),
                                           strip.text =element_text(size=20),
                                           legend.position = "right")+
                                     scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
                                     scale_fill_manual(values=c(alpha(c("#999999"),
                                                                      alpha=ALPHA),
                                                                alpha(c("#de77ae"),
                                                                      alpha=ALPHA),
                                                                alpha(c("#abd9e9"),
                                                                      alpha=ALPHA),
                                                                alpha(c("#fee08b"),
                                                                      alpha=ALPHA)), name="MAF")+
                                     scale_colour_manual(values=c(alpha(c("#999999"),
                                                                        alpha=ALPHA),
                                                                  alpha(c("#de77ae"),
                                                                        alpha=ALPHA),
                                                                  alpha(c("#abd9e9"),
                                                                        alpha=ALPHA),
                                                                  alpha(c("#fee08b"),
                                                                        alpha=ALPHA)), name="MAF")
genomewide.mean.pspn.mac.filtered


# All MAF as a function of MAC
genomewide.mean.pspn.maf.filtered <- ggplot(DATA.pnps.mean.maf.mac.selected.mac.filtered, aes(x=MAF, y=pnps.m, color=factor(MAC))) + 
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==5), aes(x=MAF,y=pNpS)) +
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==20), aes(x=MAF,y=pNpS)) +
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==50), aes(x=MAF,y=pNpS)) +
                                     geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==100), aes(x=MAF,y=pNpS)) +
                                     geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAC)), colour = NA, alpha=0.2) +
                                     geom_line(lwd=1.5) + 
                                     xlab("Minor allele frequency (MAF)") +
                                     ylab(expression(italic("p")["N"]/italic("p")["S"])) +
                                     facet_wrap(vars("PoolSeq pools < 30% missing data")) +
                                     theme_bw() +
                                     theme(axis.title.y = element_text(size = 20, angle = 90),
                                           axis.title.x = element_text(size = 20, angle = 0),
                                           axis.text=element_text(size=18),
                                           strip.text =element_text(size=20),
                                           legend.position = "right") +
                                     scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
                                     scale_fill_manual(values=c(alpha(c("#999999"),
                                                                      alpha=ALPHA),
                                                                alpha(c("#E69F00"),
                                                                      alpha=ALPHA),
                                                                alpha(c("#56B4E9"),
                                                                      alpha=ALPHA),
                                                                alpha(c("#009E73"),
                                                                      alpha=ALPHA)), name="MAC", breaks = c("5", "20", "50", "100"))+
                                     scale_colour_manual(values=c(alpha(c("#999999"),
                                                                        alpha=ALPHA),
                                                                  alpha(c("#E69F00"),
                                                                        alpha=ALPHA),
                                                                  alpha(c("#56B4E9"),
                                                                        alpha=ALPHA),
                                                                  alpha(c("#009E73"),
                                                                        alpha=ALPHA)), name="MAC", breaks = c("5", "20", "50", "100"))
genomewide.mean.pspn.maf.filtered


par(mfrow=c(3, 1))
mean.pnpsP <- ggarrange(ggarrange(genomewide.mean.pspn.mac,
                             genomewide.mean.pspn.mac.filtered,
                             labels = "A",
                             common.legend = T,
                             legend = "bottom",
                             font.label = list(size = 24, face = "bold"),
                             nrow=2),
                   ggarrange(genomewide.mean.pspn.maf,
                             genomewide.mean.pspn.maf.filtered,
                             labels = "B",
                             common.legend = T,
                             legend = "bottom",
                             font.label = list(size = 24, face = "bold"),
                             nrow=2))

mean.pnpsP


ggsave("results/aggregated_data/minmaxcov_4_99/Figure_fancy_pnps_maf_mac.pdf",
       mean.pnpsP,
       device="pdf",
       width=15,
       height=9)




###
# Identify samples to remove based on % Missing, MAC, MAF and pNpS
###

DATA.pnps.genomewide.maf0 <- macs.mafs %>%
  filter(Chrom =="genomewide" & MAF==0) %>%
  group_by(POP) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )

DATA.missing.group.maf0 <- DATA.missing.pools %>%
  filter(MAF==0) %>%
  group_by(POP) %>%
  summarise(
    n=n(),
    missing.m = mean(Missing),
    missing.sd=sd(Missing),
    missing.se =sd(Missing)/sqrt(sum(n()))
  )

DATA=merge(DATA.missing.group.maf0 ,DATA.pnps.genomewide.maf0, by.x="POP",by.y="POP")

## calculated Mean/SD and threshold based on Mean+2SD for pNpS data
Mean.pNpS=mean(DATA$pnps.m)
SD.pNpS=sd(DATA$pnps.m)
th.pNpS=Mean.pNpS-1.96*SD.pNpS

##classify
DATA$TH.pNpS <-DATA$pnps.m
DATA$TH.pNpS[DATA$TH.pNpS>th.pNpS]<-NA
DATA$TH.pNpS[DATA$TH.pNpS<=th.pNpS]<-"Exclude"
DATA$TH.pNpS[is.na(DATA$TH.pNpS)]<-"Keep"

## calculated Mean/SD and threshold based on Mean+1.96SD for %missing
Mean.missing=mean(DATA$missing.m)
SD.missing=sd(DATA$missing.m)
th.missing=Mean.missing+1.96*SD.missing

#classify
DATA$TH.missing <-DATA$missing.m
DATA$TH.missing[DATA$TH.missing < 30]<-NA
DATA$TH.missing[DATA$TH.missing >= 30]<-"Exclude"
DATA$TH.missing[is.na(DATA$TH.missing)]<-"Keep"

# make new final classification, where any pop will be excluded that is excluded based on either measurement (pNpS and/or %missing)
DATA$Status <- rep("Keep",length(DATA$TH.missing))
DATA$Status[DATA$TH.pNpS=="Exclude"]<-"Exclude"
DATA$Status[DATA$TH.private=="Exclude"]<-"Exclude"

Classify.plot<-ggplot(DATA,aes(x = pnps.m ,y = missing.m, col = TH.missing))+
               geom_point(size=10)+
               xlab(expression(italic("p")["N"]/italic("p")["S"]))+
               ylab("% of Missing Data")+
               theme_bw()+
               geom_hline(yintercept= 30, lty=2,col="black",lwd=1)+
               geom_vline(xintercept=th.pNpS,lty=2,col="black",lwd=1)+
               theme(axis.title.y = element_text(size = 20, angle = 90),
                     axis.title.x = element_text(size = 20, angle = 0),
                     axis.text=element_text(size=18),
                     strip.text =element_text(size=20),
                     legend.position = "bottom",
                     legend.title = element_text(color = "black", size = 20),
                     legend.text=element_text(size=14))+
               scale_fill_manual(values=c(alpha(c("red"),
                                                alpha=ALPHA-0.50),
                                          alpha(c("black"),
                                                alpha=ALPHA-0.50)), name = "Pool status")+
               scale_colour_manual(values=c(alpha(c("red"),
                                                alpha=ALPHA-0.50),
                                          alpha(c("black"),
                                                alpha=ALPHA-0.50)), name = "Pool status") +
               facet_wrap(vars("PoolSNP MAF = 0"))
Classify.plot


LIST<-DATA$POP[DATA$TH.missing=="Exclude"]
STATUSLIST=data.frame(POP=DATA$POP,Status=DATA$TH.missing)

## keep pnps without MAF filtering only
DATA.poolsnp.group.MAF0 <- macs.mafs %>%
  filter(MAF == 0.05 & MAC==5 & Chrom =="genomewide")

macs.mafs<-merge(macs.mafs, STATUSLIST,by.x="POP",by.y="POP")

DATA.pnps.group.mac.plot <- macs.mafs %>%
  filter(Chrom =="genomewide" & MAF==0.05) %>%
  group_by(MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )


DATA.poolsnp.plot.mac <- ggplot(DATA.pnps.group.mac.plot, aes(x=MAC,y=pnps.m))+
                         geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.05),
                                     aes(x=MAC,y=pNpS,col=Status))+
                         geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
                                     colour = NA,
                                     fill=alpha(c("blue"),alpha=0.2))+
                         geom_line(lwd=1.5, col="blue")+
                         xlab("Minor allele count (MAC)")+
                         ylab(expression(italic("p")["N"]/italic("p")["S"]))+
                         theme_bw()+
                         theme(axis.title.y = element_text(size = 20, angle = 90),
                               axis.title.x = element_text(size = 20, angle = 0),
                               axis.text=element_text(size=18),
                               strip.text =element_text(size=20),
                               legend.title = element_blank(),
                               legend.position = "none")+
                         scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                                            limits=c(0,2.1))+
                         scale_fill_manual(values=c(alpha(c("red"),
                                                          alpha=ALPHA-0.50),
                                                    alpha(c("black"),
                                                          alpha=ALPHA-0.50)))+
                         scale_colour_manual(values=c(alpha(c("red"),
                                                          alpha=ALPHA-0.50),
                                                    alpha(c("black"),
                                                          alpha=ALPHA-0.50)))+
                         facet_wrap(vars("PoolSNP MAF = 0.05"))
DATA.poolsnp.plot.mac


DATA.pnps.group.maf.plot <- macs.mafs %>%
  filter(Chrom =="genomewide" & MAC==5) %>%
  group_by(MAF) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )


DATA.poolsnp.plot.maf <- ggplot(DATA.pnps.group.maf.plot, aes(x=MAF,y=pnps.m))+
                          geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==5),
                                      aes(x=MAF,y=pNpS,col=Status))+
                          geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd),
                                      colour = NA,
                                      fill=alpha(c("blue"),alpha=0.2))+
                          geom_line(lwd=1.5,
                                    col="blue")+
                          xlab("Minor allele Frequency (MAF)")+
                          ylab(expression(italic("p")["N"]/italic("p")["S"]))+
                          theme_bw()+
                          theme(axis.title.y = element_text(size = 20, angle = 90),
                                axis.title.x = element_text(size = 20, angle = 0),
                                axis.text=element_text(size=18),
                                strip.text =element_text(size=20),
                                legend.title = element_blank(),
                                legend.position = "none")+
                          scale_y_continuous(labels = scales::number_format(accuracy = 0.05),
                                             limits=c(0,2.1))+
                          scale_fill_manual(values=c(alpha(c("red"),
                                                           alpha=ALPHA-0.50),
                                                     alpha(c("black"),
                                                           alpha=ALPHA-0.50)))+
                          scale_colour_manual(values=c(alpha(c("red"),
                                                             alpha=ALPHA-0.50),
                                                       alpha(c("black"),
                                                             alpha=ALPHA-0.50)))+
                          facet_wrap(vars("PoolSNP MAC = 5"))
DATA.poolsnp.plot.maf


FP<-ggarrange(Classify.plot,
              ggarrange(DATA.poolsnp.plot.mac,
                        DATA.poolsnp.plot.maf,
                        labels = c("B","C"),
                        font.label = list(size = 28, face = "bold"),
                        nrow=2),
              ncol=2,
              labels="A",
              common.legend = T,
              legend = "bottom",
              font.label = list(size = 28, face = "bold"))
FP

ggsave("results/aggregated_data/minmaxcov_4_99/pnps_mac_maf/Figure_missing_pnps.pdf",
       FP,
       device="pdf",
       width=15,
       height=9)


###
## Paramtest - ALL SNPS - MINCOV=5; MAXCOV=95
###

###
# TS/TV AS A FUNCTION OF MAF & MAC
##

tstv_data <- read.csv("paramTest/aggregated_data/minmaxcov_5_95/snps_all/maf_mac_tstv.mincov5.maxcov95.csv", header=TRUE, na.strings = NA)
tstv_data <- tstv_data[,-c(5:15)]

## FOR EACH POOL
# Remove Combined data
DATA.tstv.pools <- tstv_data %>%
  filter(POP !="Combined") 

# Plot ts/ts ~ MAC
pools.mac <- ggplot(DATA.tstv.pools, aes(x=MAC, y=TS_TV, color=factor(MAF))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) + 
  labs(color="MAF") + 
  ylim(0.5,3) +
  xlab("Minor allele count (MAC)") +
  ylab(expression(italic("T")["S"]/italic("T")["V"])) +
  facet_wrap(vars("Pool's Ts/Tv ratios ('min cov' = 5; 'max cov' = 95%)")) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

# Plot ts/ts ~ MAF
pools.maf <- ggplot(DATA.tstv.pools, aes(x=MAF, y=TS_TV, color=factor(MAC))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
  labs(color="MAC") +
  ylim(0.8,2.7) +
  xlab("Minor allele frequency (MAF)") +
  ylab(expression(italic("T")["S"]/italic("T")["V"])) +
  facet_wrap(vars("Pool's Ts/Tv ratios ('min cov' = 5; 'max cov' = 95%)")) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

table(DATA.tstv.pools[which(DATA.tstv.pools$TS_TV > 1.8), 3])
# NW_BIO1_S1 PA_BIO1_S1 
# 15         18

## FOR THE COMBINED POOLS
# Only Combined data
DATA.tstv.combined <- tstv_data %>%
  filter(POP =="Combined") 

# Plot ts/ts ~ MAC
comb.mac <- ggplot(DATA.tstv.combined, aes(x=MAC, y=TS_TV, color=factor(MAF))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) + 
  labs(color="MAF") + 
  ylim(0.5,3) +
  xlab("Minor allele count (MAC)") +
  ylab(expression(italic("T")["S"]/italic("T")["V"])) +
  facet_wrap(vars("Combined Pools Ts/Tv ratios ('min cov' = 5; 'max cov' = 95%)")) + 
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

# Plot ts/ts ~ MAF
comb.maf <- ggplot(DATA.tstv.combined, aes(x=MAF, y=TS_TV, color=factor(MAC))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) + 
  labs(color="MAC") + 
  ylim(0.8,2.7) +
  xlab("Minor allele frequency (MAF)") +
  ylab(expression(italic("T")["S"]/italic("T")["V"])) +
  facet_wrap(vars("Combined Pools Ts/Tv ratios ('min cov' = 5; 'max cov' = 95%)")) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))


par(mfrow=c(3, 1))
tstsP<-ggarrange(ggarrange(pools.mac,
                           comb.mac,
                           labels = "A",
                           common.legend = T,
                           legend = "bottom",
                           font.label = list(size = 24, face = "bold"),
                           nrow=2),
                 ggarrange(pools.maf,
                           comb.maf,
                           labels = "B",
                           common.legend = T,
                           legend = "bottom",
                           font.label = list(size = 24, face = "bold"),
                           nrow=2))

tstsP

ggsave("results/aggregated_data/minmaxcov_5_95/Figure_tstv_maf_mac.pdf",
       tstsP,
       device="pdf",
       width=15,
       height=9)


## FOR EACH POOL
missing_data <- read.table("paramTest/aggregated_data/minmaxcov_5_95/snps_all/maf_mac_missing.mincov5.maxcov95.txt", header=TRUE, na.strings = NA)

# Remove Overall data
DATA.missing.pools <- missing_data %>%
  #filter(POP !="Overall" & MAF == 0 | MAF == 0.01 | MAF == 0.05 | MAF == 0.1) 
  filter(POP !="Overall") 

# Plot all pools
pools.missing.mac <- ggplot(DATA.missing.pools, aes(x=MAC, y=Missing, color=factor(MAF))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
  labs(color="MAF") + 
  ylim(0,75) +
  xlab("Minor allele count (MAC)") +
  ylab("% of missing data") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

pools.missing.maf <- ggplot(DATA.missing.pools, aes(x=MAF, y=Missing, color=factor(MAC))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
  labs(color="MAC") + 
  ylim(0,75) +
  xlab("Minor allele frequency (MAF)") +
  ylab("% of missing data") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))



table(DATA.missing.pools[which(DATA.missing.pools$Missing >= 30), 3])
#NW_BIO1_S1 PA_BIO1_S1 PA_BIO4_S1 PA_BIO4_S2 WO_BIO1_S1 WO_BIO4_S3 
#18         18         12         12          3         15

# Filter pools with > 30% missing data
DATA.missing.pools.filtered <- DATA.missing.pools %>%
  filter(POP !="NW_BIO1_S1" & POP !="PA_BIO1_S1" & POP !="PA_BIO4_S1" & POP !="PA_BIO4_S2" & POP !="WO_BIO4_S3")

# Plot only kept pools
pools.missing.mac.filtered <- ggplot(DATA.missing.pools.filtered, aes(x=MAC, y=Missing, color=factor(MAF))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
  labs(color="MAF") + 
  ylim(0,75) +
  xlab("Minor allele count (MAC)") +
  ylab("% of missing data") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

pools.missing.maf.filtered <- ggplot(DATA.missing.pools.filtered, aes(x=MAF, y=Missing, color=factor(MAC))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=FALSE,fullrange=TRUE) +
  labs(color="MAC") + 
  ylim(0,75) +
  xlab("Minor allele frequency (MAF)") +
  ylab("% of missing data") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))


par(mfrow=c(3, 1))
missingP<-ggarrange(ggarrange(pools.missing.mac,
                              pools.missing.mac.filtered,
                              labels = "A",
                              common.legend = T,
                              legend = "bottom",
                              font.label = list(size = 24, face = "bold"),
                              nrow=2),
                    ggarrange(pools.missing.maf,
                              pools.missing.maf.filtered,
                              labels = "B",
                              common.legend = T,
                              legend = "bottom",
                              font.label = list(size = 24, face = "bold"),
                              nrow=2))

missingP

ggsave("results/aggregated_data/minmaxcov_5_95/Figure_missing_maf_mac_pools.allSNPs.pdf",
       missingP,
       device="pdf",
       width=15,
       height=9)



macs.mafs <- read.table("paramTest/aggregated_data/minmaxcov_5_95/snps_all/maf.mac.pnps.mincov5.maxcov95.txt", header=TRUE, na.strings = NA)


# Keep only genomewid data
DATA.pnps.pools.mafs <- macs.mafs %>%
  filter(Chrom =="genomewide")

# Plot 
genomewide.pspn.mac <- ggplot(DATA.pnps.pools.mafs, aes(x=MAC, y=pNpS, color=factor(MAF/10))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=TRUE,fullrange=TRUE) +
  labs(color="MAF") + 
  ylim(0.6,2.2) +
  xlab("Minor allele count (MAC)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

genomewide.pspn.maf <- ggplot(DATA.pnps.pools.mafs, aes(x=(MAF/10), y=pNpS, color=factor(MAC))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=TRUE,fullrange=TRUE) +
  labs(color="MAC") + 
  ylim(0.6,2.2) +
  xlab("Minor allele frequency (MAF)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))


# Keep only genomewide data for pools that passed missing data filtering
DATA.pnps.pools.mafs.filtered <- DATA.pnps.pools.mafs %>%
  filter(POP !="NW_BIO1_S1" & POP !="PA_BIO1_S1" & POP !="PA_BIO4_S1" & POP !="PA_BIO4_S2" & POP !="WO_BIO4_S3")


genomewide.pspn.mac.filtered <- ggplot(DATA.pnps.pools.mafs.filtered, aes(x=MAC, y=pNpS, color=factor(MAF/10))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=TRUE,fullrange=TRUE) +
  labs(color="MAF") + 
  ylim(0.6,2.2) +
  xlab("Minor allele count (MAC)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))

genomewide.pspn.maf.filtered <- ggplot(DATA.pnps.pools.mafs.filtered, aes(x=(MAF/10), y=pNpS, color=factor(MAC))) + 
  geom_point() +
  geom_smooth(method = 'loess', se=TRUE,fullrange=TRUE) +
  labs(color="MAC") + 
  ylim(0.6,2.2) +
  xlab("Minor allele frequency (MAF)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     legend.position = "right", axis.text = element_text(size = 14),
                     axis.title=element_text(size=16))



par(mfrow=c(3, 1))
pnpsP <- ggarrange(ggarrange(genomewide.pspn.mac,
                             genomewide.pspn.mac.filtered,
                             labels = "A",
                             common.legend = T,
                             legend = "bottom",
                             font.label = list(size = 24, face = "bold"),
                             nrow=2),
                   ggarrange(genomewide.pspn.maf,
                             genomewide.pspn.maf.filtered,
                             labels = "B",
                             common.legend = T,
                             legend = "bottom",
                             font.label = list(size = 24, face = "bold"),
                             nrow=2))

pnpsP

ggsave("results/aggregated_data/minmaxcov_5_95/Figure_pnps_maf_mac.allSNPs.pdf",
       pnpsP,
       device="pdf",
       width=15,
       height=9)


# Calculate the mean, sd and se for each MAF and MAC
DATA.pnps.mean.maf.mac <- DATA.pnps.pools.mafs %>%
  group_by(MAF,MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS, na.rm = T),
    pnps.sd=sd(pNpS, na.rm = T),
    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
  )

DATA.pnps.mean.maf.mac.selected.maf <- DATA.pnps.mean.maf.mac %>%
  filter(MAF == 0.01 | MAF == 0.1 | MAF == 0.5)

DATA.pnps.mean.maf.mac.selected.maf$MAF <- (DATA.pnps.mean.maf.mac.selected.maf$MAF/10)

DATA.pnps.mean.maf.mac.selected.mac <- DATA.pnps.mean.maf.mac %>%
  filter(MAC == 5 | MAC == 10 | MAC == 50 | MAC == 100)

# More fancy plot

# All MAC as a function of MAF
genomewide.mean.pspn.mac <- ggplot(DATA.pnps.mean.maf.mac.selected.maf, aes(x=MAC, y=pnps.m, color=factor(MAF))) + 
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.001), aes(x=MAC,y=pNpS)) +
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.01), aes(x=MAC,y=pNpS)) +
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAF==0.05), aes(x=MAC,y=pNpS)) +
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAF)), colour = NA, alpha=0.2) +
  geom_line(lwd=1.5) + 
  xlab("Minor allele count (MAC)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  facet_wrap(vars("PoolSeq all pools")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.position = "right")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
  scale_fill_manual(values=c(alpha(c("#999999"),
                                   alpha=ALPHA),
                             alpha(c("#de77ae"),
                                   alpha=ALPHA),
                             alpha(c("#abd9e9"),
                                   alpha=ALPHA)), name="MAF")+
  scale_colour_manual(values=c(alpha(c("#999999"),
                                     alpha=ALPHA),
                               alpha(c("#de77ae"),
                                     alpha=ALPHA),
                               alpha(c("#abd9e9"),
                                     alpha=ALPHA)), name="MAF")
genomewide.mean.pspn.mac

# All MAF as a function of MAC
genomewide.mean.pspn.maf <- ggplot(DATA.pnps.mean.maf.mac.selected.mac, aes(x=(MAF/10), y=pnps.m, color=factor(MAC))) + 
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==5), aes(x=MAF/10,y=pNpS)) +
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==10), aes(x=MAF/10,y=pNpS)) +
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==50), aes(x=MAF/10,y=pNpS)) +
  geom_jitter(data=filter(macs.mafs,Chrom=="genomewide" & MAC==100), aes(x=MAF/10,y=pNpS)) +
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAC)), colour = NA, alpha=0.2) +
  geom_line(lwd=1.5) + 
  xlab("Minor allele frequency (MAF)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  facet_wrap(vars("PoolSeq all pools")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.position = "right") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
  scale_fill_manual(values=c(alpha(c("#999999"),
                                   alpha=ALPHA),
                             alpha(c("#E69F00"),
                                   alpha=ALPHA),
                             alpha(c("#56B4E9"),
                                   alpha=ALPHA),
                             alpha(c("#009E73"),
                                   alpha=ALPHA)), name="MAC", breaks = c("5", "10", "50", "100"))+
  scale_colour_manual(values=c(alpha(c("#999999"),
                                     alpha=ALPHA),
                               alpha(c("#E69F00"),
                                     alpha=ALPHA),
                               alpha(c("#56B4E9"),
                                     alpha=ALPHA),
                               alpha(c("#009E73"),
                                     alpha=ALPHA)), name="MAC", breaks = c("5", "10", "50", "100"))
genomewide.mean.pspn.maf

# Filtered data
# Calculate the mean, sd and se for each MAF and MAC
DATA.pnps.mean.maf.mac.filtered <- DATA.pnps.pools.mafs.filtered %>%
  group_by(MAF,MAC) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS, na.rm = T),
    pnps.sd=sd(pNpS, na.rm = T),
    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
  )

DATA.pnps.mean.maf.mac.selected.maf.filtered <- DATA.pnps.mean.maf.mac.filtered %>%
  filter(MAF == 0.01 | MAF == 0.1 | MAF == 0.5)

DATA.pnps.mean.maf.mac.selected.maf.filtered$MAF <- (DATA.pnps.mean.maf.mac.selected.maf.filtered$MAF/10)
DATA.pnps.pools.mafs.filtered$MAF <- (DATA.pnps.pools.mafs.filtered$MAF/10)

DATA.pnps.mean.maf.mac.selected.mac.filtered <- DATA.pnps.mean.maf.mac.filtered %>%
  filter(MAC == 5 | MAC == 10 | MAC == 50 | MAC == 100)

# More fancy plot
# All MAC as a function of MAF
genomewide.mean.pspn.mac.filtered <- ggplot(DATA.pnps.mean.maf.mac.selected.maf.filtered, aes(x=MAC, y=pnps.m, color=factor(MAF))) + 
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.001), aes(x=MAC,y=pNpS)) +   
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.1), aes(x=MAC,y=pNpS)) +
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAF==0.5), aes(x=MAC,y=pNpS)) +
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAF)), colour = NA, alpha=0.2) +
  geom_line(lwd=1.5) + 
  xlab("Minor allele count (MAC)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  facet_wrap(vars("PoolSeq pools < 30% missing data")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.position = "right")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
  scale_fill_manual(values=c(alpha(c("#999999"),
                                   alpha=ALPHA),
                             alpha(c("#de77ae"),
                                   alpha=ALPHA),
                             alpha(c("#abd9e9"),
                                   alpha=ALPHA)), name="MAF")+
  scale_colour_manual(values=c(alpha(c("#999999"),
                                     alpha=ALPHA),
                               alpha(c("#de77ae"),
                                     alpha=ALPHA),
                               alpha(c("#abd9e9"),
                                     alpha=ALPHA)), name="MAF")
genomewide.mean.pspn.mac.filtered


# All MAF as a function of MAC
genomewide.mean.pspn.maf.filtered <- ggplot(DATA.pnps.mean.maf.mac.selected.mac.filtered, aes(x=MAF/10, y=pnps.m, color=factor(MAC))) + 
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==5), aes(x=MAF,y=pNpS)) +
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==10), aes(x=MAF,y=pNpS)) +
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==50), aes(x=MAF,y=pNpS)) +
  geom_jitter(data=filter(DATA.pnps.pools.mafs.filtered, MAC==100), aes(x=MAF,y=pNpS)) +
  geom_ribbon(aes(ymin = pnps.m-pnps.sd, ymax = pnps.m+pnps.sd, fill=factor(MAC)), colour = NA, alpha=0.2) +
  geom_line(lwd=1.5) + 
  xlab("Minor allele frequency (MAF)") +
  ylab(expression(italic("p")["N"]/italic("p")["S"])) +
  facet_wrap(vars("PoolSeq pools < 30% missing data")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20, angle = 90),
        axis.title.x = element_text(size = 20, angle = 0),
        axis.text=element_text(size=18),
        strip.text =element_text(size=20),
        legend.position = "right") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05), limits=c(0.6,2.1)) +
  scale_fill_manual(values=c(alpha(c("#999999"),
                                   alpha=ALPHA),
                             alpha(c("#E69F00"),
                                   alpha=ALPHA),
                             alpha(c("#56B4E9"),
                                   alpha=ALPHA),
                             alpha(c("#009E73"),
                                   alpha=ALPHA)), name="MAC", breaks = c("5", "10", "50", "100"))+
  scale_colour_manual(values=c(alpha(c("#999999"),
                                     alpha=ALPHA),
                               alpha(c("#E69F00"),
                                     alpha=ALPHA),
                               alpha(c("#56B4E9"),
                                     alpha=ALPHA),
                               alpha(c("#009E73"),
                                     alpha=ALPHA)), name="MAC", breaks = c("5", "10", "50", "100"))
genomewide.mean.pspn.maf.filtered


par(mfrow=c(3, 1))
mean.pnpsP <- ggarrange(ggarrange(genomewide.mean.pspn.mac,
                                  genomewide.mean.pspn.mac.filtered,
                                  labels = "A",
                                  common.legend = T,
                                  legend = "bottom",
                                  font.label = list(size = 24, face = "bold"),
                                  nrow=2),
                        ggarrange(genomewide.mean.pspn.maf,
                                  genomewide.mean.pspn.maf.filtered,
                                  labels = "B",
                                  common.legend = T,
                                  legend = "bottom",
                                  font.label = list(size = 24, face = "bold"),
                                  nrow=2))

mean.pnpsP


ggsave("results/aggregated_data/minmaxcov_5_95/Figure_fancy_pnps_maf_mac.allSNPs.pdf",
       mean.pnpsP,
       device="pdf",
       width=15,
       height=9)



