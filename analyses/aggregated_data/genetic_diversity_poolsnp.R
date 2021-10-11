##########################################################
#       SNP-SPECIFIC GENETIC DIVERSITY, GENE DIVERSITY   # 
#         pi_N and pi_S, and GENETIC CONSTRAIN           #
##########################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

###
###
### ---- SETUP ANALYZES ENVIRONMENT ----
###
###

### Recover R-renv environment
setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()

### Remove last features
rm(list=ls())
ls()

### Load R libraries
library(tidyverse)
library(gridExtra)
library(scales)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(xtable)

### Import auxiliary R functions  
source("AglyPoolseq/analyses/aggregated_data/aux_func.R")

ALPHA=0.75

poolnames <- c("MN-Av.1", "MN-V.1",
               "ND-Av.1", "ND-V.1", 
               "NW-Av.1", "NW-Av.2", "NW-V.1", "NW-V.2", "NW-V.3", 
               "PA-Av.1", "PA-V.1", "PA-V.2", 
               "WI-Av.1", "WI-Av.2", "WI-V.1", "WI-V.2", "WI-V.3", 
               "WO-Av.1", "WO-V.1", "WO-V.2", "WO-V.3")

###
###
### --- LOAD THE DATASETS ---
###
###

### LOAD THE MATRIX WITH ML-ESTIMATED REFERENCE ALLELE FREQUENCIES 
imputedRefMLFreq = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_pools.txt'
dt.1.imputedRefMLFreq.dataFrame <- read.table(file=imputedRefMLFreq, sep = '\t', header = TRUE)
dim(dt.1.imputedRefMLFreq.dataFrame) # 262866     25

# Check the variable classes
class(dt.1.imputedRefMLFreq.dataFrame$Chromosome) # character
class(dt.1.imputedRefMLFreq.dataFrame$Position) # numeric?

### LOAD POOLFSTAT CALCULATED LOCUS-SPECIFIC FSTs
poolfstat_intralocusFst = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/intralocus_fst_data.txt'
poolfstat_intralocusFst.dataFrame <- read.table(file=poolfstat_intralocusFst, sep = ' ', header = TRUE)
dim(poolfstat_intralocusFst.dataFrame) # 262866      5

# Check the variable classes
class(poolfstat_intralocusFst.dataFrame$Chromosome) # character
class(poolfstat_intralocusFst.dataFrame$Position) # numeric?

###
###
### --- CALCULATE SUMMARY STATISTICS ---
###
###

### COMPUTE OVERALL HETEROZYGOSITY ACROSS SAMPLES - TOTAL HETEROZYGOSITY H_T
###-------------------------------------------------------------------------

## FOR ALL LOCI - ACROSS ALL POOLS

# check the % of missing data of each marker
perc_indMissing_snp <- apply(dt.1.imputedRefMLFreq.dataFrame[,-c(1:4)], 1,  
                             FUN = function(x) 100*(sum(is.na(x))/length(colnames(dt.1.imputedRefMLFreq.dataFrame[,-c(1:4)]))))

max(perc_indMissing_snp) # 47.61905

# Calculate total HE for SNPs with < 50% missing data
# with all samples, no markers has > 50% missing data
total_he <- apply(dt.1.imputedRefMLFreq.dataFrame[,-c(1:4)], 1, totalHE) # 1 indicates rows: locus
length(total_he) # 262866
sum(!is.na(total_he)) # 262866

### FOR ALL LOCI - ACROSS ALL POOLS WITHIN BIOTYPES

# Subset the allele frequency matrix for each biotype
dt.1.imputedRefMLFreq.dataFrame.b1 <- dt.1.imputedRefMLFreq.dataFrame[,-c(6,8,11:13,15,16,19:21,23:25)]
dt.1.imputedRefMLFreq.dataFrame.b4 <- dt.1.imputedRefMLFreq.dataFrame[,-c(5,7,9,10,14,17,18,22)]

# check the % of missing data of each marker
# Biotype 1
perc_indMissing_snp_b1 <- apply(dt.1.imputedRefMLFreq.dataFrame.b1[,-c(1:4)], 1,  
                                FUN = function(x) 100*(sum(is.na(x))/length(colnames(dt.1.imputedRefMLFreq.dataFrame.b1[,-c(1:4)]))))

max(perc_indMissing_snp_b1) # 87.5%

# Biotype 4
perc_indMissing_snp_b4 <- apply(dt.1.imputedRefMLFreq.dataFrame.b4[,-c(1:4)], 1,  
                                FUN = function(x) 100*(sum(is.na(x))/length(colnames(dt.1.imputedRefMLFreq.dataFrame.b4[,-c(1:4)]))))

max(perc_indMissing_snp_b4) # 76.92308%

# Calculate the total HE for all SNPs
total_he_b1 <- apply(dt.1.imputedRefMLFreq.dataFrame.b1[,-c(1:4)], 1, totalHE) # 1 indicates rows: locus
total_he_b4 <- apply(dt.1.imputedRefMLFreq.dataFrame.b4[,-c(1:4)], 1, totalHE) # 1 indicates rows: locus

# Then imput NA for SNPs with > 50% of missing within each biotype
# Biotype 1
total_he_b1[which(perc_indMissing_snp_b1 > 50)] <- NA

length(total_he_b1) # 262866
sum(!is.na(total_he_b1)) # 255592
((262866-255592)/262866)*100 # 2.767189% of the SNPS has >50% missing data

# Biotype 4
total_he_b4[which(perc_indMissing_snp_b4 > 50)] <- NA

length(total_he_b4) # 262866
sum(!is.na(total_he_b4)) # 253743
((262866-253743)/262866)*100 # 3.47059% of the SNPS has >50% missing data

average_total_he <- data.frame(pools  =colnames(cbind(total_he, total_he_b1, total_he_b4)),
                              average=colMeans(cbind(total_he, total_he_b1, total_he_b4), na.rm = T),
                              sd     =apply(cbind(total_he, total_he_b1, total_he_b4), 2, function(x) sd(x, na.rm = TRUE)),
                              se     =apply(cbind(total_he, total_he_b1, total_he_b4), 2, function(x) sd(x, na.rm = TRUE)/sum(!is.na(x))))

### EXPORT AS A LATEX TABLE                              
#print(xtable(average_total_he, type = "latex"), file = "results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/average_total_he_table.tex")


### BARPLOT WITH AVERAGE VALUES
# Prepare the error bar to place in the barplot 
limits_total_he <- aes(ymax = average + sd,
                        ymin = average - sd)

p1 <- ggplot(data = average_total_he, aes(x = pools, y = average))
p1 <- p1 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
            geom_errorbar(limits_total_he, position = position_dodge(0.9), width = 0.25) + 
            labs(x = NULL, y =expression(paste("Average"," ", italic(H)[E]))) + 
            ylim(0, 0.5) + 
            scale_x_discrete(labels = c("Global", "Avirulent", "Virulent")) +
            theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                               axis.line = element_line(colour = "black"),
                               axis.text = element_text(size = 14),
                               axis.title=element_text(size=20),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1

### BOXPLOT ALL VALUES
#total_he_data4plot <- data.frame(pool_id   = c(rep("a", length(total_he)),
#                                                rep("b", length(total_he_b1)),
#                                                rep("c", length(total_he_b4))),
#                                  he_across = c(total_he, 
#                                                total_he_b1, 
#                                                total_he_b4))
#
##p2 <- ggplot(total_he_data4plot, aes(x=pool_id, y=he_across))
#p2 <- p2 + geom_boxplot(width=0.1, fill = "gray", outlier.shape = NA) +
#           scale_x_discrete(labels = c("Global", "Avirulent", "Virulent")) +
#           labs(x = NULL, y =expression(paste("Average"," ", italic(H)[E]))) + 
#           #ylim(-0.01, 0.31) +
#           theme_bw() + theme(panel.border = element_rect(colour = "black"), 
#                              axis.line = element_line(colour = "black"),
#                              axis.text = element_text(size = 14),
#                              axis.title=element_text(size=20))
#p2

#save_checkpoint_1
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

### COMPUTE THE WITHIN-SAMPLE HETEROZYGOSITY AS PI WHITHIN SAMPLE
###--------------------------------------------------------------

### CALCULATE THE WITHIN-SAMPLE HETEROZYGOSITY FOR EACH SAMPLE
pi_within_snp_sample <- apply(dt.1.imputedRefMLFreq.dataFrame[,-c(1:4)], 2, withinLocusPi)

### CALCULATE THE AVERAGES, SD AND SE 
average_pi_within_samples <- data.frame(pools  =poolnames,
                                        biotypes = ifelse(as.numeric(do.call(rbind, strsplit(colnames(pi_within_snp_sample), "", fixed=TRUE))[,5]) == 1, "Avirulent", "Virulent" ),
                                        average=colMeans(pi_within_snp_sample, na.rm = T),
                                        sd     =apply(pi_within_snp_sample, 2, function(x) sd(x, na.rm = TRUE)),
                                        se     =apply(pi_within_snp_sample, 2, function(x) sd(x, na.rm = TRUE)/sum(!is.na(x)))
                                        ) 

### EXPORT AS A LATEX TABLE                              
#print(xtable(average_pi_within_samples, type = "latex"), file = "results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/average_pi_within-samples_table.tex")

### NON-PARAMETRIC TEST FOR DIFFERENCES ON THE WITHIN-SAMPLE HETEROZYGOSISTY BETWEEN SAMPLES
## kruskal.test
# FROM A MATRIX TO A LIST OF VECTORS CONTAINING WITHIN-SAMPLE HETEROZYGOSITIES
pi_within_snp_sample.list <- split(pi_within_snp_sample, rep(1:ncol(pi_within_snp_sample), each = nrow(pi_within_snp_sample)))

# REMOVE MISSING VALUES WITHIN EACH LIST OBJECT
pi_within_snp_sample.list.filteredna <- Map(na.omit, pi_within_snp_sample.list)

# KRUSKAL.TEST ON K SAMPLES
kruskal.test(pi_within_snp_sample.list.filteredna)
# Kruskal-Wallis chi-squared = 45842, df = 20, p-value < 2.2e-16

### NON-PARAMETRIC PAIRWISE TEST WITHIN STATES
## wilcox.test
# MN
wilcox.test(x=pi_within_snp_sample[,1],
            y=pi_within_snp_sample[,2])
#W = 3.1788e+10, p-value < 2.2e-16

# ND
wilcox.test(x=pi_within_snp_sample[,3],
            y=pi_within_snp_sample[,4])
#W = 2.7878e+10, p-value < 2.2e-16

# NW
wilcox.test(x=rowMeans(cbind(pi_within_snp_sample[,5], pi_within_snp_sample[,6]), na.rm = T),
            y=rowMeans(cbind(pi_within_snp_sample[,7], pi_within_snp_sample[,8], pi_within_snp_sample[,9]), na.rm = T))
#W = 3.0135e+10, p-value < 2.2e-16

# PA
wilcox.test(x=pi_within_snp_sample[,10],
            y=rowMeans(cbind(pi_within_snp_sample[,11], pi_within_snp_sample[,12]), na.rm = T))
#W = 8554550830, p-value < 2.2e-16

# WI
wilcox.test(x=rowMeans(cbind(pi_within_snp_sample[,13], pi_within_snp_sample[,14]), na.rm = T),
            y=rowMeans(cbind(pi_within_snp_sample[,15], pi_within_snp_sample[,16], pi_within_snp_sample[,17]), na.rm = T))
#W = 3.0958e+10, p-value < 2.2e-16

# WO
wilcox.test(x=pi_within_snp_sample[,18],
            y=rowMeans(cbind(pi_within_snp_sample[,19], pi_within_snp_sample[,20], pi_within_snp_sample[,21]), na.rm = T))
#W = 2.6046e+10, p-value < 2.2e-16

### BARPLOT
limits_pi_within_samples <- aes(ymax = average + sd,
                                ymin = average - sd)

p3 <- ggplot(data = average_pi_within_samples, aes(x = pools, y = average, fill = factor(biotypes)))
p3 <- p3 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
            geom_errorbar(limits_pi_within_samples , position = position_dodge(0.9), width = 0.25) + 
            labs(x = NULL, y = expression(paste("Within-sample Heterozygosity"," (", italic(pi)["within"], ")"))) + 
            #ylim(0, 0.5) +
            scale_fill_manual(values = c("#247F00","#AB1A53"),
                              name   = "Biotypes",
                              breaks = c("Avirulent", "Virulent"),
                              labels = c("Avirulent", "Virulent")) +
            theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                               axis.line = element_line(colour = "black"),
                               axis.text = element_text(size = 14),
                               axis.title=element_text(size=20),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3
          
### CALCULATE THE AVERAGE WITHIN-SAMPLE HETEROZYGOSITY ACROSS SAMPLES

### FOR ALL LOCI - ACROSS ALL POOLS
average_pi_within_snp_sample <- apply(pi_within_snp_sample, 1, function(x) mean(x, na.rm = T))
length(average_pi_within_snp_sample) # 262866
sum(!is.na(average_pi_within_snp_sample)) #262866

### FOR ALL LOCI - ACROSS ALL POOLS WITHIN BIOTYPES
average_pi_within_snp_sample_b1 <- apply(pi_within_snp_sample[,-c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) mean(x, na.rm = T))
average_pi_within_snp_sample_b4 <- apply(pi_within_snp_sample[, c(2,4,7,8,9,11,12,15,16,17,19,20,21)], 1, function(x) mean(x, na.rm = T))

# CORRECT FOR % OF MISSING DATA WITHIN BIOTYPES
average_pi_within_snp_sample_b1[which(perc_indMissing_snp_b1 > 50)] <- NA

length(average_pi_within_snp_sample_b1) # 262866
sum(!is.na(average_pi_within_snp_sample_b1)) # 255592

average_pi_within_snp_sample_b4[which(perc_indMissing_snp_b4 > 50)] <- NA

length(average_pi_within_snp_sample_b4) # 262866
sum(!is.na(average_pi_within_snp_sample_b4)) # 253743

### PAIRWISE TEST WITHIN STATES BETWEEN BIOTYPES
wilcox.test(x=average_pi_within_snp_sample_b1,
            y=average_pi_within_snp_sample_b4,  na.action="na.omit")
#W = 3.2571e+10, p-value = 0.006107

wilcox.test(x=average_pi_within_snp_sample,
            y=average_pi_within_snp_sample_b1,  na.action="na.omit")
#W = 3.3139e+10, p-value < 2.2e-16

wilcox.test(x=average_pi_within_snp_sample,
            y=average_pi_within_snp_sample_b4,  na.action="na.omit")
#W = 3.3068e+10, p-value = 1.426e-07

## PREPARE THE DATA FOR BARPLOT BIOTYPES
average_pi_within_snp_sample_table <- data.frame(pools  =colnames(cbind(average_pi_within_snp_sample, average_pi_within_snp_sample_b1, average_pi_within_snp_sample_b4)),
                                                 average=colMeans(cbind(average_pi_within_snp_sample, average_pi_within_snp_sample_b1, average_pi_within_snp_sample_b4), na.rm = T),
                                                 sd     =apply(cbind(average_pi_within_snp_sample, average_pi_within_snp_sample_b1, average_pi_within_snp_sample_b4), 2, function(x) sd(x, na.rm = TRUE)),
                                                 se     =apply(cbind(average_pi_within_snp_sample, average_pi_within_snp_sample_b1, average_pi_within_snp_sample_b4), 2, function(x) sd(x, na.rm = TRUE)/sum(!is.na(x)))
                                                 )        
### EXPORT AS A LATEX TABLE                              
#print(xtable(average_pi_within_snp_sample_table, type = "latex"), file = "results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/average_pi_within-snp-samples_table.tex")

### BARPLOT WITH AVERAGE VALUES
# Prepare the error bar to place in the barplot 
limits_pi_within_snp_samples <- aes(ymax = average + sd,
                                    ymin = average - sd)

p4 <- ggplot(data = average_pi_within_snp_sample_table, aes(x = pools, y = average))
p4 <- p4 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
            geom_errorbar(limits_pi_within_snp_samples, position = position_dodge(0.9), width = 0.25) + 
            labs(x = NULL, y =expression(paste("Average"," ", italic(pi)["within"]))) + 
            ylim(0, 0.5) + 
            scale_x_discrete(labels = c("Global", "Avirulent", "Virulent")) +
            theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                               axis.line = element_line(colour = "black"),
                               axis.text = element_text(size = 14),
                               axis.title=element_text(size=20),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4

## BOXPLOT ALL VALUES
#pi_within_snp_samples_data4plot <- data.frame(pool_id   = c(rep("a", length(average_pi_within_snp_sample)),
#                                                            rep("b", length(average_pi_within_snp_sample_b1)),
#                                                            rep("c", length(average_pi_within_snp_sample_b4))),
#                                          avg_pi_within = c(average_pi_within_snp_sample, 
#                                                            average_pi_within_snp_sample_b1, 
#                                                            average_pi_within_snp_sample_b4))
#         
#p5 <- ggplot(pi_within_snp_samples_data4plot, aes(x=pool_id, y=avg_pi_within))
#p5 <- p5 + geom_boxplot(width=0.1, fill = "gray", outlier.shape = NA) +
#           scale_x_discrete(labels = c("Global", "Avirulent", "Virulent")) +
#           labs(x = NULL, y =expression(paste("Average"," ", italic(pi)["within"]))) + 
#           ylim(0, 0.5) +
#           theme_bw() + theme(panel.border = element_rect(colour = "black"), 
#                              axis.line = element_line(colour = "black"),
#                              axis.text = element_text(size = 14),
#                              axis.title=element_text(size=20))
#p5

### COMBINE PLOTS 3 AND 4
par(mfrow=c(1, 2))
avg_pi_within_combined_plot <- ggarrange(ggarrange(p3,
                                         labels = "A",
                                         font.label = list(size = 24, face = "bold"),
                                         nrow=1, legend = "top"),
                               ggarrange(p4,
                                         labels = "B",
                                         font.label = list(size = 24, face = "bold"),
                                         nrow=1))

avg_pi_within_combined_plot

ggsave("results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/avg_pi_within_samples_biotypes.pdf",
       avg_pi_within_combined_plot,
       device="pdf",
       width=15,
       height=8.5)

#save_checkpoint_2
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

###
###
### --- GENETICS DIVERSITY AND FST ---
###
###

### FOR ALL LOCI - ACROSS ALL POOLS
genetic_diversity <- geneticDiversity(x=dt.1.imputedRefMLFreq.dataFrame, Nindv = 5) 

### Average for biotype samples
genetic_diversity.b1 <- geneticDiversity(x=dt.1.imputedRefMLFreq.dataFrame.b1, Nindv = 5) 
genetic_diversity.b4 <- geneticDiversity(x=dt.1.imputedRefMLFreq.dataFrame.b4, Nindv = 5) 

### CONSERVATIVE APPROACH: IMPUTE NA IN LOCUS WITH > 50% MISSING DATA 

# FOR BIOTYPE 1
genetic_diversity.b1$pi_within[which(perc_indMissing_snp_b1 > 50)] <- NA
genetic_diversity.b1$pi_between[which(perc_indMissing_snp_b1 > 50)] <- NA
genetic_diversity.b1$H_T[which(perc_indMissing_snp_b1 > 50)] <- NA
genetic_diversity.b1$fst[which(perc_indMissing_snp_b1 > 50)] <- NA

# Check if NA imputation worked
length(genetic_diversity.b1$pi_within) # 262866
sum(!is.na(genetic_diversity.b1$pi_within)) #255592

# FOR BIOTYPE 4
genetic_diversity.b4$pi_within[which(perc_indMissing_snp_b4 > 50)] <- NA
genetic_diversity.b4$pi_between[which(perc_indMissing_snp_b4 > 50)] <- NA
genetic_diversity.b4$H_T[which(perc_indMissing_snp_b4 > 50)] <- NA
genetic_diversity.b4$fst[which(perc_indMissing_snp_b4 > 50)] <- NA

# Check if NA imputation worked
length(genetic_diversity.b4$pi_within) # 262866
sum(!is.na(genetic_diversity.b4$pi_within)) #253743

### COMBINE ALL SUMMARY STATISTICS IN A TABLE
summary_statistics_table <- data.frame(genetic_diversity$snpInfo, 
                                       pi_within=genetic_diversity$pi_within, pi_between=genetic_diversity$pi_between,
                                       H_T=genetic_diversity$H_T, WC_Beta=genetic_diversity$fst, 
                                       POOLFSTAT_FST=poolfstat_intralocusFst.dataFrame$fst)

summary_statistics_table.b1 <- data.frame(genetic_diversity.b1$snpInfo, 
                                          pi_within=genetic_diversity.b1$pi_within, pi_between=genetic_diversity.b1$pi_between,
                                          H_T=genetic_diversity.b1$H_T, WC_Beta=genetic_diversity.b1$fst)

summary_statistics_table.b4 <- data.frame(genetic_diversity.b4$snpInfo, 
                                          pi_within=genetic_diversity.b4$pi_within, pi_between=genetic_diversity.b4$pi_between,
                                          H_T=genetic_diversity.b4$H_T, WC_Beta=genetic_diversity.b4$fst)

# change Position to integer
summary_statistics_table$Position <- as.integer(summary_statistics_table$Position)
summary_statistics_table.b1$Position <- as.integer(summary_statistics_table.b1$Position)
summary_statistics_table.b4$Position <- as.integer(summary_statistics_table.b4$Position)

#save_checkpoint_3
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

###
###
### --- CORRELATION BETWEEN FST ESTIMATES ---
###
###

### CHECK IF POOLFSTAT WEIR & COCKERHAM's FST (1984) CORRELATES 
### WITH COCKERHAM & WEIR BETA (1993) (From the average between and within locus pi)

fst_data4plot <- data.frame(x=summary_statistics_table$POOLFSTAT_FST,
                            y=summary_statistics_table$WC_Beta) 


fst.cor <- cor(fst_data4plot$x, fst_data4plot$y, use = "pairwise.complete.obs")
print(fst.cor)# 0.9071797

par(mar=c(5,5,4,1)+.1)
breaks.plot <- c(54.598185, 7.389056, 1)
fst.plot   <- ggplot(fst_data4plot, aes(x,y))
fst.plot <- fst.plot + 
            stat_bin2d(bins=100, na.rm = T, show.legend=T) + 
            scale_fill_gradientn(colours=viridis::viridis(32), trans="log", name="counts"
                                 ,
                                 breaks = round(breaks.plot),labels = round(breaks.plot),
            ) +
            xlab(expression(paste("Weir & Cockerham's", " ", italic(F)[ST]))) +
            ylab(expression(paste("Cockeram & Weir's", " ", italic(beta))))+
            xlim(-0.2,0.6) +
            ylim(-0.2,0.6) +
            geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") +
            annotate(geom="text", x=-0.16, y=0.6, size= 5, label=expression(paste(italic(r), " = ", "0.907"))) + 
            theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                               legend.position = "right", axis.text = element_text(size = 14),
                               axis.title=element_text(size=16))
fst.plot

# CORRELATION PLOT
ggsave("results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/correlation_among_fsts_estimates.pdf",
       fst.plot,
       device="pdf",
       width=7.6,
       height=6.6)

#save_checkpoint_4
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")


###
###
### --- SNP ANNOTATION ---
###
###

### LOAD SNP ANNOTATION TABLE
### IT CONTAINS ALL SNPS PRODUCED BY THE PIPELINE BUT BEFORE POOLFSTAT FILTERING
snpsANN          = 'vcf/aggregated_data/minmaxcov_4_99/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.vcf.gz.info.tsv'
snpsANN_header <- c('CHROM', 'POS', 'REF_ALT', 'ADP', 'DP', 'NC', 'AC', 'AF', 'ANN')

snps_ann_table <- read.table(file = snpsANN, sep = '\t', col.names = snpsANN_header)

## Extract only the first annotation of each SNP in the snps_ann_table
raw_annotation <- do.call('rbind', strsplit(snps_ann_table$ANN, '[|]' ))
annotation_names <- c('MUTATION', 'TYPE', 'IMPACT', 'GENE', "REGION", "FUNCTION", "NT_CHANGE", "AA_CHANGE")
annotation <- raw_annotation[,c(1,2,3,4,6,8,10,11)]
colnames(annotation) <- annotation_names

## Add a NS or SS label for missense_variant OR synonymous_variant
fitness_effect  <- ifelse(annotation[,2] == "missense_variant" | annotation[,2] == "missense_variant&splice_region_variant", "NS", 
                          ifelse(annotation[,2] == "synonymous_variant" | annotation[,2] == "splice_region_variant&synonymous_variant", "SS", NA))

## Combine this annotation with the rest of the res table
final_annotation <- data.frame(snps_ann_table[,-c(9)], annotation, FITNES_EFFECT=fitness_effect)
dim(final_annotation) # 268589     17

## KEEP THE SAME SNPS AS IN THE ML IMPUTED ALLELE FREQUENCY TABLE
final_annotation <- final_annotation[paste0(final_annotation$CHROM, "_", final_annotation$POS) 
                                     %in% paste0(summary_statistics_table$Chromosome, "_", summary_statistics_table$Position), ]

dim(final_annotation) # 262866     17

# SANITY CHECK: CHROM AND POS ARE THE SAME IN summary_statistics_table and commonsnps_annotation
sum(summary_statistics_table$Chromosome == final_annotation$CHROM)  #262866
sum(summary_statistics_table$Position == final_annotation$POS) # 262866

sum(paste0(summary_statistics_table$Chromosome, "_", summary_statistics_table$Position
           ) == paste0(final_annotation$CHROM,"_",final_annotation$POS)) # 262866

setdiff(paste0(summary_statistics_table$Chromosome, "_", summary_statistics_table$Position),
        paste0(final_annotation$CHROM,"_",final_annotation$POS)) # 0

###
###
### --- COMBINED SNP SUMMARY STATISTICS WITH SNP ANNOTATION ---
###
###

summary_statistics_annotation_table <- data.frame(final_annotation, summary_statistics_table[,-c(1:4)])

summary_statistics_annotation_table.b1 <- data.frame(final_annotation, summary_statistics_table.b1[,-c(1:4)])
summary_statistics_annotation_table.b4 <- data.frame(final_annotation, summary_statistics_table.b4[,-c(1:4)])

#save_checkpoint_5
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

###
###
### --- GENE pi_N and pi_S ---
###
###

## TAKE ONLY SYNONYMOUS AND NON-SYNONYMOUS MUTATIONS
ss_ns_mutations <- summary_statistics_annotation_table %>%
  filter(FITNES_EFFECT == "SS" | FITNES_EFFECT == "NS")

ss_ns_mutations.b1 <- summary_statistics_annotation_table.b1 %>%
  filter(FITNES_EFFECT == "SS" | FITNES_EFFECT == "NS")

ss_ns_mutations.b4 <- summary_statistics_annotation_table.b4 %>%
  filter(FITNES_EFFECT == "SS" | FITNES_EFFECT == "NS")
  
## CALCULATE MEAN GENE pi_N and pi_S

## FOR ALL LOCI - ACROSS ALL POOLS
# pN/piS from the average within pi for all samples
gene_piN_piS <- ss_ns_mutations %>%
  group_by(GENE) %>%
  summarise(
    # GLOBAL
    nNS = length(pi_within[FITNES_EFFECT=="NS"]),
    nSS = length(pi_within[FITNES_EFFECT=="SS"]),
    piNwi = mean(pi_within[FITNES_EFFECT=="NS"], na.rm = T),
    piSwi = mean(pi_within[FITNES_EFFECT=="SS"], na.rm = T),
    piNbt =  mean(pi_between[FITNES_EFFECT=="NS"], na.rm = T),
    piSbt =  mean(pi_between[FITNES_EFFECT=="SS"], na.rm = T)
   )

dim(gene_piN_piS) # 7473    7

# Make a copy to export the table with all genes (not only the filtered)
# ADDED 17/SEPT/21
gene_piN_piS_copy <- gene_piN_piS

# Remove genes with < 2 nonsynonymous AND < 4 synonymous mutations
gene_piN_piS <- gene_piN_piS[which(gene_piN_piS$nNS >=2 & gene_piN_piS$nSS >= 4), ]
dim(gene_piN_piS) # 433 7 

# Remove genes with piN within AND piS within equals to NA
gene_piN_piS <- gene_piN_piS[!(is.na(gene_piN_piS$piNwi) & is.na(gene_piN_piS$piSwi)), ]
dim(gene_piN_piS) # 433  7

# Remove genes with piN between AND piS between equals to NA
gene_piN_piS <- gene_piN_piS[!(is.na(gene_piN_piS$piNbt) & is.na(gene_piN_piS$piSbt)), ]
dim(gene_piN_piS) # # 433  7

# Add NA to piS within = 0
gene_piN_piS[which(gene_piN_piS$piSwi == 0), 'piSwi'] <- NA

# CALCULATE piN/piS within 
gene_piN_piS <- mutate(gene_piN_piS, piNpiSwi = piNwi/piSwi)

# Add NA to piS between = 0
gene_piN_piS[which(gene_piN_piS$piSbt == 0), 'piSbt'] <- NA

## CALCULATE piN/piS between 
gene_piN_piS <- mutate(gene_piN_piS, piNpiSbt = piNbt/piSbt)

# CALCULATE NI and alpha
gene_piN_piS <- mutate(gene_piN_piS, NI = piNpiSbt/piNpiSwi, alpha = (1 - NI))

## BIOTYPE 1
## pN/piS from the average within pi for all samples
gene_piN_piS.b1 <- ss_ns_mutations.b1 %>%
  group_by(GENE) %>%
  summarise(
    # GLOBAL
    nNS = length(pi_within[FITNES_EFFECT=="NS"]),
    nSS = length(pi_within[FITNES_EFFECT=="SS"]),
    piNwi = mean(pi_within[FITNES_EFFECT=="NS"], na.rm = T),
    piSwi = mean(pi_within[FITNES_EFFECT=="SS"], na.rm = T),
    piNbt =  mean(pi_between[FITNES_EFFECT=="NS"], na.rm = T),
    piSbt =  mean(pi_between[FITNES_EFFECT=="SS"], na.rm = T),
  )

dim(gene_piN_piS.b1) # 7473    7

# Remove genes with < 2 nonsynonymous AND < 4 synonymous mutations
gene_piN_piS.b1 <- gene_piN_piS.b1[which(gene_piN_piS.b1$nNS >=2 & gene_piN_piS.b1$nSS >= 4), ]
dim(gene_piN_piS.b1) # 433  10

# Remove genes with piN within AND piS within equals to NA
gene_piN_piS.b1 <- gene_piN_piS.b1[!(is.na(gene_piN_piS.b1$piNwi) & is.na(gene_piN_piS.b1$piSwi)), ]
dim(gene_piN_piS.b1) # 433   7

# Remove genes with piN between AND piS within equals to NA
gene_piN_piS.b1 <- gene_piN_piS.b1[!(is.na(gene_piN_piS.b1$piNbt) & is.na(gene_piN_piS.b1$piSbt)), ]
dim(gene_piN_piS.b1) # 433   7

# Add NA to piS within = 0
gene_piN_piS.b1[which(gene_piN_piS.b1$piSwi == 0), 'piSwi'] <- NA

# CALCULATE piN/piS within 
gene_piN_piS.b1 <- mutate(gene_piN_piS.b1, piNpiSwi = piNwi/piSwi)

# Add NA to piS between = 0
gene_piN_piS.b1[which(gene_piN_piS.b1$piSbt == 0), 'piSbt'] <- NA

# CALCULATE piN/piS between 
gene_piN_piS.b1 <- mutate(gene_piN_piS.b1, piNpiSbt = piNbt/piSbt)

# CALCULATE NI
gene_piN_piS.b1 <- mutate(gene_piN_piS.b1, NI = piNpiSbt/piNpiSwi, alpha = (1 - NI))

## BIOTYPE 4
## pN/piS from the average within pi for all samples
gene_piN_piS.b4 <- ss_ns_mutations.b4 %>%
  group_by(GENE) %>%
  summarise(
    # GLOBAL
    nNS = length(pi_within[FITNES_EFFECT=="NS"]),
    nSS = length(pi_within[FITNES_EFFECT=="SS"]),
    piNwi = mean(pi_within[FITNES_EFFECT=="NS"], na.rm = T),
    piSwi = mean(pi_within[FITNES_EFFECT=="SS"], na.rm = T),
    piNbt =  mean(pi_between[FITNES_EFFECT=="NS"], na.rm = T),
    piSbt =  mean(pi_between[FITNES_EFFECT=="SS"], na.rm = T),
  )

dim(gene_piN_piS.b4) # 7473    11

# Remove genes with < 2 nonsynonymous AND < 4 synonymous mutations
gene_piN_piS.b4 <- gene_piN_piS.b4[which(gene_piN_piS.b4$nNS >=2 & gene_piN_piS.b4$nSS >= 4), ]
dim(gene_piN_piS.b4) # 433  10

# Remove genes with piN within AND piS within equals to NA
gene_piN_piS.b4 <- gene_piN_piS.b4[!(is.na(gene_piN_piS.b4$piNwi) & is.na(gene_piN_piS.b4$piSwi)), ]
dim(gene_piN_piS.b4) # 433  10

## Remove genes with piN between AND piS within equals to NA
gene_piN_piS.b4 <- gene_piN_piS.b4[!(is.na(gene_piN_piS.b4$piNbt) & is.na(gene_piN_piS.b4$piSbt)), ]
dim(gene_piN_piS.b4) # 433  10

# Add NA to piS within = 0
gene_piN_piS.b4[which(gene_piN_piS.b4$piSwi == 0), 'piSwi'] <- NA

# CALCULATE piN/piS within 
gene_piN_piS.b4 <- mutate(gene_piN_piS.b4, piNpiSwi = piNwi/piSwi)

# Add NA to piS between = 0
gene_piN_piS.b4[which(gene_piN_piS.b4$piSbt == 0), 'piSbt'] <- NA

# CALCULATE piN/piS between 
gene_piN_piS.b4 <- mutate(gene_piN_piS.b4, piNpiSbt = piNbt/piSbt)

# CALCULATE NI
gene_piN_piS.b4 <- mutate(gene_piN_piS.b4, NI = piNpiSbt/piNpiSwi, alpha = (1 - NI))

## PREPARE DATA FOR BARPLOT
dim(gene_piN_piS) # 433  11
dim(gene_piN_piS.b1) # 433  11
dim(gene_piN_piS.b4) # 433  11


### PREPARE THE piN/piS TABLE FOR BOXPLOT
piNpiS_NI_data4plot <- data.frame(pool_id = c(rep("b1", length(gene_piN_piS.b1$GENE)),
                                              rep("b4", length(gene_piN_piS.b4$GENE))),
                                  piSwi  = c(gene_piN_piS.b1$piSwi,
                                             gene_piN_piS.b4$piSwi),
                                  piNpiS = c(gene_piN_piS.b1$piNpiSwi,
                                             gene_piN_piS.b4$piNpiSwi),
                                  NI     = c(gene_piN_piS.b1$NI,
                                             gene_piN_piS.b4$NI),
                                  alpha  = c(gene_piN_piS.b1$alpha,
                                             gene_piN_piS.b4$alpha)
                                  )

## piS within
p6 <- ggplot(piNpiS_NI_data4plot, aes(x=pool_id, y=piSwi, fill=pool_id))
p6 <- p6 + geom_violin(trim=TRUE, scale = "width", adjust = .5)
p6 <- p6 + geom_boxplot(width=0.1, color="black", alpha=0.2) +  
           scale_fill_manual(values = c("#247F00","#AB1A53"),
                             name   = "Biotypes",
                             breaks = c("b1", "b4"),
                             labels = c(expression("Avirulent"), "Virulent")) + 
           scale_x_discrete(labels = c("", "")) +
           xlab(expression(paste(italic(pi)["within"]))) +
           ylab("") +
           ylim(c(0, 0.6)) +
           theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                              axis.line = element_line(colour = "black"),
                              axis.text = element_text(size = 14),
                              axis.title=element_text(size=20),
                              legend.text.align = 0) +
          stat_compare_means(na.rm = TRUE)
p6

## piNpiS within
p7 <- ggplot(piNpiS_NI_data4plot, aes(x=pool_id, y= piNpiS, fill=pool_id))
p7 <- p7 + geom_violin(trim=TRUE, scale = "width", adjust = .5)
p7 <- p7 + geom_boxplot(width=0.1, color="black", alpha=0.2) +  
           scale_fill_manual(values = c("#247F00","#AB1A53"),
                             name   = "Biotypes",
                             breaks = c("b1", "b4"),
                             labels = c(expression("Avirulent"), "Virulent")) + 
           scale_x_discrete(labels = c("", "")) +
           xlab(expression(paste(italic(pi)[N]/italic(pi)[S]))) +
           ylab("") +
           theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                              axis.line = element_line(colour = "black"),
                              axis.text = element_text(size = 14),
                              axis.title=element_text(size=20),
                              legend.text.align = 0) +
           stat_compare_means(na.rm = TRUE)
p7

## NI
p8 <- ggplot(piNpiS_NI_data4plot, aes(x=pool_id, y= NI, fill=pool_id))
p8 <- p8 + geom_violin(trim=TRUE, scale = "width", adjust = .5)
p8 <- p8 + geom_boxplot(width=0.1, color="black", alpha=0.2) +  
      scale_fill_manual(values = c("#247F00","#AB1A53"),
                        name   = "Biotypes",
                        breaks = c("b1", "b4"),
                        labels = c(expression("Avirulent"), "Virulent")) + 
      scale_x_discrete(labels = c("", "")) +
      xlab(expression(italic(NI))) +
      ylab("") +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 14),
                         axis.title=element_text(size=20),
                         legend.text.align = 0) +
      stat_compare_means(na.rm = TRUE)
p8

## ALPHA
p9 <- ggplot(piNpiS_NI_data4plot, aes(x=pool_id, y= alpha, fill=pool_id))
p9 <- p9 + geom_violin(trim=TRUE, scale = "width", adjust = .5)
p9 <- p9 + geom_boxplot(width=0.1, color="black", alpha=0.2) +  
           scale_fill_manual(values = c("#247F00","#AB1A53"),
                             name   = "Biotypes",
                             breaks = c("b1", "b4"),
                             labels = c(expression("Avirulent"), "Virulent")) + 
           scale_x_discrete(labels = c("", "")) +
           xlab(expression(italic(alpha))) +
           ylab("") +
           theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                              axis.line = element_line(colour = "black"),
                              axis.text = element_text(size = 14),
                              axis.title=element_text(size=20),
                              legend.text.align = 0) +
           stat_compare_means(na.rm = TRUE)
p9

## COMBINE PLOTS 6, 7, 8, 9
par(mfrow=c(1, 2))
evolution_plot <- ggarrange(p6,p7,p8, p9,
                          #labels = c("A", "B", "C", "D"),
                          font.label = list(size = 24, face = "bold"),
                          nrow=2, ncol=2, legend = "bottom", common.legend = T)
      
evolution_plot

ggsave("results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/piSwithin_piNpiS_NI_alpha_biotypes.pdf",
       evolution_plot,
       device="pdf",
       width=11,
       height=8.5)

#save_checkpoint_6
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

### Scatterplo pN/pS
#plot(gene_piN_piS.b1$piNpiSwi, gene_piN_piS.b4$piNpiSwi,
#     xlim = c(0,6), ylim = c(0,6))
#abline(a=0, b=1, col="green")
#
#gene_piN_piS.b1[which(gene_piN_piS.b1$piNpiSwi > 4), ]

###
###
### --- EXPORT TABLES ---
###
###

### COMBINED GENOME-WIDE SNP ANNOTATION AND GENETIC DIVERSITY
dim(summary_statistics_annotation_table) #262866     22
dim(summary_statistics_annotation_table.b1) #262866     22
dim(summary_statistics_annotation_table.b4) #262866     22

## SANIT CHECK IF ALL TABLES HAVE THE SAME SNP ORDER
sum(paste0(summary_statistics_annotation_table$CHROM,"_", summary_statistics_annotation_table$POS) == 
      paste0(summary_statistics_annotation_table.b1$CHROM,"_", summary_statistics_annotation_table.b1$POS))

sum(paste0(summary_statistics_annotation_table$CHROM,"_", summary_statistics_annotation_table$POS) == 
      paste0(summary_statistics_annotation_table.b4$CHROM,"_", summary_statistics_annotation_table.b4$POS))


sum(paste0(summary_statistics_annotation_table.b1$CHROM,"_", summary_statistics_annotation_table.b1$POS) == 
      paste0(summary_statistics_annotation_table.b4$CHROM,"_", summary_statistics_annotation_table.b4$POS))

## CHANGE VARIABLE NAMES FOR THE DIVERSITY STATISTICS OF EACH BIOTYPE TABLE
# first make a copy of the each table
summary_statistics_annotation_table.b1.copy <- summary_statistics_annotation_table.b1[,c(18:21)]
summary_statistics_annotation_table.b4.copy <- summary_statistics_annotation_table.b4[,c(18:21)]

# then change the stats names (add .bx to pi_within, pi_between, H_T, and WC_Beta)
colnames(summary_statistics_annotation_table.b1.copy) <- c("pi_within_Av", "pi_between_Av", "H_T_Av", "WC_Beta_Av")
colnames(summary_statistics_annotation_table.b4.copy) <- c("pi_within_Vi", "pi_between_Vi", "H_T_Vi", "WC_Beta_Vi")

combined_summary_statistics_annotation_table <- cbind(summary_statistics_annotation_table[,c(1:17)],
                                                      summary_statistics_annotation_table.b1.copy,
                                                      summary_statistics_annotation_table.b4.copy,
                                                      summary_statistics_annotation_table[,c(18:22)])
head(combined_summary_statistics_annotation_table)

# Save to a file
#write.table(combined_summary_statistics_annotation_table,
#            file="results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/combined_summary_statistics_annotation_table.txt", sep = "\t",
#            quote = F, row.names = F, col.names = TRUE)

### COMBINED GENIC GENETIC DIVERSITY AND ADAPTIVE CONSTRAINS
dim(gene_piN_piS) # 433  11
dim(gene_piN_piS.b1) # 433  11
dim(gene_piN_piS.b4) # 433  11

## SANIT CHECK IF ALL TABLES HAVE THE SAME GENE ORDER
sum(gene_piN_piS$GENE == gene_piN_piS.b1$GENE)
sum(gene_piN_piS$GENE == gene_piN_piS.b4$GENE)
sum(gene_piN_piS.b1$GENE == gene_piN_piS.b4$GENE)

## CHANGE VARIABLE NAMES FOR THE DIVERSITY STATISTICS OF EACH BIOTYPE TABLE
# first make a copy of the each table
gene_piN_piS.b1.copy <- gene_piN_piS.b1[,c(4:11)]
gene_piN_piS.b4.copy <- gene_piN_piS.b4[,c(4:11)]

# then change the stats names 
colnames(gene_piN_piS.b1.copy) <- c("piNwi_Av", "piSwi_Av", "piNbt_Av", "piSbt_Av", "piNpiSwi_Av", "piNpiSbt_Av", "NI_Av", "alpha_Av")
colnames(gene_piN_piS.b4.copy) <- c("piNwi_Vi", "piSwi_Vi", "piNbt_Vi", "piSbt_Vi", "piNpiSwi_Vi", "piNpiSbt_Vi", "NI_Vi", "alpha_Vi")

combined_gene_piN_piS <- cbind(gene_piN_piS[,c(1:3)],
                               gene_piN_piS.b1.copy,
                               gene_piN_piS.b4.copy,
                               gene_piN_piS[,c(4:11)])

head(combined_gene_piN_piS)

# Save to a file
#write.table(combined_gene_piN_piS,
#            file="results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/combined_gene_piN_piS.txt", sep = "\t",
#            quote = F, row.names = F, col.names = TRUE)


### GLOBAL GENIC GENETIC DIVERSITY AND ADAPTIVE CONSTRAINS FOR ALL POSSIBLE GENES
##### REPETE ABOVE TO EXPORT THE FULL/NULL TABLE (COMMENT AFTER EXPORT TABLE)
# ADDED 17/SEPT/21
# THE PART THAT KEEP GENES WITH # OF SS > 4 & NS > 2 WAS REMOVED
# Remove genes with piN within AND piS within equals to NA
gene_piN_piS_copy <- gene_piN_piS_copy[!(is.na(gene_piN_piS_copy$piNwi) & is.na(gene_piN_piS_copy$piSwi)), ]
dim(gene_piN_piS_copy) # 7473    7

# Remove genes with piN between AND piS between equals to NA
gene_piN_piS_copy <- gene_piN_piS_copy[!(is.na(gene_piN_piS_copy$piNbt) & is.na(gene_piN_piS_copy$piSbt)), ]
dim(gene_piN_piS_copy) # 7473    7

# Add NA to piS within = 0
gene_piN_piS_copy[which(gene_piN_piS_copy$piSwi == 0), 'piSwi'] <- NA

# CALCULATE piN/piS within 
gene_piN_piS_copy <- mutate(gene_piN_piS_copy, piNpiSwi = piNwi/piSwi)

# Add NA to piS between = 0
gene_piN_piS_copy[which(gene_piN_piS_copy$piSbt == 0), 'piSbt'] <- NA

### CALCULATE piN/piS between 
gene_piN_piS_copy <- mutate(gene_piN_piS_copy, piNpiSbt = piNbt/piSbt)

## CALCULATE NI and alpha
gene_piN_piS_copy <- mutate(gene_piN_piS_copy, NI = piNpiSbt/piNpiSwi, alpha = (1 - NI))
dim(gene_piN_piS_copy) #  7473   11

# Save to a file
#write.table(gene_piN_piS_copy,
#            file="results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/global_gene_piN_piS.txt", sep = "\t",
#            quote = F, row.names = F, col.names = TRUE)
#####

#save_checkpoint_7
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")


###
###
### --- GENETIC LOAD ---
###
### 

### KEEP ONLY SNPS WITH FUNCTION == protein_coding

# ANNOTATE SNP SUMMARY STATISTICS TABLE
summary_statistics_annotation_table.codingsnps <- summary_statistics_annotation_table[which(summary_statistics_annotation_table$FUNCTION == "protein_coding"), ]

head(summary_statistics_annotation_table.codingsnps)


# REFERENCE ALLELE FREQUENCY TABLE
dt.1.imputedRefMLFreq.dataFrame.codingsnps <- dt.1.imputedRefMLFreq.dataFrame[which(summary_statistics_annotation_table$FUNCTION == "protein_coding"), ]
dt.1.imputedRefMLFreq.dataFrame.codingsnps$Position <- as.integer(dt.1.imputedRefMLFreq.dataFrame.codingsnps$Position)

head(dt.1.imputedRefMLFreq.dataFrame.codingsnps)

# Sanity check if the snps are in the same order in both datasets
dim(summary_statistics_annotation_table.codingsnps) # 151322     22
dim(dt.1.imputedRefMLFreq.dataFrame.codingsnps) # 151322     25

sum(paste0(summary_statistics_annotation_table.codingsnps$CHROM,"_", summary_statistics_annotation_table.codingsnps$POS) == 
      paste0(dt.1.imputedRefMLFreq.dataFrame.codingsnps$Chromosome,"_", dt.1.imputedRefMLFreq.dataFrame.codingsnps$Position))

### CALCULATE POTENTIAL LOAD
## Potential Load: total number of mutations of impact class i in pool j / total number of genic mutations in population k

# FIRST DEFINE A MAF THRESHOULD TO DEFINE IF THE ALTERNATIVE ALLELE IS PRESENT IN THE POOL
maf_potential = 0.05

# Then assign 1 for MAF >= maf_potential, 0 otherwise
potential_load <- apply(dt.1.imputedRefMLFreq.dataFrame.codingsnps[,-c(1:4)], 2, function(x) ifelse(1-x >= maf_potential, 1, 0))

# Combine with the SNP ANNOTATION IMPACT COLUMN
potential_load <- data.frame(IMPACT= summary_statistics_annotation_table.codingsnps$IMPACT,
                             potential_load)

low_potential_load <- apply(potential_load[, -c(1)], 2, function(x) sum(x[potential_load$IMPACT == "LOW"], na.rm = T)/sum(!is.na(x)))
moderate_potential_load <- apply(potential_load[, -c(1)], 2, function(x) sum(x[potential_load$IMPACT == "MODERATE"], na.rm = T)/sum(!is.na(x)))
high_potential_load <- apply(potential_load[, -c(1)], 2, function(x) sum(x[potential_load$IMPACT == "HIGH"], na.rm = T)/sum(!is.na(x)))
modifier_potential_load <- apply(potential_load[, -c(1)], 2, function(x) sum(x[potential_load$IMPACT == "MODIFIER"], na.rm = T)/sum(!is.na(x)))

potential_load_pool <- data.frame(pool=names(low_potential_load),
                                  biotype=do.call(rbind, strsplit(names(low_potential_load), split = '.', fixed = T))[,2],
                                  LOW=low_potential_load,
                                  MODERATE=moderate_potential_load,
                                  HIGH=high_potential_load,
                                  NOLOAD=modifier_potential_load)

potential_load_table <- potential_load_pool %>%
  group_by(biotype) %>%
  summarise(
    # GLOBAL
    n=n(),
    # LOW
    low.mean = mean(LOW),
    low.sd=sd(LOW),
    low.se =sd(LOW)/sqrt(sum(n())),
    # MODERATE
    moderate.mean = mean(MODERATE),
    moderate.sd=sd(MODERATE),
    moderate.se =sd(MODERATE)/sqrt(sum(n())),
    # HIGH
    high.mean = mean(HIGH),
    high.sd=sd(HIGH),
    high.se =sd(HIGH)/sqrt(sum(n())),
    # HIGH
    noload.mean = mean(NOLOAD),
    noload.sd=sd(NOLOAD),
    noload.se =sd(NOLOAD)/sqrt(sum(n()))
  )

potential_load4plot <- data.frame(mean=c(potential_load_table$low.mean, potential_load_table$moderate.mean, potential_load_table$high.mean, potential_load_table$noload.mean),
                                  sd=c(potential_load_table$low.sd, potential_load_table$moderate.sd, potential_load_table$high.sd, potential_load_table$noload.sd),
                                  se=c(potential_load_table$low.se, potential_load_table$moderate.se, potential_load_table$high.se, potential_load_table$noload.se),
                                  BIOTYPE=c(potential_load_table$biotype, potential_load_table$biotype, potential_load_table$biotype, potential_load_table$biotype),
                                  IMPACT=c(rep("A", 2), rep("B", 2), rep("C", 2), rep("D", 2))
                                  )

# Prepare the error bar to place in the barplot 
limits_potential_load <- aes(ymax = mean + sd,
                             ymin = mean - sd)


# Potential load plot
p10 <- ggplot(data = potential_load4plot, aes(x = BIOTYPE, y = mean, fill=IMPACT))
p10 <- p10 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
             geom_errorbar(limits_potential_load, position = position_dodge(0.9), width = 0.25) + 
             labs(x = NULL, y = "Potential Load") + 
             ylim(0, 0.75) + 
             scale_x_discrete(labels = c("Avirulent", "Virulent")) +
             scale_fill_manual(values = c("#fee090", "#fdae61", "#f46d43","#4575b4"),
                               name   = "Impact",
                               breaks = c("A", "B", "C", "D"),
                               labels = c("Low", "Moderate", "High", "Modifier")
                               ) + 
             theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                                axis.line = element_line(colour = "black"),
                                axis.text = element_text(size = 14),
                                axis.title=element_text(size=20),
                                axis.text.x = element_text(angle = 0))
p10

### CALCULATE REALIZED LOAD
# Realized Load: total number of homozygous mutations of impact class i in pools j / 2 x total number of sites of impact class i in pools j

# FIRST DEFINE A MAF THRESHOULD TO DEFINE IF THE ALTERNATIVE ALLELE IS PRESENT IN THE POOL
maf_homozygosity = 0.75

# Then assign 1 for MAF >= maf_potential, 0 otherwise
realized_load <- apply(dt.1.imputedRefMLFreq.dataFrame.codingsnps[,-c(1:4)], 2, function(x) ifelse(1-x >= maf_homozygosity, 1, 0))

# Combine with the SNP ANNOTATION IMPACT COLUMN
realized_load <- data.frame(IMPACT= summary_statistics_annotation_table.codingsnps$IMPACT,
                            realized_load)

low_realized_load <- apply(realized_load[, -c(1)], 2, function(x) sum(x[realized_load$IMPACT == "LOW"], na.rm = T)/sum(!is.na(x)))
moderate_realized_load <- apply(realized_load[, -c(1)], 2, function(x) sum(x[realized_load$IMPACT == "MODERATE"], na.rm = T)/sum(!is.na(x)))
high_realized_load <- apply(realized_load[, -c(1)], 2, function(x) sum(x[realized_load$IMPACT == "HIGH"], na.rm = T)/sum(!is.na(x)))
modifier_realized_load <- apply(realized_load[, -c(1)], 2, function(x) sum(x[realized_load$IMPACT == "MODIFIER"], na.rm = T)/sum(!is.na(x)))

realized_load_pool <- data.frame(pool=names(low_realized_load),
                                  biotype=do.call(rbind, strsplit(names(low_realized_load), split = '.', fixed = T))[,2],
                                  LOW=low_realized_load,
                                  MODERATE=moderate_realized_load,
                                  HIGH=high_realized_load,
                                  NOLOAD=modifier_realized_load)

realized_load_table <- realized_load_pool %>%
  group_by(biotype) %>%
  summarise(
    # GLOBAL
    n=n(),
    # LOW
    low.mean = mean(LOW),
    low.sd=sd(LOW),
    low.se =sd(LOW)/sqrt(sum(n())),
    # MODERATE
    moderate.mean = mean(MODERATE),
    moderate.sd=sd(MODERATE),
    moderate.se =sd(MODERATE)/sqrt(sum(n())),
    # HIGH
    high.mean = mean(HIGH),
    high.sd=sd(HIGH),
    high.se =sd(HIGH)/sqrt(sum(n())),
    # HIGH
    noload.mean = mean(NOLOAD),
    noload.sd=sd(NOLOAD),
    noload.se =sd(NOLOAD)/sqrt(sum(n()))
  )

realized_load4plot <- data.frame(mean=c(realized_load_table$low.mean, realized_load_table$moderate.mean, realized_load_table$high.mean, realized_load_table$noload.mean),
                                  sd=c(realized_load_table$low.sd, realized_load_table$moderate.sd, realized_load_table$high.sd, realized_load_table$noload.sd),
                                  se=c(realized_load_table$low.se, realized_load_table$moderate.se, realized_load_table$high.se, realized_load_table$noload.se),
                                  BIOTYPE=c(realized_load_table$biotype, realized_load_table$biotype, realized_load_table$biotype, realized_load_table$biotype),
                                  IMPACT=c(rep("A", 2), rep("B", 2), rep("C", 2), rep("D", 2))
)

# Prepare the error bar to place in the barplot 
limits_realized_load <- aes(ymax = mean + sd,
                             ymin = mean - sd)


# realized load plot
p11 <- ggplot(data = realized_load4plot, aes(x = BIOTYPE, y = mean, fill=IMPACT))
p11 <- p11 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
       geom_errorbar(limits_realized_load, position = position_dodge(0.9), width = 0.25) + 
       labs(x = NULL, y = "Realized Load") + 
       ylim(0, 0.75) + 
       scale_x_discrete(labels = c("Avirulent", "Virulent")) +
       scale_fill_manual(values = c("#fee090", "#fdae61", "#f46d43","#4575b4"),
                         name   = "Impact",
                         breaks = c("A", "B", "C", "D"),
                         labels = c("Low", "Moderate", "High", "Modifier")
       ) + 
       theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                          axis.line = element_line(colour = "black"),
                          axis.text = element_text(size = 14),
                          axis.title=element_text(size=20),
                          axis.text.x = element_text(angle = 0))
p11

## COMBINE PLOTS 10 AND 11
par(mfrow=c(1, 2))
load_plot <- ggarrange(p10,p11,
                       #labels = c("A", "B", "C", "D"),
                       font.label = list(size = 24, face = "bold"),
                       nrow=1, ncol=2, legend = "bottom", common.legend = T)

load_plot

ggsave("results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/genetic_load_biotypes.pdf",
       load_plot,
       device="pdf",
       width=15,
       height=8.5)

#save_checkpoint_8
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

###
###
### --- PROPORTION OF SYNONYMOUS AND NONSYNONYMOUS MUTATIONS ACROSS GENOME ---
###
### 

### DATA FOR THIS ANALYSIS WAS OBTAINED WITH:
### AglyPoolseq/utils/PNPS4VCF.py
### This script only considers missense_variants as nonsynonymous mutations and synonymous_variants as synonymous mutations,
### and counts only if the genotype in vcf file is 0/1 (polymorphic?)


pnps.data <- read.table(file = "results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.vcf.pnps.gz", 
                        header = T, na.strings = NA)

### FOR ALL SCAFFOLDS
# Remove genomewide estimates
scaffolds.pnps <- pnps.data %>%
  filter(Chrom !="genomewide") 

scaffolds.pnps <- scaffolds.pnps %>%
  mutate(scaffolds.pnps, 
         BIOTYPES=do.call(rbind, strsplit(scaffolds.pnps$POP, 
                                          split = '_', fixed = T))[,2])

## POOLS: Calculate the mean, sd, se pnps for each pool for all scaffolds
mean.scaffolds.pnps <- scaffolds.pnps %>%
  group_by(POP) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS, na.rm = T),
    pnps.sd=sd(pNpS, na.rm = T),
    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
  )

# Add a column with the biotype 
mean.scaffolds.pnps <- mutate(mean.scaffolds.pnps, 
                              BIOTYPES=do.call(rbind, strsplit(mean.scaffolds.pnps$POP,  split = '_', fixed = T))[,2],
                              POOLNAMES = poolnames)

# Prepare the error bar to place in the barplot 
limits_pnps <- aes(ymax = pnps.m + pnps.se,
                   ymin = pnps.m - pnps.se)


# pnps plot
p12 <- ggplot(data = mean.scaffolds.pnps, aes(x = POOLNAMES, y = pnps.m, fill=BIOTYPES))
p12 <- p12 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
              geom_errorbar(limits_pnps, position = position_dodge(0.9), width = 0.25) + 
              labs(x = NULL, y = expression(italic(p)[N]/italic(p)[S])) + 
              ylim(0, 2) + 
              #scale_x_discrete(labels = c("Avirulent", "Virulent")) +
              scale_fill_manual(values = c("#247F00","#AB1A53"),
                                name   = "Biotypes",
                                breaks = c("BIO1", "BIO4"),
                                labels = c(expression("Avirulent"), "Virulent")) + 
              theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                                 axis.line = element_line(colour = "black"),
                                 axis.text = element_text(size = 14),
                                 axis.title=element_text(size=20),
                                 axis.text.x = element_text(angle = 90))
p12

## BIOTYPES: Calculate the mean, sd, se pnps for each pool for BIOTYPES
mean.scaffolds.pnps.biotypes <- scaffolds.pnps %>%
  group_by(BIOTYPES) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS, na.rm = T),
    pnps.sd=sd(pNpS, na.rm = T),
    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
  )

# Prepare the error bar to place in the barplot 
limits_pnps.biotypes <- aes(ymax = pnps.m + pnps.se,
                            ymin = pnps.m - pnps.se)


# pnps plot
p13 <- ggplot(data = mean.scaffolds.pnps.biotypes, aes(x = BIOTYPES, y = pnps.m))
p13 <- p13 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
              geom_errorbar(limits_pnps.biotypes, position = position_dodge(0.9), width = 0.25) + 
              labs(x = NULL, y = expression(italic(p)[N]/italic(p)[S])) + 
              ylim(0, 2) + 
              scale_x_discrete(labels = c("Avirulent", "Virulent")) +
              #scale_fill_manual(values = c("#247F00","#AB1A53"),
              #                  name   = "Biotypes",
              #                  breaks = c("BIO1", "BIO4"),
              #                  labels = c(expression("Avirulent"), "Virulent")) + 
              theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                                 axis.line = element_line(colour = "black"),
                                 axis.text = element_text(size = 14),
                                 axis.title=element_text(size=20),
                                 axis.text.x = element_text(angle = 0))
p13

par(mfrow=c(3, 1))
pnps_plot <- ggarrange(ggarrange(p13,p12,
                                 labels = c("A", "B"),
                                 common.legend = T,
                                 legend = "bottom",
                                 font.label = list(size = 24, face = "bold"),
                                 nrow=2),
                       common.legend = T,
                       legend = "bottom"
)

pnps_plot

ggsave("results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/pnps_biotypes.pdf",
       pnps_plot,
       device="pdf",
       width=4.6,
       height=6.7)

#save_checkpoint_9
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/diversity.poolsnp.workspace27Ago21.RData")

### FOR SIGNIFICANT SCAFFOLDS
#significant_scaffolds <- read.table(file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst.bed", 
#                                    header = F, sep = "\t")
#
#significant_scaffolds <- unique(significant_scaffolds[,1])
#
## Remove genomewide estimates AND keep only significant scaffolds
#sign.scaffolds.pnps <- pnps.data %>%
#  filter(Chrom !="genomewide") 
#
#sign.scaffolds.pnps <- sign.scaffolds.pnps %>%
#  filter(Chrom %in% significant_scaffolds) 
#
#sign.scaffolds.pnps <- sign.scaffolds.pnps %>%
#  mutate(sign.scaffolds.pnps, 
#         BIOTYPES=do.call(rbind, strsplit(sign.scaffolds.pnps$POP, 
#                                          split = '_', fixed = T))[,2])
#
### POOLS: Calculate the mean, sd, se pnps for each pool for SIGNIFICANT scaffolds
#mean.sign.scaffolds.pnps <- sign.scaffolds.pnps %>%
#  group_by(POP) %>%
#  summarise(
#    n=n(),
#    pnps.m = mean(pNpS, na.rm = T),
#    pnps.sd=sd(pNpS, na.rm = T),
#    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
#  )
#
## Add a column with the biotype 
#mean.sign.scaffolds.pnps <- mutate(mean.sign.scaffolds.pnps, 
#                                   BIOTYPES=do.call(rbind, strsplit(mean.sign.scaffolds.pnps$POP,  split = '_', fixed = T))[,2],
#                                   POOLNAMES = poolnames)
#
## Prepare the error bar to place in the barplot 
#limits_sign.pnps <- aes(ymax = pnps.m + pnps.se,
#                       ymin = pnps.m - pnps.se)
#
#
## pnps plot
#p14 <- ggplot(data = mean.sign.scaffolds.pnps, aes(x = POOLNAMES, y = pnps.m, fill=BIOTYPES))
#p14 <- p14 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
#              geom_errorbar(limits_sign.pnps , position = position_dodge(0.9), width = 0.25) + 
#              labs(x = NULL, y = expression(italic(p)[N]/italic(p)[S])) + 
#              ylim(0, 2) + 
#              #scale_x_discrete(labels = c("Avirulent", "Virulent")) +
#              scale_fill_manual(values = c("#247F00","#AB1A53"),
#                                name   = "Biotypes",
#                                breaks = c("BIO1", "BIO4"),
#                                labels = c(expression("Avirulent"), "Virulent")) + 
#              theme_bw() + theme(panel.border = element_rect(colour = "black"), 
#                                 axis.line = element_line(colour = "black"),
#                                 axis.text = element_text(size = 14),
#                                 axis.title=element_text(size=20),
#                                 axis.text.x = element_text(angle = 90))
#p14
#
### BIOTYPES: Calculate the mean, sd, se pnps for each pool for BIOTYPES
#mean.sign.scaffolds.pnps.biotypes <- sign.scaffolds.pnps %>%
#  group_by(BIOTYPES) %>%
#  summarise(
#    n=n(),
#    pnps.m = mean(pNpS, na.rm = T),
#    pnps.sd=sd(pNpS, na.rm = T),
#    pnps.se =sd(pNpS, na.rm = T)/sqrt(sum(n()))
#  )
#
## Prepare the error bar to place in the barplot 
#limits_sign_pnps.biotypes <- aes(ymax = pnps.m + pnps.se,
#                                 ymin = pnps.m - pnps.se)
#
#
## pnps plot
#p15 <- ggplot(data = mean.sign.scaffolds.pnps.biotypes, aes(x = BIOTYPES, y = pnps.m))
#p15 <- p15 +  geom_bar(stat = "identity", position = position_dodge(0.9), orientation = 90) +
#              geom_errorbar(limits_sign_pnps.biotypes, position = position_dodge(0.9), width = 0.25) + 
#              labs(x = NULL, y = expression(italic(p)[N]/italic(p)[S])) + 
#              ylim(0, 2) + 
#              scale_x_discrete(labels = c("Avirulent", "Virulent")) +
#              #scale_fill_manual(values = c("#247F00","#AB1A53"),
#              #                  name   = "Biotypes",
#              #                  breaks = c("BIO1", "BIO4"),
#              #                  labels = c(expression("Avirulent"), "Virulent")) + 
#              theme_bw() + theme(panel.border = element_rect(colour = "black"), 
#                                 axis.line = element_line(colour = "black"),
#                                 axis.text = element_text(size = 14),
#                                 axis.title=element_text(size=20),
#                                 axis.text.x = element_text(angle = 0))
#p15



