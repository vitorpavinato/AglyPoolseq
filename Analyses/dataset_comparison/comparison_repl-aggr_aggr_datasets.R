##########################################################
#             COMPARISON BETWEEN TWO DATASETS            #
#     1) SNPs REPLICATED-AGGREGATED SNPs READs POOLS     #
#           2) SNPs AGGREGATED BAMs - POOLS              #
##########################################################

## WHY THEY PRODUCE DIFFERENT PCAs?
## POPULATION CLUSTER IN 1 IS STABLE? OR IT IS JUST RANDOM?
## IN WHAT WAY THE TWO WAYS OF PROCESSING THE DATA TO PRODUCE THE SNP DATA DIFFER?

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
library(ggplot2)
library(ggpubr)

# Import auxiliary R functions  
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

ALPHA=0.75

## POOL'S INFORMATION
# Pool's names
poolnames <- c("MN-B1.1", "MN-B4.1",
               "ND-B1.1", "ND-B4.1", 
               "NW-B1.1", "NW-B1.2", "NW-B4.1", "NW-B4.2", "NW-B4.3", 
               "PA-B1.1", "PA-B4.1", "PA-B4.2", 
               "WI-B1.1", "WI-B1.2", "WI-B4.1", "WI-B4.2", "WI-B4.3", 
               "WO-B1.1", "WO-B4.1", "WO-B4.2", "WO-B4.3")

# Pool's/biotypes colors - for PCA
biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53")


###
###
### ---- INVESTIGATE WHY PCAs FROM DIFFERENT DATASETS ARE DIFFERENT PART 1: EXPLORATORY ----
###
###

### DATASET 1: ALLELE FREQUENCY DATA OF TECHNICAL REPLICATES AGGREGATED FROM VCF READ COUNTS

## LOAD DATAFRAME WITH INFORMATION OF POLYMORPHIC SNPS USED FOR THE PCA
snpInfo.dt.repl <- read.table(file="results/technical_replicates/minmaxcov_3_99/snpinfo_replicates_aggregated_data.txt", header = T)

### DATASET 2: ALLELE FREQUENCY DATA OF TECHNICAL REPLICATES AGGREGATED FROM READ ALIGNMENT

## LOAD DATAFRAME WITH INFORMATION OF POLYMORPHIC SNPS USED FOR THE PCA
snpInfo.dt.aggr <- read.table(file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/snpinfo_aggregated_data_pools.txt", header = T)

## MEAN COVERAGE
# create a vector of breaks to use across datasets
br <- seq(0, 1300, 50)

# quick histogram to compare the breaks counts
mean.cov.hist.repl <- hist(snpInfo.dt.repl$MEAN_COV, br, plot = FALSE)
mean.cov.hist.aggr <- hist(snpInfo.dt.aggr$MEAN_COV, br, plot = FALSE)

mean.cov.data <- merge(x=cbind(mean.cov.hist.repl$breaks[-c(27)], mean.cov.hist.repl$counts/sum(mean.cov.hist.repl$counts)),
                       y=cbind(mean.cov.hist.aggr$breaks[-c(27)], mean.cov.hist.aggr$counts/sum(mean.cov.hist.aggr$counts)),
                       by=1)

plot(mean.cov.data$V2.x, mean.cov.data$V2.y, xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1)

# FANCY HISOTRAMS
mean.cov.data.hist <- data.frame(values=c(snpInfo.dt.repl$MEAN_COV, 
                                          snpInfo.dt.aggr$MEAN_COV),
                                 groups=c(rep("A", length(snpInfo.dt.repl$MEAN_COV)),
                                          rep("B", length(snpInfo.dt.aggr$MEAN_COV)))
                                 )

mean.cov.hists <- ggplot(mean.cov.data.hist, aes(x = values, fill = factor(groups))) + 
                  geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                  facet_wrap(vars("Mean COV")) +
                  scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                   alpha=ALPHA),
                                             alpha(c("#56B4E9"),
                                                   alpha=ALPHA)), 
                                    name="DATASET", 
                                    breaks= c("A", "B"),
                                    labels = c("Replicated-pools", "Aggregated-pools"))


## MIN COVERAGE
# get the frequency of each minimun coverage bin values
# compare % of min coverage bin values lower than the defined threshold for each dataset

# MIN COV OF DATASET 1 = 3
freq.mean.cov.bins.repl <- table(snpInfo.dt.repl$MIN_COV)/length(snpInfo.dt.repl$MIN_COV)*100
#           1            2            3  
# 0.019947637  0.713128039 51.100860242

# MIN COV OF DATASET 2 = 4
freq.mean.cov.bins.aggr <- table(snpInfo.dt.aggr$MIN_COV)/length(snpInfo.dt.aggr$MIN_COV)*100

# Compare the frequency of the first 10 classes in each dataset
plot(as.vector(freq.mean.cov.bins.repl[1:20]), as.vector(freq.mean.cov.bins.aggr[1:20]), pch=19,
     xlim=c(0,70), ylim=c(0,70), xlab = "Bins frequencies for the 'repl' dataset", ylab = "Bins frequencies for the 'aggr' dataset")
abline(a=0, b=1)


# FANCY HISOTRAMS
min.cov.data.hist <- data.frame(values=c(snpInfo.dt.repl$MIN_COV, 
                                          snpInfo.dt.aggr$MIN_COV),
                                 groups=c(rep("A", length(snpInfo.dt.repl$MIN_COV)),
                                          rep("B", length(snpInfo.dt.aggr$MIN_COV)))
)

min.cov.hists <- ggplot(min.cov.data.hist, aes(x = values, fill = factor(groups))) + 
                 geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                 facet_wrap(vars("Min COV")) +
                 scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                  alpha=ALPHA),
                                            alpha(c("#56B4E9"),
                                                  alpha=ALPHA)), 
                                   name="DATASET", 
                                   breaks= c("A", "B"),
                                   labels = c("Replicated-pools", "Aggregated-pools"))

## MAX COVERAGE

# FANCY HISOTRAMS
max.cov.data.hist <- data.frame(values=c(snpInfo.dt.repl$MAX_COV, 
                                         snpInfo.dt.aggr$MAX_COV),
                                groups=c(rep("A", length(snpInfo.dt.repl$MAX_COV)),
                                         rep("B", length(snpInfo.dt.aggr$MAX_COV)))
)

max.cov.hists <- ggplot(max.cov.data.hist, aes(x = values, fill = groups)) + 
                 facet_wrap(vars("Max COV")) +
                 geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                 scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                  alpha=ALPHA),
                                            alpha(c("#56B4E9"),
                                                  alpha=ALPHA)), 
                                   name="DATASET", 
                                   breaks= c("A", "B"),
                                   labels = c("Replicated-pools", "Aggregated-pools"))
               

## MEAN_REFFREQ
# FANCY HISOTRAMS
mean.reffreq.data.hist <- data.frame(values=c(snpInfo.dt.repl$MEAN_REFFREQ, 
                                             snpInfo.dt.aggr$MEAN_REFFREQ),
                                    groups=c(rep("A", length(snpInfo.dt.repl$MEAN_REFFREQ)),
                                             rep("B", length(snpInfo.dt.aggr$MEAN_REFFREQ)))
)

mean.reffreq.hists <- ggplot(mean.reffreq.data.hist, aes(x = values, fill = groups)) + 
                      geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                      facet_wrap(vars("Mean Ref Freq")) +
                      scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                       alpha=ALPHA),
                                                 alpha(c("#56B4E9"),
                                                       alpha=ALPHA)), 
                                        name="DATASET", 
                                        breaks= c("A", "B"),
                                        labels = c("Replicated-pools", "Aggregated-pools"))
                    
## MIN_REFFREQ
# FANCY HISOTRAMS
min.reffreq.data.hist <- data.frame(values=c(snpInfo.dt.repl$MIN_REFFREQ, 
                                             snpInfo.dt.aggr$MIN_REFFREQ),
                                    groups=c(rep("A", length(snpInfo.dt.repl$MIN_REFFREQ)),
                                             rep("B", length(snpInfo.dt.aggr$MIN_REFFREQ)))
)

min.reffreq.hists <- ggplot(min.reffreq.data.hist, aes(x = values, fill = groups)) + 
                     geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                     facet_wrap(vars("Min Ref Freq")) +
                     scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                      alpha=ALPHA),
                                                alpha(c("#56B4E9"),
                                                      alpha=ALPHA)), 
                                       name="DATASET", 
                                       breaks= c("A", "B"),
                                       labels = c("Replicated-pools", "Aggregated-pools"))

## MAX_REFFREQ
# FANCY HISOTRAMS
max.reffreq.data.hist <- data.frame(values=c(snpInfo.dt.repl$MAX_REFFREQ, 
                                             snpInfo.dt.aggr$MAX_REFFREQ),
                                    groups=c(rep("A", length(snpInfo.dt.repl$MAX_REFFREQ)),
                                             rep("B", length(snpInfo.dt.aggr$MAX_REFFREQ)))
)

max.reffreq.hists <- ggplot(max.reffreq.data.hist, aes(x = values, fill = groups)) + 
                     geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                     facet_wrap(vars("Max Ref Freq")) +
                     scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                      alpha=ALPHA),
                                                alpha(c("#56B4E9"),
                                                      alpha=ALPHA)), 
                                       name="DATASET", 
                                       breaks= c("A", "B"),
                                       labels = c("Replicated-pools", "Aggregated-pools"))
                   

## MISSING %
# FANCY HISOTRAMS
missing.data.hist <- data.frame(values=c(snpInfo.dt.repl$MISSING, 
                                         snpInfo.dt.aggr$MISSING),
                                groups=c(rep("A", length(snpInfo.dt.repl$MISSING)),
                                         rep("B", length(snpInfo.dt.aggr$MISSING)))
)

missing.hists <- ggplot(missing.data.hist, aes(x = values, fill = groups)) + 
                 geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
                 facet_wrap(vars("% Missing")) +
                 scale_fill_manual(values=c(alpha(c("#E69F00"),
                                                  alpha=ALPHA),
                                            alpha(c("#56B4E9"),
                                                  alpha=ALPHA)), 
                                   name="DATASET", 
                                   breaks= c("A", "B"),
                                   labels = c("Replicated-pools", "Aggregated-pools"))
               

par(mfrow=c(3, 1))
ggarrange(mean.cov.hists, min.cov.hists, max.cov.hists,
          mean.reffreq.hists, min.reffreq.hists, max.reffreq.hists,
          common.legend = T,
          legend = "bottom",
          font.label = list(size = 24, face = "bold"),
          nrow=2,
          ncol=3)


missing.hists

## Savepoint_1
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/dataset_comparison/repl-aggr-reads_aggr-bams_datasets.workspace28Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/dataset_comparison/repl-aggr-reads_aggr-bams_datasets.workspace28Jul21.RData")
############## END THIS PART

###
###
### ---- INVESTIGATE WHY PCAs FROM DIFFERENT DATASETS ARE DIFFERENT PART 2: COMMOM SNPS ----
###
###

## SUBSET THE DATASET 2 AND GET THE SAME SNPS IN DATASET 1

# Take SNPs in 2 similar in 1
snpInfo.dt.aggr.subset <- snpInfo.dt.aggr[paste0(snpInfo.dt.aggr$Chromosome,"_", snpInfo.dt.aggr$Position) %in% 
                                               paste0(snpInfo.dt.repl$Chromosome,"_",snpInfo.dt.repl$Position), ]

# Take SNPs in 1 similar to 2 subset
snpInfo.dt.repl.subset <- snpInfo.dt.repl[paste0(snpInfo.dt.repl$Chromosome,"_", snpInfo.dt.repl$Position) %in% 
                                            paste0(snpInfo.dt.aggr.subset$Chromosome,"_",snpInfo.dt.aggr.subset$Position), ]


# DOES THE TWO SUBSETS HAVE THE SAME SNPS
length(snpInfo.dt.aggr.subset$Chromosome) # 37132
length(snpInfo.dt.repl.subset$Chromosome) # 37132
sum(paste0(snpInfo.dt.aggr.subset$Chromosome,"_", snpInfo.dt.aggr.subset$Position) == 
      paste0(snpInfo.dt.repl.subset$Chromosome,"_", snpInfo.dt.repl.subset$Position)
    ) # 37132

## EXPORT THE LIST OF COMMON SNPS BETWEEN DATASETS
list_common_snps_names <- paste0(snpInfo.dt.repl.subset$Chromosome,"_", snpInfo.dt.repl.subset$Position)

# Save to a file only the combined chromosome and SNP position
write.table(list_common_snps_names,
            file="results/dataset_comparison/repl-aggr_aggr_common_snps_names.txt", sep = "\t",
            quote = F, row.names = F, col.names = F)
  
# MAKE SCATTERPLOTS FOR ALL VARIABLES AND LOOK FOR CORRELATIONS
plot(snpInfo.dt.aggr.subset[,5] ~ snpInfo.dt.repl.subset[,5], 
     xlim=c(0,1250), ylim=c(0,1250), xlab="MEAN COV REPL-AGGR", ylab="MEAN COV AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,5], snpInfo.dt.repl.subset[,5]) #0.9838574

plot(snpInfo.dt.aggr.subset[,6] ~ snpInfo.dt.repl.subset[,6], 
     xlim=c(0,290), ylim=c(0,290), xlab="MIN COV REPL-AGGR", ylab="MIN COV AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,6], snpInfo.dt.repl.subset[,6]) # 0.7824844

plot(snpInfo.dt.aggr.subset[,7] ~ snpInfo.dt.repl.subset[,7], 
     xlim=c(0,3500), ylim=c(0,3500), xlab="MAX COV REPL-AGGR", ylab="MAX COV AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,7], snpInfo.dt.repl.subset[,7]) # 0.9525403

plot(snpInfo.dt.aggr.subset[,8] ~ snpInfo.dt.repl.subset[,8], 
     xlim=c(0,1), ylim=c(0,1), xlab="MEAN_REFFREQ REPL-AGGR", ylab="MEAN_REFFREQ COV AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,8], snpInfo.dt.repl.subset[,8]) # 0.9942916

plot(snpInfo.dt.aggr.subset[,9] ~ snpInfo.dt.repl.subset[,9], 
     xlim=c(0,1), ylim=c(0,1), xlab="MIN_REFFREQ REPL-AGGR", ylab="MIN_REFFREQ COV AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,9], snpInfo.dt.repl.subset[,9]) # 0.9322479

plot(snpInfo.dt.aggr.subset[,10] ~ snpInfo.dt.repl.subset[,10], 
     xlim=c(0,1), ylim=c(0,1), xlab="MAX_REFFREQ REPL-AGGR", ylab="MAX_REFFREQ COV AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,10], snpInfo.dt.repl.subset[,10]) # 0.9016029

plot(snpInfo.dt.aggr.subset[,11] ~ snpInfo.dt.repl.subset[,11], 
     xlim=c(0,50), ylim=c(0,50), xlab="MISSING % REPL-AGGR", ylab="MISSING % AGGR")
abline(a=0, b=1)
cor(snpInfo.dt.aggr.subset[,11], snpInfo.dt.repl.subset[,11]) #0.3860252

## FINDING COMMON THRESHOLDS
# MEAN_COV
min(snpInfo.dt.aggr$MEAN_COV); max(snpInfo.dt.aggr$MEAN_COV)
min(snpInfo.dt.aggr.subset[,5]); max(snpInfo.dt.aggr.subset[,5])

min(snpInfo.dt.repl$MEAN_COV); max(snpInfo.dt.repl$MEAN_COV)
min(snpInfo.dt.repl.subset[,5]); max(snpInfo.dt.repl.subset[,5])
#min MEAN_COV = 7.666667

# MAX_COV
min(snpInfo.dt.aggr$MAX_COV); max(snpInfo.dt.aggr$MAX_COV)
min(snpInfo.dt.aggr.subset[,7]); max(snpInfo.dt.aggr.subset[,7])

min(snpInfo.dt.repl$MAX_COV); max(snpInfo.dt.repl$MAX_COV)
min(snpInfo.dt.repl.subset[,7]); max(snpInfo.dt.repl.subset[,7])
#min MAX_COV = 14

## Savepoint_2
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/dataset_comparison/repl-aggr-reads_aggr-bams_datasets.workspace28Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/dataset_comparison/repl-aggr-reads_aggr-bams_datasets.workspace28Jul21.RData")
############## END THIS PART

###
###
### ---- INVESTIGATE WHY PCAs FROM DIFFERENT DATASETS ARE DIFFERENT PART 3: COMPARATIVE PCAs ----
###
###

### USE THE COMMON SET OF SNPS TO PRODUCE THE PCA FROM THE RESPECTIVE ML ALLELE FREQUENCIES DATA

## UPLOAD THE ALLELE FREQUENCIES MATRIX

# DATASET 1: ALLELE FREQUENCY DATA OF TECHNICAL REPLICATES AGGREGATED FROM VCF READ COUNTS
dt.aggregated.imputedRefMLFreq.dataFrame <- read.table(file="results/technical_replicates/minmaxcov_3_99/imputedML_REFAlleleFrequencies_repl-aggregated_pools.txt", 
                                                       sep = "\t", header = T)
dim(dt.aggregated.imputedRefMLFreq.dataFrame)
# 40105; 25

# DATASET 2: ALLELE FREQUENCY DATA OF TECHNICAL REPLICATES AGGREGATED FROM READ ALIGNMENT
dt.1.imputedRefMLFreq.dataFrame <- read.table(file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_pools.txt", 
                                              sep = "\t", header = T)
dim(dt.1.imputedRefMLFreq.dataFrame)
# 262866; 25

## SUBSET EACH DATASET TO KEEP ONLY COMMON SNPS

# DATASET 1:
dt.aggregated.imputedRefMLFreq.subset <- dt.aggregated.imputedRefMLFreq.dataFrame[paste0(dt.aggregated.imputedRefMLFreq.dataFrame$Chromosome,"_", 
                                                                                         dt.aggregated.imputedRefMLFreq.dataFrame$Position) %in% 
                                                                                    paste0(snpInfo.dt.repl.subset$Chromosome,"_",
                                                                                           snpInfo.dt.repl.subset$Position),] 

# DATASET 2:
dt.1.imputedRefMLFreq.subset <-dt.1.imputedRefMLFreq.dataFrame[paste0(dt.1.imputedRefMLFreq.dataFrame$Chromosome, "_",
                                                                      dt.1.imputedRefMLFreq.dataFrame$Position) %in% 
                                                                 paste0(snpInfo.dt.repl.subset$Chromosome,"_",
                                                                        snpInfo.dt.repl.subset$Position), ]

## PCAs
# DATASET 1 - WITH ALL POOLS:
W_1.1 <- scale(t(dt.aggregated.imputedRefMLFreq.subset[,-c(1:4)]), scale=TRUE) #centering
W_1.1[1:10,1:10]
W_1.1[is.na(W_1.1)]<-0
cov.W_1.1<-cov(t(W_1.1))
eig.result_1.1<-eigen(cov.W_1.1)
eig.vec_1.1<-eig.result_1.1$vectors
lambda_1.1<-eig.result_1.1$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_1.1/sum(lambda_1.1),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.1,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_1.1/sum(lambda_1.1), col="red")

l1_1.1 <- 100*lambda_1.1[1]/sum(lambda_1.1) 
l2_1.1 <- 100*lambda_1.1[2]/sum(lambda_1.1)

# PCA plot 
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_1.1[,1],eig.vec_1.1[,2], col=biotype.col,
     xlim = c(-0.65, 0.35), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1_1.1,2), "%)"), ylab=paste0("PC1 (", round(l2_1.1,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_1.1[,1],eig.vec_1.1[,2], 
     poolnames, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.1)
abline(v=0,h=0,col="grey",lty=3)

# DATASET 1 - WITHOUT OUTLIERS POOLS:
W_1.2 <- scale(t(dt.aggregated.imputedRefMLFreq.subset[,-c(1:4, 9,25)]), scale=TRUE) #centering
W_1.2[1:10,1:10]
W_1.2[is.na(W_1.2)]<-0
cov.W_1.2<-cov(t(W_1.2))
eig.result_1.2<-eigen(cov.W_1.2)
eig.vec_1.2<-eig.result_1.2$vectors
lambda_1.2<-eig.result_1.2$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_1.2/sum(lambda_1.2),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_1.2/sum(lambda_1.2), col="red")

l1_1.2 <- 100*lambda_1.2[1]/sum(lambda_1.2) 
l2_1.2 <- 100*lambda_1.2[2]/sum(lambda_1.2)

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_1.2[,1],eig.vec_1.2[,2], col=biotype.col[-c(5,21)],
     xlim = c(-0.45, 0.75), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1_1.2,2), "%)"), ylab=paste0("PC1 (", round(l2_1.2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_1.2[,1],eig.vec_1.2[,2], 
     poolnames[-c(5,21)], pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# DATASET 2 - WITH ALL POOLS:
W_2.1 <- scale(t(dt.1.imputedRefMLFreq.subset[,-c(1:4)]), scale=TRUE) #centering
W_2.1[1:10,1:10]
W_2.1[is.na(W_2.1)]<-0
cov.W_2.1<-cov(t(W_2.1))
eig.result_2.1<-eigen(cov.W_2.1)
eig.vec_2.1<-eig.result_2.1$vectors
lambda_2.1<-eig.result_2.1$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_2.1/sum(lambda_2.1),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.1,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_2.1/sum(lambda_2.1), col="red")

l1_2.1 <- 100*lambda_2.1[1]/sum(lambda_2.1) 
l2_2.1 <- 100*lambda_2.1[2]/sum(lambda_2.1)

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_2.1[,1],eig.vec_2.1[,2], col=biotype.col,
     xlim = c(-0.55, 0.45), ylim = c(-0.65, 0.45),
     xlab=paste0("PC1 (", round(l1_2.1,2), "%)"), ylab=paste0("PC1 (", round(l2_2.1,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_2.1[,1],eig.vec_2.1[,2], 
     poolnames, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.1)
abline(v=0,h=0,col="grey",lty=3)

# DATASET 2 - WITHOUT OUTLIERS POOLS:
W_2.2 <- scale(t(dt.1.imputedRefMLFreq.subset[,-c(1:4, 9,25)]), scale=TRUE) #centering
W_2.2[1:10,1:10]
W_2.2[is.na(W_2.2)]<-0
cov.W_2.2<-cov(t(W_2.2))
eig.result_2.2<-eigen(cov.W_2.2)
eig.vec_2.2<-eig.result_2.2$vectors
lambda_2.2<-eig.result_2.2$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_2.2/sum(lambda_2.2),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_2.2/sum(lambda_2.2), col="red")

l1_2.2 <- 100*lambda_2.2[1]/sum(lambda_2.2) 
l2_2.2 <- 100*lambda_2.2[2]/sum(lambda_2.2)

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_2.2[,1],eig.vec_2.2[,2], col=biotype.col[-c(5,21)],
     xlim = c(-0.7, 0.45), ylim = c(-0.55, 0.65),
     xlab=paste0("PC1 (", round(l1_2.2,2), "%)"), ylab=paste0("PC1 (", round(l2_2.2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_2.2[,1],eig.vec_2.2[,2], 
     poolnames[-c(5,21)], pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

### COMBINED PLOTS
# PCA plot
pdf("results/dataset_comparison/pca_imputed_allelefreq_repl-aggr_aggr_common_snps.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)

par(mar=c(5,5,4,1)+.1, mfrow=c(2, 2))
# PCA plot DATASET 1 - WITH ALL POOLS
plot(eig.vec_1.1[,1],eig.vec_1.1[,2], col=biotype.col, main="repl-aggr all pools",
     xlim = c(-0.65, 0.35), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1_1.1,2), "%)"), ylab=paste0("PC1 (", round(l2_1.1,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_1.1[,1],eig.vec_1.1[,2], 
     poolnames, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.1)
abline(v=0,h=0,col="grey",lty=3)

# PCA plot DATASET 1 - WITHOUT OUTLIERS POOLS
plot(eig.vec_1.2[,1],eig.vec_1.2[,2], col=biotype.col[-c(5,21)], main="repl-aggr without high % missing pools",
     xlim = c(-0.45, 0.75), ylim = c(-0.55, 0.55),
     xlab=paste0("PC1 (", round(l1_1.2,2), "%)"), ylab=paste0("PC1 (", round(l2_1.2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_1.2[,1],eig.vec_1.2[,2], 
     poolnames[-c(5,21)], pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PCA plot DATASET 2 - WITH ALL POOLS:
plot(eig.vec_2.1[,1],eig.vec_2.1[,2], col=biotype.col, main="aggr all pools",
     xlim = c(-0.55, 0.45), ylim = c(-0.65, 0.45),
     xlab=paste0("PC1 (", round(l1_2.1,2), "%)"), ylab=paste0("PC1 (", round(l2_2.1,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_2.1[,1],eig.vec_2.1[,2], 
     poolnames, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.1)
abline(v=0,h=0,col="grey",lty=3)

# PCA plot DATASET 2 - WITHOUT OUTLIERS POOLS:
plot(eig.vec_2.2[,1],eig.vec_2.2[,2], col=biotype.col[-c(5,21)], main="aggr without high % missing pools",
     xlim = c(-0.7, 0.45), ylim = c(-0.55, 0.65),
     xlab=paste0("PC1 (", round(l1_2.2,2), "%)"), ylab=paste0("PC1 (", round(l2_2.2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_2.2[,1],eig.vec_2.2[,2], 
     poolnames[-c(5,21)], pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

dev.off()

## Savepoint_3
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/dataset_comparison/repl-aggr-reads_aggr-bams_datasets.workspace28Jul21.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/dataset_comparison/repl-aggr-reads_aggr-bams_datasets.workspace28Jul21.RData")
############## END THIS PART