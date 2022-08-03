##########################################################
#                   BAYPASS package                      # 
#              DEMOGRAPHY AND GENOME SCAN                #
##########################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

## Working directory
working_dir <- "/fs/project/PAS1554/aphidpool"

setwd(working_dir)
renv::restore()
#renv::snapshot()
## Remove last features
rm(list=ls())
ls()

library(corrplot)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ape)
library(geigen)
#library(epiR)
library(qvalue)

# Auxiliary functions from Baypass
source("/users/PAS1554/vitorpavinato/local/baypass/2.2.0/utils/baypass_utils.R")

# Auxiliary functions from aux_func.R
source("AglyPoolseq/analyses/aggregated_data/aux_func.R")

# Sample's names
poolnames <- c("MN-Av.1", "MN-V.1",
               "ND-Av.1", "ND-V.1", 
               "NW-Av.1", "NW-Av.2", "NW-V.1", "NW-V.2", "NW-V.3", 
               "PA-Av.1", "PA-V.1", "PA-V.2", 
               "WI-Av.1", "WI-Av.2", "WI-V.1", "WI-V.2", "WI-V.3", 
               "WO-Av.1", "WO-V.1", "WO-V.2", "WO-V.3")

biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53")

# Get the observed data parameters
n.snps = 262866
n.pools = 21

###
###
### ---- CHECK RUN CONVERGENCY WITH PILOT RUNS----
###
###Pilot runs were removed from the results folder

### The idea is to check the run convergency to later increase the size of the runs
### for the core model (and any other model?). The first check is going to be between
### the reference allele frequency and counts from the maximum likelihood imputation
### perfomed on the read counts loaded with poolfstat and the baypass posterior estimates
### of the reference allele frequency and counts. Because the baypass dataset was
### generated with a poolfstat function, the marker order should be the same.
### The second test is going to be with posterior estimates (frequency) between
### independe baypass runs.
#
### DATASETS
#
## FILE PATH
# Alelle frequencies
#ml.frequencies.file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_pools.txt"
#baypass.core.posteriors_0.file <- "/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core_summary_yij_pij.out"
#baypass.core.posteriors_1.file <- "/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core1_summary_yij_pij.out"
#baypass.core.posteriors_2.file <- "/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core2_summary_yij_pij.out"
#
# Marker list
#baypass.snps.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.snpdet"
#
# Omega matrices
#baypass.core.omega_0.file <- "/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core_mat_omega.out"
#baypass.core.omega_1.file <- "/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core1_mat_omega.out"
#baypass.core.omega_2.file <- "/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core2_mat_omega.out"
#
## IMPORT DATA
# Alelle frequencies
#ml.frequencies <- read.table(ml.frequencies.file, header = T)
#baypass.core.posteriors_0 <- read.table(baypass.core.posteriors_0.file, header = T)
#baypass.core.posteriors_1 <- read.table(baypass.core.posteriors_1.file, header = T)
#baypass.core.posteriors_2 <- read.table(baypass.core.posteriors_2.file, header = T)
#
# Marker list
#baypass.snps <- read.table(baypass.snps.file, header = F)
#colnames(baypass.snps) <- names(ml.frequencies[1:4])
#
# Omega matrices
#baypass.core.omega_0 <- as.matrix(read.table(baypass.core.omega_0.file, header = F))
#baypass.core.omega_1 <- as.matrix(read.table(baypass.core.omega_1.file, header = F))
#baypass.core.omega_2 <- as.matrix(read.table(baypass.core.omega_2.file, header = F))
#
#dimnames(baypass.core.omega_0)=list(poolnames, poolnames)
#dimnames(baypass.core.omega_1)=list(poolnames, poolnames)
#dimnames(baypass.core.omega_2)=list(poolnames, poolnames)
#
# Sanity check markers order between files
#sum(paste0(ml.frequencies$Chromosome,"_",ml.frequencies$Position) == 
#      paste0(baypass.snps$Chromosome, "_", baypass.snps$Position)) == n.snps
#
## Convert the table with BAYPASS posteriors to a matrix of nsnps x npops
#refFreq.matrix  <- matrix(baypass.core.posteriors_0$M_P, nrow=n.snps, ncol=n.pools, byrow=TRUE)
#refCount.matrix <- matrix(baypass.core.posteriors_0$M_Y, nrow=n.snps, ncol=n.pools, byrow=TRUE)
#
# check read count approximation with round()
#round(refFreq.matrix, 3)
#round(refCount.matrix)
#
### AMONG ML REFERENCE ALLELE AND POSTERIOR REFERENCE ALLELE
#
## FREQUENCY
# Create a list with Reference AF matrices - ML and Posteriors
#list.frequencies <- list(x=round(ml.frequencies[,-c(1:4)], 3), # ML Estimates on read count
#                         y=round(refFreq.matrix, 3)            # Baypass posterior estimates
#                         )                
#
# Calculate the correlation statistics for the allele frequencies
#af.corr <- calculate.multiple.correlations(list.frequencies)
#row.names(af.corr) <- poolnames
#print(af.corr)
#
## COUNT
# Create a list with Reference Read Count matrices - ML and Posteriors
#list.counts <- list(x=round(ml.frequencies[,-c(1:4)]*10), # ML Estimates on read count
#                    y=round(refCount.matrix)) # Posterior estimates
#
# Calculate the correlation statistics for the allele counts
#counts.corr <- calculate.multiple.correlations(list.counts)
#row.names(counts.corr) <- poolnames
#print(counts.corr)
#
# Plot Posterior Reference allale count against ML Reference allele count
#list.plots <- list()
#for (i in 1:n.pools)
#{
#  x=round(ml.frequencies[,-c(1:4)]*10)[,i]
#  y=round(refCount.matrix)[,i]
#  
#  df <- data.frame(x = x,y = y)
#  list.plots[[i]] <- ggplot(data = df,aes(x = x,y = y)) + stat_sum() + 
#    ggtitle (poolnames[i]) +
#    xlab("ML") +
#    ylab("Post") +
#    xlim(0,10) +
#    ylim(0,10) +
#    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") + 
#    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
#                       panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#                       legend.position = "none", axis.text = element_text(size = 12),
#                       axis.title=element_text(size=12))
#}
#
#p1 <- cowplot::plot_grid(list.plots[[1]],list.plots[[2]],list.plots[[3]],list.plots[[4]],
#                         list.plots[[5]],list.plots[[6]],list.plots[[7]],list.plots[[8]],
#                         list.plots[[9]],list.plots[[10]],list.plots[[11]],list.plots[[12]],
#                         list.plots[[13]],list.plots[[14]],list.plots[[15]],list.plots[[16]],
#                         list.plots[[17]],list.plots[[18]],list.plots[[19]],list.plots[[20]],
#                         list.plots[[21]], axis = "none",
#                         nrow = 6,
#                         ncol = 4,
#                         labels = NULL)
#
#save_plot(filename = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/refcount_correlation.pdf", 
#          p1, ncol = 2, base_asp = 1.1,
#          base_width=8.5, base_height=22)
#
### AMONG POSTERIOR ALLELE FREQUENCIES ESTIMATES OF 3 INDEPENDENT BAYPASS RUNS
#
## COUNT
#cor(baypass.core.posteriors_0$M_Y, baypass.core.posteriors_1$M_Y) 
#cor(baypass.core.posteriors_0$M_Y, baypass.core.posteriors_2$M_Y)
#cor(baypass.core.posteriors_1$M_Y, baypass.core.posteriors_2$M_Y)
#
## FREQUENCIES
#cor(baypass.core.posteriors_0$M_P, baypass.core.posteriors_1$M_P) 
#cor(baypass.core.posteriors_0$M_P, baypass.core.posteriors_2$M_P)
#cor(baypass.core.posteriors_1$M_P, baypass.core.posteriors_2$M_P)
#
### AMONG OMEGA MATRICES ESTIMATES OF 3 INDEPENDENT BAYPASS RUNS
#
## COMPARING SVD pFST and DISTANCE MATRIX
# Using SVD decomposition values
#omega_0.svd <- svd(baypass.core.omega_0)
#omega_1.svd <- svd(baypass.core.omega_1)
#omega_2.svd <- svd(baypass.core.omega_2)
#
#list.dist.svf <- list(x=omega_0.svd$u, y=omega_2.svd$u)
#svd.corr <- calculate.multiple.correlations(list.dist.svf)
#
# Using SVD plots
#par(mfrow=c(1, 3))
#omega_0.svd.plots <- plot.omega.mod(omega=baypass.core.omega_0, pop.names=poolnames, col=biotype.col)
#abline(v=0, h=0, lty=3,col = "gray60")
#omega_1.svd.plots <- plot.omega.mod(omega=baypass.core.omega_1, pop.names=poolnames, col=biotype.col)
#abline(v=0, h=0, lty=3,col = "gray60")
#omega_2.svd.plots <- plot.omega.mod(omega=baypass.core.omega_2, pop.names=poolnames, col=biotype.col)
#abline(v=0, h=0, lty=3,col = "gray60")

######START HEHE
###
###
### ---- CHECK RUN CONVERGENCY WITH INDEPEND RUNS ----
###
###

### After checking convergency between independet BAYPASS runs for the same mode but
### with different seeds (core model only). I decided the increase the run length 
### without comprimising running time. The same run length was used to produce results
### with the core mode (XtX statistics) and with the standard covariate model (STD) with 
### the Important Sampling (IS) algorithm. Because IS could be inconsistent, I run 3 different
### analyses with different seeds. The idea here is check the convergency between core, 
### STD models and the ML allele frequencies, and betweeen IS runs.

### DATASETS

## FILE PATH
# Alelle frequencies
ml.frequencies.file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_pools.txt"
baypass.core.posteriors_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core0/core0_summary_yij_pij.out"
baypass.stdis.posteriors_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis0/stdis0_summary_yij_pij.out"
baypass.stdis.posteriors_1.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis1/stdis1_summary_yij_pij.out"
baypass.stdis.posteriors_2.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis2/stdis2_summary_yij_pij.out"
baypass.aux.posteriors_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aux0/aux0_summary_yij_pij.out"

# Marker list
baypass.snps.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.snpdet"

# Omega matrices
baypass.core.omega_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core0/core0_mat_omega.out"
baypass.stdis.omega_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis0/stdis0_mat_omega.out"
baypass.stdis.omega_1.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis1/stdis1_mat_omega.out"
baypass.stdis.omega_2.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis2/stdis2_mat_omega.out"

## IMPORT DATA
# Alelle frequencies
ml.frequencies <- read.table(ml.frequencies.file, header = T)
baypass.core.posteriors_0 <- read.table(baypass.core.posteriors_0.file, header = T)
baypass.stdis.posteriors_0 <- read.table(baypass.stdis.posteriors_0.file, header = T)
baypass.stdis.posteriors_1 <- read.table(baypass.stdis.posteriors_1.file, header = T)
baypass.stdis.posteriors_2 <- read.table(baypass.stdis.posteriors_2.file, header = T)
baypass.aux.posteriors_0 <- read.table(baypass.aux.posteriors_0.file, header = T)

# Marker list
baypass.snps <- read.table(baypass.snps.file, header = F)
colnames(baypass.snps) <- names(ml.frequencies[1:4])

# Omega matrices
baypass.core.omega_0 <- as.matrix(read.table(baypass.core.omega_0.file, header = F))
baypass.stdis.omega_0 <- as.matrix(read.table(baypass.stdis.omega_0.file, header = F))
baypass.stdis.omega_1 <- as.matrix(read.table(baypass.stdis.omega_1.file, header = F))
baypass.stdis.omega_2 <- as.matrix(read.table(baypass.stdis.omega_2.file, header = F))

dimnames(baypass.core.omega_0)=list(poolnames, poolnames)
dimnames(baypass.stdis.omega_0)=list(poolnames, poolnames)
dimnames(baypass.stdis.omega_1)=list(poolnames, poolnames)
dimnames(baypass.stdis.omega_2)=list(poolnames, poolnames)

# Sanity check markers order between files
sum(paste0(ml.frequencies$Chromosome,"_",ml.frequencies$Position) == 
      paste0(baypass.snps$Chromosome, "_", baypass.snps$Position)) == n.snps

## Convert the table with BAYPASS posteriors to a matrix of nsnps x npops

# For the core model
core0.refFreq.matrix  <- matrix(baypass.core.posteriors_0$M_P, nrow=n.snps, ncol=n.pools, byrow=TRUE)
core0.refCount.matrix <- matrix(baypass.core.posteriors_0$M_Y, nrow=n.snps, ncol=n.pools, byrow=TRUE)

# For AUX model
aux0.refFreq.matrix  <- matrix(baypass.aux.posteriors_0$M_P, nrow=n.snps, ncol=n.pools, byrow=TRUE)
aux0.refCount.matrix <- matrix(baypass.aux.posteriors_0$M_Y, nrow=n.snps, ncol=n.pools, byrow=TRUE)

# For one STD model
stdis0.refFreq.matrix  <- matrix(baypass.stdis.posteriors_0$M_P, nrow=n.snps, ncol=n.pools, byrow=TRUE)
stdis0.refCount.matrix <- matrix(baypass.stdis.posteriors_0$M_Y, nrow=n.snps, ncol=n.pools, byrow=TRUE)

## check read count approximation with round()
#round(core0.refFreq.matrix, 3)
#round(core0.refCount.matrix)

## Savepoint_1
##-------------
#save.image("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
#load("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
############## END THIS PART

### CONVERGENCE: AMONG ML REFERENCE ALLELE AND POSTERIOR REFERENCE ALLELE

## FREQUENCY
# Create a list with Reference AF matrices - ML and Posteriors
list.frequencies <- list(x=round(ml.frequencies[,-c(1:4)], 3),       # ML Estimates on read count
                         y=round(aux0.refFreq.matrix, 3)) 
                         #y=round(core0.refFreq.matrix, 3))            # Baypass posterior estimates for core model

# Calculate the correlation statistics for the allele frequencies
af.corr <- calculate.multiple.correlations(list.frequencies)
row.names(af.corr) <- poolnames
print(af.corr)
rm(list.frequencies)
rm(af.corr)

## COUNT
# Create a list with Reference Read Count matrices - ML and Posteriors
list.counts <- list(x=round(ml.frequencies[,-c(1:4)]*10), # ML Estimates on read count
                    y=round(aux0.refCount.matrix))
                    #y=round(core0.refCount.matrix))       # Posterior estimates

# Calculate the correlation statistics for the allele counts
counts.corr <- calculate.multiple.correlations(list.counts)
row.names(counts.corr) <- poolnames
print(counts.corr)
rm(list.counts)
rm(counts.corr)

### AMONG POSTERIOR ALLELE FREQUENCIES ESTIMATES OF 3 INDEPENDENT BAYPASS RUNS

## COUNT
cor(baypass.core.posteriors_0$M_Y, baypass.aux.posteriors_0$M_Y)

#cor(baypass.core.posteriors_0$M_Y, baypass.stdis.posteriors_0$M_Y) 
#cor(baypass.core.posteriors_0$M_Y, baypass.stdis.posteriors_1$M_Y)
#cor(baypass.core.posteriors_0$M_Y, baypass.stdis.posteriors_2$M_Y)
#cor(baypass.stdis.posteriors_0$M_Y, baypass.stdis.posteriors_1$M_Y)
#cor(baypass.stdis.posteriors_0$M_Y, baypass.stdis.posteriors_2$M_Y)
#cor(baypass.stdis.posteriors_1$M_Y, baypass.stdis.posteriors_2$M_Y)

## FREQUENCIES
cor(baypass.core.posteriors_0$M_P, baypass.aux.posteriors_0$M_P)

#cor(baypass.core.posteriors_0$M_P, baypass.stdis.posteriors_0$M_P) 
#cor(baypass.core.posteriors_0$M_P, baypass.stdis.posteriors_1$M_P)
#cor(baypass.core.posteriors_0$M_P, baypass.stdis.posteriors_2$M_P)
#cor(baypass.stdis.posteriors_0$M_P, baypass.stdis.posteriors_1$M_P)
#cor(baypass.stdis.posteriors_0$M_P, baypass.stdis.posteriors_2$M_P)
#cor(baypass.stdis.posteriors_1$M_P, baypass.stdis.posteriors_2$M_P)

### AMONG OMEGA MATRICES ESTIMATES OF 3 INDEPENDENT BAYPASS RUNS

## COMPARING SVD pFST and DISTANCE MATRIX
# Using SVD decomposition values
omega_core0.svd <- svd(baypass.core.omega_0)
omega_stdis0.svd <- svd(baypass.stdis.omega_0)
omega_stdis1.svd <- svd(baypass.stdis.omega_1)
omega_stdis2.svd <- svd(baypass.stdis.omega_2)

list.dist.svf <- list(x=omega_core0.svd$u, y=omega_stdis0.svd$u)
svd.corr <- calculate.multiple.correlations(list.dist.svf)
print(svd.corr)
rm(list.dist.svf)
rm(svd.corr)

# Using SVD plots
# SAVE MANHATTAN PLOT TO A PDF FILE
pdf("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/omega_convergency.pdf", # File name
    width = 11, height = 8.50,                                                               # Width and height in inches
    bg = "white",                                                                            # Background color
    colormodel = "cmyk",                                                                     # Color model (cmyk is required for most publications)
)
par(mfrow=c(2, 2))
omega_core0.svd.plots <- plot.omega.mod(omega=baypass.core.omega_0, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")
omega_stdis0.svd.plots <- plot.omega.mod(omega=baypass.stdis.omega_0, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")
omega_stdis1.svd.plots <- plot.omega.mod(omega=baypass.stdis.omega_1, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")
omega_stdis2.svd.plots <- plot.omega.mod(omega=baypass.stdis.omega_2, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")
dev.off()

## Savepoint_2
##-------------
#save.image("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
#load("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
############## END THIS PART

######
### BAYPASS ANALYSES
######

###
###
### --- DEMOGRAPHIC INFERENCE ---
###
###
## With core model omega matrix (already calculated)

# Visualization of the OMEGA matrix
# As a PCoA with SVD decomposition
pdf("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/omega_core0_model_pca.pdf",  # File name
    width = 11, height = 8.50,                                                                    # Width and height in inches
    bg = "white",                                                                                 # Background color
    colormodel = "cmyk",                                                                          # Color model (cmyk is required for most publications)
)
omega_core0.svd.plots <- plot.omega.mod(omega=baypass.core.omega_0, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")
dev.off()

# as a correlation plot
omega_core0.cor_matrix = cov2cor(baypass.core.omega_0)
dimnames(omega_core0.cor_matrix) = list(poolnames, poolnames)

pdf("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/omega_core0_model_corrplot.pdf",  # File name
    width = 11, height = 8.50,                                                                             # Width and height in inches
    bg = "white",                                                                                          # Background color
    colormodel = "cmyk",                                                                                   # Color model (cmyk is required for most publications)
)
corrplot(omega_core0.cor_matrix, method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))
dev.off()

# as a heatmap and hierarchical clustering tree (using the average agglomeration method)
hclust.ave <- function(x) hclust(x, method="average")
pdf("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/omega_core0_model_heatmap.pdf",  # File name
    width = 11, height = 8.50,                                                                             # Width and height in inches
    bg = "white",                                                                                          # Background color
    colormodel = "cmyk",                                                                                   # Color model (cmyk is required for most publications)
)
heatmap(1-omega_core0.cor_matrix, hclustfun = hclust.ave,
        main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

# as a tree
pdf("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/omega_core0_model_tree.pdf",  # File name
    width = 11, height = 8.50,                                                                             # Width and height in inches
    bg = "white",                                                                                          # Background color
    colormodel = "cmyk",                                                                                   # Color model (cmyk is required for most publications)
)
plot(nj(as.dist(1-omega_core0.cor_matrix)))
dev.off()

## Savepoint_3
##-------------
#save.image("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
#load("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
############## END THIS PART

###
###
### --- SELECTION ---
###
###

##############
# CORE MODEL #
##############

## BAYPASS XtX ESTIMATES
# UploadXtX differentiation estimates (using the calibrated XtXst estimator)
baypass.core.xtx_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core0/core0_summary_pi_xtx.out"
baypass.core.xtx_0 <- read.table(baypass.core.xtx_0.file, header = T)

# check estimated XtXst distributions
hist(baypass.core.xtx_0$XtXst, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "XtX",
     main = "")
lines(density(baypass.core.xtx_0$XtXst), # density plot
      lwd = 2, # thickness of line
      col = "red")

# check behavior of the p-values associated to the estimated XtXst
hist(10**(-1*baypass.core.xtx_0$log10.1.pval.),freq=F,breaks=50)
abline(h=1) 
## CONCLUSION :: P-values don't behave well

# Manhattan plots of XtXst and p-values
layout(matrix(1:2,2,1))
plot(baypass.core.xtx_0$XtXst)
plot(baypass.core.xtx_0$log10.1.pval.,ylab="XtX P-value (-log10 scale)")
abline(h=3,lty=2, col="red") #0.001 p--value theshold

##
## -- CALIBRATE THE XtX POSTERIOR ESTIMATES WITH NORMALIZE TRANSFORMATION --
## 

# Posterior mean of XtX
xtx_core0 = baypass.core.xtx_0$M_XtX

hist(xtx_core0, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "XtX",
     main = "")
lines(density(xtx_core0), # density plot
      lwd = 2, # thickness of line
      col = "red")

# Posterior mean of the XtX NORMALIZED
xtx_core0.std = (xtx_core0 - mean(xtx_core0))/(sd(xtx_core0))

hist(xtx_core0.std, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "Standardized XtX",
     main = "")
lines(density(xtx_core0.std), # density plot
      lwd = 2, # thickness of line
      col = "red")

# Calculate empirical p-value
xtx_core0_empval = empPvals(stat = xtx_core0.std, stat0 = xtx_core0.std[which(xtx_core0.std <= 4)]) #With qvalue package
qqplot(x=xtx_core0_empval, y=10**(-1*baypass.core.xtx_0$log10.1.pval.)) 
abline(a = 0, b = 1, col="red")

# Threshold based on lowest empirical pvalue - 1%
thr_pval1ppm_core0 = min(xtx_core0[xtx_core0_empval < 0.01]) #XtX
thr_pval1ppm_core0.std = min(xtx_core0.std[xtx_core0_empval < 0.01]) #standardized XtX

# Threshold based on lowest empirical pvalue - 5%
thr_pval5ppm_core0 = min(xtx_core0[xtx_core0_empval < 0.05]) #XtX
thr_pval5ppm_core0.std = min(xtx_core0.std[xtx_core0_empval < 0.05]) #standardized XtX

# Threshold based on qvalue (FDR)
xtx_core0_empval.qvalue = qvalue(xtx_core0_empval)

# 1% FDR
thr_fdr1pcent_core0 = min(xtx_core0[xtx_core0_empval.qvalue$qvalue < 0.01])#XtX
thr_fdr1pcent_core0.std = min(xtx_core0.std[xtx_core0_empval.qvalue$qvalue < 0.01])

# 5% FDR
thr_fdr5pcent_core0 = min(xtx_core0[xtx_core0_empval.qvalue$qvalue < 0.05])#XtX
thr_fdr5pcent_core0.std = min(xtx_core0.std[xtx_core0_empval.qvalue$qvalue < 0.05]) #standardized XtX

# 10% FDR
thr_fdr10pcent_core0 = min(xtx_core0[xtx_core0_empval.qvalue$qvalue < 0.1])#XtX
thr_fdr10pcent_core0.std = min(xtx_core0.std[xtx_core0_empval.qvalue$qvalue < 0.1]) #standardized XtX

# compute the 1% quantile threshold from the observed distribution
thr_quantile99_core0 = quantile(xtx_core0, probs=c(0.99))
thr_quantile01_core0 = quantile(xtx_core0, probs=c(0.01))

thr_quantile99_core0.std = quantile(xtx_core0.std, probs=c(0.99))
thr_quantile01_core0.std = quantile(xtx_core0.std, probs=c(0.01))

##
## -- CALIBRATE THE XtX STATISTICS WITH PODs --
##

# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution
baypass.core.beta_0.file <- "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core0/core0_summary_beta_params.out"
baypass.core.beta_0 <- read.table(baypass.core.beta_0.file, header = T)$Mean

# upload the original BAYPASS data to obtain the total allele count
baypass_genobaypass.file <- "/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.genobaypass"
baypass_genobaypass <- geno2YN(baypass_genobaypass.file)

# Create the POD
#setwd("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/")
sim.core.data <- simulate.baypass(omega.mat=baypass.core.omega_0, 
                                  nsnp=100000, sample.size = 10, 
                                  coverage = baypass_genobaypass$NN, 
                                  beta.pi=baypass.core.beta_0, pi.maf=0, suffix="core_pods")


# RUN BAYPASS WITH POD DATA
#setwd("/fs/scratch/PAS1715/aphidpool")
suffix_="core_pods"
g_file <- paste0("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/", "Gpool.", suffix_)
poolsize_file <- paste0("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/","poolsize.", suffix_)
outprefix <- paste0("core.pods")
logfile <- paste0("core.pod.log")
outdir <- paste0("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/")
wd <- paste0("/fs/scratch/PAS1715/aphidpool")

core_pods_baypass <- make_baypassrun_slurm_pods(n_hours=10, g_file = g_file, 
                                                poolsize_file = poolsize_file, 
                                                outprefix = outprefix, 
                                                logfile = logfile, 
                                                outdir = outdir)

## TO RUN BAYPASS IN PODS, RUN IN THE TERMINAL:
# wd=/fs/scratch/PAS1715/aphidpool # to run on scratch
# sbatch ${wd}/run_baypass_pods.sh

#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################

# get estimate of omega from the POD analysis
pod.core.omega_0 = as.matrix(read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/pods/corePods/core.pods_mat_omega.out"))
plot(pod.core.omega_0, baypass.core.omega_0) ; abline(a=0,b=1)
fmd.dist(pod.core.omega_0, baypass.core.omega_0)

layout(matrix(1:2,1,2))
omega_core0.svd.plots <- plot.omega.mod(omega=baypass.core.omega_0, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")

omega_core0_pod.svd.plots <- plot.omega.mod(omega=pod.core.omega_0, pop.names=poolnames, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")

# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/pods/corePods/core.pods_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef, baypass.core.beta_0)  ; abline(a=0,b=1)
cor(pod.pi.beta.coef, baypass.core.beta_0)

#######################################################
# XtX calibration
#######################################################
#get the pod XtX
xtx_core0_pod = read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/pods/corePods/core.pods_summary_pi_xtx.out",h=T)$M_XtX

# compute the quantile 1% threshold
thr_quantile99_core0_pod = quantile(xtx_core0_pod, probs=c(0.99))
thr_quantile01_core0_pod = quantile(xtx_core0_pod, probs=c(0.01))

# Compute the p-values of XtX calibrated with PODS XtX (From Gautier 2017 presentation)
hist(xtx_core0, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "Standardized XtX",
     main = "")
lines(density(xtx_core0_pod), # density plot
      lwd = 2, # thickness of line
      col = "red")

# Posterior mean of the XtX NORMALIZED
xtx_core0_pod.std = (xtx_core0_pod - mean(xtx_core0_pod))/(sd(xtx_core0_pod))

hist(xtx_core0.std, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "Standardized XtX",
     main = "")
lines(density(xtx_core0_pod.std), # density plot
      lwd = 2, # thickness of line
      col = "red")

# compute the quantile 1% threshold
thr_quantile99_core0_pod.std = quantile(xtx_core0_pod.std, probs=c(0.99))
thr_quantile01_core0_pod.std = quantile(xtx_core0_pod.std, probs=c(0.01))

# Calculate empirical p-value
xtx_core0_pod_empval = empPvals(stat = xtx_core0.std, stat0 = xtx_core0_pod.std[which(xtx_core0_pod.std <= 4)]) #With qvalue package
qqplot(x=xtx_core0_pod_empval, y=10**(-1*baypass.core.xtx_0$log10.1.pval.)) 
abline(a = 0, b = 1, col="red")

xtx_core0_pod_empval.qvalue = qvalue(xtx_core0_pod_empval)

# Threshold based on qvalue (FDR)
# FDR 1%
thr_fdr1pcent_core0_pod = min(xtx_core0[xtx_core0_pod_empval.qvalue$qvalue < 0.01])# XtX
thr_fdr1pcent_core0_pod.std = min(xtx_core0.std[xtx_core0_empval.qvalue$qvalue < 0.01]) #standardized XtX

# FDR 5%
thr_fdr5pcent_core0_pod = min(xtx_core0[xtx_core0_pod_empval.qvalue$qvalue < 0.05])# XtX
thr_fdr5pcent_core0_pod.std = min(xtx_core0.std[xtx_core0_empval.qvalue$qvalue < 0.05]) #standardized XtX

# FDR 10%
thr_fdr10pcent_core0_pod = min(xtx_core0[xtx_core0_pod_empval.qvalue$qvalue < 0.1])# XtX
thr_fdr10pcent_core0_pod.std = min(xtx_core0.std[xtx_core0_empval.qvalue$qvalue < 0.1]) #standardized XtX

### Make plot of XtX and XtX normalized
### Add all thresholds

## Plot XtX with empirically and pod calculate threshold
## lines are empirically calculated
## dashed lines are calculated with pods
#layout(matrix(1:2,2,1))
plot(xtx_core0)
#abline(h=thr_pval1ppm_core0,  lty=1, col="#6fb157") #green  **
#abline(h=thr_pval5ppm_core0,  lty=1, col="#736cb9") #blue   *
#abline(h=thr_fdr1pcent_core0, lty=1, col="#bd8f3a") #orange ***
#abline(h=thr_fdr5pcent_core0, lty=1, col="#ac568b") #pink   ***
#abline(h=thr_fdr10pcent_core0,lty=1, col="#4cb6ad") #blue   ***
#abline(h=thr_quantile99_core0,lty=1, col="#b3534b") #red    **
#abline(h=thr_quantile01_core0,lty=1, col="#77874f") #green  balancing selection?

#abline(h=thr_fdr1pcent_core0_pod, lty=2, col="#bd8f3a")  #orange ***
#abline(h=thr_fdr5pcent_core0_pod, lty=2, col="#ac568b")  #pink   ***
#abline(h=thr_fdr10pcent_core0_pod,lty=2, col="#4cb6ad")  #blue   ***
abline(h=thr_quantile99_core0_pod, lty=2, col="#b3534b")  #red    **
abline(h=thr_quantile01_core0_pod,lty=2, col="#77874f")   #green  balancing selection?

## Plot normalized XtX with empirically and pod calculate threshold
## lines are empirically calculated
## dashed lines are calculated with pods
plot(xtx_core0.std)
#abline(h=thr_pval1ppm_core0.std,  lty=1, col="#6fb157") #green  **
#abline(h=thr_pval5ppm_core0.std,  lty=1, col="#736cb9") #blue   *
#abline(h=thr_fdr1pcent_core0.std, lty=1, col="#bd8f3a") #orange ***
#abline(h=thr_fdr5pcent_core0.std, lty=1, col="#ac568b") #pink   ***
#abline(h=thr_fdr10pcent_core0.std,lty=1, col="#4cb6ad") #blue   ***
#abline(h=thr_quantile99_core0.std,lty=1, col="#b3534b") #red    **
#abline(h=thr_quantile01_core0.std,lty=1, col="#77874f") #green  balancing selection?

#abline(h=thr_fdr1pcent_core0_pod.std, lty=2, col="#bd8f3a")  #orange ***
#abline(h=thr_fdr5pcent_core0_pod.std, lty=2, col="#ac568b")  #pink   ***
#abline(h=thr_fdr10pcent_core0_pod.std,lty=2, col="#4cb6ad")  #blue   ***
abline(h=thr_quantile99_core0_pod.std, lty=2, col="#b3534b")  #red    **
abline(h=thr_quantile01_core0_pod.std,lty=2, col="#77874f")   #green  balancing selection?


## Identify significant high values of XtX statistics 
## Use the 99% quantile from the PODs distribution
dim(baypass.snps[which(xtx_core0 > thr_quantile99_core0_pod), ]) # 4645 positive selection SNPs
dim(baypass.snps[which(xtx_core0 < thr_quantile01_core0_pod ), ]) # 7463 balanced selection?

## Export candidate outliers: core model XtX 
write.table(baypass.snps[which(xtx_core0 > thr_quantile99_core0_pod), c(1:2)],
            file = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/signf_snps_core0_podsQuant99.txt",
            col.names = F, row.names = F, quote = F)


## Savepoint_4
##-------------
#save.image("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
#load("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
############## END THIS PART

################
#   AUX MODEL  #
################
baypass.aux.betai_reg_0.file <- "results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aux0/aux0_summary_betai.out"
baypass.aux.betai_reg_0 <- read.table(baypass.aux.betai_reg_0.file, header = T)

baypass.aux.xtx_0.file <- "results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aux0/aux0_summary_pi_xtx.out"
baypass.aux.xtx_0 <- read.table(baypass.aux.xtx_0.file, header = T)

## Check correlation between core0 XtX and aux0 XtX
cor(baypass.aux.xtx_0$M_XtX, xtx_core0)
plot(baypass.aux.xtx_0$M_XtX, xtx_core0)

## Take a preliminary view on the AUX results
graphics.off()
layout(matrix(1:4,4,1))
plot(baypass.aux.betai_reg_0$BF.dB.,xlab="SNP",ylab="BFmc (in dB)") 
plot(baypass.aux.betai_reg_0$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient")) 
plot(baypass.aux.xtx_0$M_XtX,xlab="SNP",ylab="XtX corrected for SMS")
plot(xtx_core0,xlab="SNP",ylab="core XtX")

## Preliminar plot only with BFmc Jeffrey's rule
graphics.off()
plot(baypass.aux.betai_reg_0$BF.dB.,xlab="SNP",ylab="BFmc (in dB)")
abline(h=c(15,20),  lty=2, col=c("green","red"))

thr_BFmc_decisive_aux0 = 20
thr_BFmc_strongEvi_aux0 = 15

##
## -- CALIBRATE THE BFmc STATISTICS WITH PODs --
##

## With the same PODs generated for the core model, run baypass AUX this time with a covariate file

#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################

# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod.pi.beta.coef_aux=read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/pods/auxPods/core_aux.pods_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef_aux, baypass.core.beta_0)  ; abline(a=0,b=1)
cor(pod.pi.beta.coef_aux, baypass.core.beta_0)

#######################################################
# BFmc calibration
#######################################################

# the estimated BFmc
bfmc_aux0 = baypass.aux.betai_reg_0$BF.dB.

#get the pod BFmc
bfmc_aux0_pod = read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/pods/auxPods/core_aux.pods_summary_betai.out",h=T)$BF.dB.

# compute the quantile 1% threshold
thr_quantile99_aux0_pod = quantile(bfmc_aux0_pod, probs=c(0.99))

# Compute the p-values of XtX calibrated with PODS XtX (From Gautier 2017 presentation)
hist(bfmc_aux0, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "BFmc",
     main = "")
lines(density(bfmc_aux0_pod), # density plot
      lwd = 2, # thickness of line
      col = "red")

## Calculate empirical p-value and qvalue
bfmc_aux0_pod_empval = empPvals(stat = bfmc_aux0, stat0 = bfmc_aux0_pod[which(bfmc_aux0_pod <= 30)]) #With qvalue package

# Threshold based on lowest empirical pvalue - 1%
thr_pval1ppm_aux0_pod = min(bfmc_aux0[bfmc_aux0_pod_empval < 0.01]) 
thr_pval5ppm_aux0_pod = min(bfmc_aux0[bfmc_aux0_pod_empval < 0.05])

## Threshold based on qvalue (FDR)
bfmc_aux0_pod_empval.qvalue = qvalue(bfmc_aux0_pod_empval)

# FDR 1%
thr_fdr1pcent_aux0_pod = min(bfmc_aux0[bfmc_aux0_pod_empval.qvalue$qvalue < 0.01])
thr_fdr5pcent_aux0_pod = min(bfmc_aux0[bfmc_aux0_pod_empval.qvalue$qvalue < 0.05])
thr_fdr10pcent_aux0_pod = min(bfmc_aux0[bfmc_aux0_pod_empval.qvalue$qvalue < 0.1])

## Plot XtX with empirically and pod calculate threshold
## lines are empirically calculated
## dashed lines are calculated with pods
plot(bfmc_aux0)
abline(h=thr_BFmc_decisive_aux0,  lty=1, col="red") #**
abline(h=thr_BFmc_strongEvi_aux0,  lty=1, col="green")# *
abline(h=thr_quantile99_aux0_pod,  lty=2, col="red") #red  ***
#abline(h=thr_pval1ppm_aux0_pod,  lty=2, col="#736cb9") #blue   **
#abline(h=thr_pval5ppm_aux0_pod, lty=2, col="#bd8f3a") #orange *
#abline(h=thr_fdr1pcent_aux0_pod , lty=2, col="#ac568b") #pink   ****
#abline(h=thr_fdr5pcent_aux0_pod ,lty=2, col="#4cb6ad") #blue   ****
#abline(h=thr_fdr10pcent_aux0_pod ,lty=2, col="#b3534b") #red    ****

## Identify significant high values of XtX statistics 
## Use the 99% quantile from the PODs distribution
dim(baypass.snps[which(bfmc_aux0 > thr_quantile99_aux0_pod), ]) # 2054 positive selected SNPs associated with covariates

## Export candidate outliers: aux model BFmc 
# quantile 1%
write.table(baypass.snps[which(bfmc_aux0 > thr_quantile99_aux0_pod), c(1:2)],
            file = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/signf_snps_aux0_podsQuant99.txt",
            col.names = F, row.names = F, quote = F)

# BFmc > 20
write.table(baypass.snps[which(bfmc_aux0 > thr_BFmc_decisive_aux0), c(1:2)],
            file = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/signf_snps_aux0_BFmc20.txt",
            col.names = F, row.names = F, quote = F)

### Combine evidence an plot significant XtX, BFmc and multi-locus FST
final.table <- data.frame(baypass.snps,
                          XtX=xtx_core0,
                          BFmc=bfmc_aux0,
                          Beta_is=baypass.aux.betai_reg_0$M_Beta
                          )
colnames(final.table) <- c("CHROM", "POS", "Ref", "Alt", "XtX", "BFmc", "Beta_is")

final.table[,"dots_core_aux"] <- ifelse(final.table$XtX > thr_quantile99_core0_pod & final.table$BFmc > thr_quantile99_aux0_pod, 1, 20)
final.table[,"colors_core_aux"] <- ifelse(final.table$XtX > thr_quantile99_core0_pod & final.table$BFmc > thr_quantile99_aux0_pod, "#fca311", "#e5e5e5")

## BAYPASS XtX AND BFmc
pdf("results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/XtX_BFmc_sign.pdf",  # File name
    width = 11, height = 8.50,                                                                             # Width and height in inches
    bg = "white",                                                                                          # Background color
    colormodel = "cmyk",                                                                                   # Color model (cmyk is required for most publications)
)
plot(
  final.table$XtX,
  final.table$BFmc,
  ylab="BFmc (in dB)",xlab="XtX", pch = final.table$dots_core_aux, col= final.table$colors_core_aux)
abline(v=thr_quantile99_core0_pod,lty=2, col="black") 
abline(h=20,lty=2, col="black") #BF threshold for decisive evidence (according to Jeffreys’ rule)
abline(h=thr_quantile99_aux0_pod,lty=2, col="red") 
dev.off()

##  BAYPASS XtX AND BFmc AND FST
lociids_signf_WC_beta_window <-scan(file = "results/aggregated_data/minmaxcov_4_99/signf_multilocusfst_poolsnp/lociIDs_signf_WC_beta_window.txt",what = "character")
lociids_signf_WC_beta_window <- unique(lociids_signf_WC_beta_window)
length(lociids_signf_WC_beta_window)

final.table[,"colores_fst"] <- NA
final.table[paste0(final.table$CHROM,"_",final.table$POS) %in% lociids_signf_WC_beta_window, "colores_fst"] <- "#14213d"

pdf("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/XtX_BFmc_WC_beta_sign.pdf",  # File name
    width = 11, height = 8.50,                                                                             # Width and height in inches
    bg = "white",                                                                                          # Background color
    colormodel = "cmyk",                                                                                   # Color model (cmyk is required for most publications)
)
plot(
  final.table$XtX,
  final.table$BFmc,
  ylab="BFmc (in dB)",xlab="XtX", pch = final.table$dots_core_aux, col= final.table$colors_core_aux)
points(final.table$XtX, final.table$BFmc, col=final.table$colores_fst)
abline(v=thr_quantile99_core0_pod,lty=2, col="black") 
abline(h=20,lty=2, col="black") #BF threshold for decisive evidence (according to Jeffreys’ rule)
abline(h=thr_quantile99_aux0_pod,lty=2, col="#e63946")
text(x=18.8, y=c(20, thr_quantile99_aux0_pod)-0.1,
     c("decisive BF",
       "quantile 99%"),
     cex=0.75, pos=3,col=c("black", "#e63946"))

dev.off()

## Savepoint_5
##-------------
#save.image("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
#load("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/analysis/baypass.poolsnp.workspace01Jan22.RData")
############## END THIS PART

################
# STD IS MODEL #
################
baypass.stdis.betai_reg_0.file <- "results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis0/stdis0_summary_betai_reg.out"
baypass.stdis.betai_reg_1.file <- "results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis1/stdis1_summary_betai_reg.out"
baypass.stdis.betai_reg_2.file <- "results/analysis_osc//aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis2/stdis2_summary_betai_reg.out"

baypass.stdis.betai_reg_0 <- read.table(baypass.stdis.betai_reg_0.file, header = T)
baypass.stdis.betai_reg_1 <- read.table(baypass.stdis.betai_reg_1.file, header = T)
baypass.stdis.betai_reg_2 <- read.table(baypass.stdis.betai_reg_2.file, header = T)

plot(baypass.stdis.betai_reg_0$BF.dB., baypass.stdis.betai_reg_2$BF.dB.)
## high correlation among Beta and eBPis among different runs
## BFis are not correlated
## Use only Betas and eBPis of one run, forget about BFis
## Alleles are sample randomly?

baypass.stdis.xtx_0.file <- "results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis0/stdis0_summary_pi_xtx.out"
baypass.stdis.xtx_0 <- read.table(baypass.stdis.xtx_0.file, header = T)

baypass.stdis.xtx_1.file <- "results/analysis_osc/aggregated_data/minmaxcov_4_99/baypass_poolsnp/stdis1/stdis1_summary_pi_xtx.out"
baypass.stdis.xtx_1 <- read.table(baypass.stdis.xtx_1.file, header = T)

plot(baypass.core.xtx_0$M_XtX, baypass.stdis.xtx_0$M_XtX)
plot(baypass.core.xtx_0$M_XtX, baypass.stdis.xtx_1$M_XtX)

stdi_BFis_table <- cbind(baypass.stdis.betai_reg_0$BF.dB., baypass.stdis.betai_reg_1$BF.dB., baypass.stdis.betai_reg_2$BF.dB.)
stdi_BFis <- apply(stdi_BFis_table, 1, FUN = function(x) median(x, na.rm = T))

plot(density(stdi_BFis_table[,1]), lwd = 2, col="red")
lines(density(stdi_BFis_table[,2]), # density plot
      lwd = 2, # thickness of line
      col = "blue")
lines(density(stdi_BFis_table[,3]), # density plot
      lwd = 2, # thickness of line
      col = "green")

plot(stdi_BFis_table[,1], stdi_BFis_table[,2])

plot(baypass.stdis.betai_reg_0$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
abline(h=20,lty=2, col="red")

baypass.snps[which(xtx_core0 < 18 ), ]

graphics.off()
layout(matrix(1:3,3,1))
plot(baypass.stdis.betai_reg_0$BF.dB.,xlab="SNP",ylab="BFis (in dB)") 
plot(baypass.stdis.betai_reg_0$eBPis,xlab="SNP",ylab="eBPis") 
plot(baypass.stdis.betai_reg_0$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))

plot(
     baypass.core.xtx_0$XtXst,
     baypass.stdis.xtx_0$XtXst,
     ylab="BF",xlab="XtX p-value (-log10 scale)")
abline(h=20,lty=2) #BF threshold for decisive evidence (according to Jeffreys’ rule)
abline(v=15,lty=2) 


final.table <- data.frame(baypass.snps,
                          XtX=baypass.core.xtx_0$XtXst,
                          BF=baypass.stdis.betai_reg_0$BF.dB.,
                          eBPi=baypass.stdis.betai_reg_0$eBPis,
                          Beta_is=baypass.stdis.betai_reg_0$Beta_is)

colnames(final.table) <- c("CHROM", "POS", "Ref", "Alt", "XtX", "BF", "eBPis", "Beta_is")
final.table[which(final.table$XtX > 100 & final.table$CHROM == "scaffold_255"), ]
