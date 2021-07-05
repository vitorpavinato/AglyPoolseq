##########################################################
#                   BAYPASS package                      # 
#              DEMOGRAPHY AND GENOME SCAN                #
##########################################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES

setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()
## Remove last features
rm(list=ls())
ls()

library(corrplot)
library(ggplot2)
library(cowplot)
library(qvalue)

# Install Bioconductor qvalue package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("qvalue")

# Auxiliary functions from Baypass
source("/users/PAS1554/vitorpavinato/local/baypass/2.2.0/utils/baypass_utils.R")

# Auxiliary functions from aux_func.R
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

## POOL NAMES - 19
pool.names <- c("MN-B1.1", "MN-B4.1",
                "ND-B1.1", "ND-B4.1", 
                "NW-B1.2", "NW-B4.1", "NW-B4.2", "NW-B4.3", 
                "PA-B1.1", "PA-B4.1", "PA-B4.2", 
                "WI-B1.1", "WI-B1.2", "WI-B4.1", "WI-B4.2", "WI-B4.3", 
                "WO-B1.1", "WO-B4.1", "WO-B4.2")

biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53")


######
### CHECK RUN CONVERGENCY
######

### WITH ALLELE FREQUENCY ESTIMATES
## ML -based Estimate vs BAYPASS' Posterior Estimates

# Upload Maximum Likelihood Allele Frequency estimates from POOLFSTAT-POOLSNP Analyses
load("results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/reduced.dt.1.imputedRefMLFreq.filtered03.RData")

# check read count approximation with round()
#round(head(reduced.dt.1.imputedRefMLFreq.filtered03)*10)

# Get the observed data parameters
n.snps = dim(reduced.dt.1.imputedRefMLFreq.filtered03)[1]
n.pools = length(pool.names)

# Upload BAYPASS *_yij_pij (Reference allele frequency) posterior estimates  
dt.1.posteriorRefEstimates <- read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30_summary_yij_pij.out",
                                          header=TRUE)
# Check BAYPASS Posterior estimates matrix
#head(dt.1.posteriorRefEstimates)

# Convert the table with BAYPASS posteriors to a matrix of nsnps x npops
dt.1.posteriorRefFreq.matrix <- matrix(dt.1.posteriorRefEstimates$M_P, nrow=n.snps, ncol=n.pools, byrow=TRUE)
dt.1.posteriorRefCount.matrix <- matrix(dt.1.posteriorRefEstimates$M_Y, nrow=n.snps, ncol=n.pools, byrow=TRUE)

# check read count approximation with round()
round(dt.1.posteriorRefFreq.matrix, 3)
round(dt.1.posteriorRefCount.matrix)

# REF ALLELE FREQUENCY
# Create a list with Reference AF matrices - ML and Posteriors
list.frequencies <- list(x=reduced.dt.1.imputedRefMLFreq.filtered03, # ML Estimates on read count
                         y=dt.1.posteriorRefFreq.matrix) # posterior estimates

# Calculate the correlation statistics
af.corr <- calculate.multiple.correlations(list.frequencies)
row.names(af.corr) <- pool.names

#           cor   mes rsquared   bias
# MN-B1.1 0.839 0.016    0.701 -0.003
# MN-B4.1 0.793 0.024    0.626 -0.005
# ND-B1.1 0.795 0.025    0.627 -0.003
# ND-B4.1 0.790 0.032    0.606 -0.001
# NW-B1.2 0.794 0.025    0.625  0.001
# NW-B4.1 0.820 0.019    0.670  0.006
# NW-B4.2 0.793 0.029    0.617 -0.004
# NW-B4.3 0.788 0.025    0.617 -0.012
# PA-B1.1 0.773 0.031    0.587  0.003
# PA-B4.1 0.744 0.029    0.545  0.023
# PA-B4.2 0.764 0.028    0.580  0.007
# WI-B1.1 0.854 0.023    0.707 -0.012
# WI-B1.2 0.833 0.017    0.688 -0.016
# WI-B4.1 0.793 0.023    0.627  0.004
# WI-B4.2 0.808 0.022    0.648  0.004
# WI-B4.3 0.773 0.026    0.594  0.010
# WO-B1.1 0.747 0.032    0.555  0.002
# WO-B4.1 0.764 0.032    0.572 -0.019
# WO-B4.2 0.782 0.025    0.607  0.011

# REF ALLELE COUNTS
# Create a list with Reference Read Count matrices - ML and Posteriors
list.counts <- list(x=round(reduced.dt.1.imputedRefMLFreq.filtered03*10), # ML Estimates on read count
                    y=round(dt.1.posteriorRefCount.matrix)) # Posterior estimates

# Calculate the correlation statistics
counts.corr <- calculate.multiple.correlations(list.counts)
row.names(counts.corr) <- pool.names

#           cor   mes rsquared   bias
#MN-B1.1 0.972 0.315    0.942  0.016
#MN-B4.1 0.948 0.719    0.889 -0.005
#ND-B1.1 0.947 0.787    0.884  0.005
#ND-B4.1 0.939 1.177    0.857  0.018
#NW-B1.2 0.933 0.973    0.856  0.019
#NW-B4.1 0.961 0.475    0.917  0.037
#NW-B4.2 0.942 1.029    0.869  0.008
#NW-B4.3 0.938 0.876    0.867 -0.051
#PA-B1.1 0.932 1.152    0.851  0.033
#PA-B4.1 0.928 1.037    0.842  0.104
#PA-B4.2 0.929 1.025    0.850  0.065
#WI-B1.1 0.959 0.736    0.907 -0.022
#WI-B1.2 0.966 0.399    0.928 -0.031
#WI-B4.1 0.941 0.790    0.875  0.028
#WI-B4.2 0.947 0.730    0.886  0.035
#WI-B4.3 0.927 1.009    0.846  0.066
#WO-B1.1 0.924 1.207    0.837  0.035
#WO-B4.1 0.928 1.207    0.841 -0.076
#WO-B4.2 0.940 0.828    0.871  0.064

# Plot Posterior Reference allale count against ML Reference allele count
list.plots <- list()
for (i in 1:19)
{
  x=round(reduced.dt.1.imputedRefMLFreq.filtered03*10)[,i]
  y=round(dt.1.posteriorRefCount.matrix)[,i]
  
  df <- data.frame(x = x,y = y)
  list.plots[[i]] <- ggplot(data = df,aes(x = x,y = y)) + stat_sum() + 
    ggtitle (pool.names[i]) +
    xlab("ML Ref Count") +
    ylab("Posterior Ref Count") +
    xlim(0,10) +
    ylim(0,10) +
    geom_abline(intercept = 0, slope = 1, color = "black",  linetype="dashed") + 
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       legend.position = "right", axis.text = element_text(size = 12),
                       axis.title=element_text(size=12))
  
  #plot(x=round(reduced.dt.1.imputedRefMLFreq.filtered03*10)[,i],
  #     y=round(dt.1.posteriorRefCount.matrix)[,i], 
  #     xlab="ML Ref AF", ylab="Posterior Ref AF", main=pool.names[i],
  #     xlim=c(0,10), ylim=c(0,10), pch=19,
  #)
  #abline(a=0, b=1, col = "black", lty = 3)
}

p1 <- cowplot::plot_grid(list.plots[[1]],list.plots[[2]],
                         list.plots[[3]],list.plots[[4]],
                         list.plots[[5]],list.plots[[6]],
                         list.plots[[7]],list.plots[[8]],
                         list.plots[[9]],list.plots[[10]],
                         nrow = 5,
                         ncol = 2,
                         labels = NULL)
save_plot(filename = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/refcount_correlation_p1.pdf", 
          p1, ncol = 2, base_asp = 1.1,
          base_width=8, base_height=14)

p2 <- cowplot::plot_grid(list.plots[[11]],list.plots[[12]],
                         list.plots[[13]],list.plots[[14]],
                         list.plots[[15]],list.plots[[16]],
                         list.plots[[17]],list.plots[[18]],
                         list.plots[[19]],
                         nrow = 5,
                         ncol = 2,
                         labels = NULL)

save_plot(filename = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/refcount_correlation_p2.pdf", 
          p2, ncol = 2, base_asp = 1.1,
          base_width=8, base_height=14)


### WITH DISTANCE ESTIMATES
## Pairwise FST(POOLFSTAT) vs Omega BAYPASS'

# upload poolfstat pairwise fst matrix
load(file = "results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/s.dt.1.pfst.RData")
pfst.matrix <- s.dt.1.pfst@PairwiseFSTmatrix
pfst.matrix[is.na(pfst.matrix)] <- 0

# upload the estimated Omega matrix
omega.matrix=as.matrix(read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30_mat_omega.out"))
dimnames(omega.matrix)=list(pool.names,pool.names)

## COMPARING SVD pFST and DISTANCE MATRIX
# Using SVD decomposition values
pfst.svd <- svd(pfst.matrix)
omega.svd <- svd(omega.matrix)

list.dist.svf <- list(x=pfst.svd$u, y=omega.svd$u)
svd.corr <- calculate.multiple.correlations(list.dist.svf)

# Using SVD plots
pfst.svd <- plot.omega.mod(omega=pfst.matrix,pop.names=pool.names, col=biotype.col)
omega.svd <- plot.omega.mod(omega=omega.matrix,pop.names=pool.names, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")

######
### BAYPASS ANALYSES
######

### DEMOGRAPHY

## Upload OMEGA MATRIX (repeated from above)
omega.matrix=as.matrix(read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30_mat_omega.out"))
dimnames(omega.matrix)=list(pool.names,pool.names)

## Visualization of the OMEGA matrix
# Using SVD decomposition
plot.omega.mod(omega=omega.matrix, pop.names=pool.names, col=biotype.col)
abline(v=0, h=0, lty=3,col = "gray60")

# as a correlation plot
cor.matrix=cov2cor(omega.matrix)
dimnames(cor.matrix)=list(pool.names,pool.names)

corrplot(cor.matrix, method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

# as a heatmap and hierarchical clustering tree (using the average agglomeration method)
hclust.ave <- function(x) hclust(x, method="average")
heatmap(1-cor.matrix, hclustfun = hclust.ave,
        main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

#plot.omega.mod(omega=(1-cor.matrix), pop.names=pool.names, col=biotype.col)
#abline(v=0, h=0, lty=3,col = "gray60")


### SELECTION

## BAYPASS XtX ESTIMATES
# UploadXtX differentiation estimates (using the calibrated XtXst estimator)
core.snp.xtx=read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30_summary_pi_xtx.out", 
                        header=TRUE)

# check estimated XtXst distribution
hist(core.snp.xtx$XtXst, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "XtX",
     main = "")
lines(density(core.snp.xtx$XtXst), # density plot
      lwd = 2, # thickness of line
      col = "red")

# check behavior of the p-values associated to the estimated XtXst
hist(10**(-1*core.snp.xtx$log10.1.pval.),freq=F,breaks=50)
abline(h=1) 
## CONCLUSION :: P-values don't behave well

# Manhattan plots of XtXst and p-values
layout(matrix(1:2,2,1))
plot(core.snp.xtx$XtXst)
plot(core.snp.xtx$log10.1.pval.,ylab="XtX P-value (-log10 scale)")
abline(h=3,lty=2, col="red") #0.001 p--value theshold

## CALIBRATE THE XtX STATISTICS WITH NORMALIZE TRANSFORMATION

# Posterior mean of XtX
xtx=core.snp.xtx$M_XtX

hist(xtx, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "XtX",
     main = "")
lines(density(xtx), # density plot
      lwd = 2, # thickness of line
      col = "red")

# Posterior mean of the XtX NORMALIZED
xtx.std= (core.snp.xtx$M_XtX - mean(core.snp.xtx$M_XtX))/(sd(core.snp.xtx$M_XtX))

hist(xtx.std, # histogram
     col="green", # column color
     border="black",
     breaks = 100,
     prob = TRUE, # show densities instead of frequencies
     xlab = "Standardized XtX",
     main = "")
lines(density(xtx.std), # density plot
      lwd = 2, # thickness of line
      col = "red")

# Calculate empirical p-value
#empval.xtx = empPvals(stat = xtx.std, stat0 = xtx.std[which(xtx.std <= 4)]) #With qvalue package
#qqplot(x=empval.xtx, y=10**(-1*core.snp.xtx$log10.1.pval.))
#abline(a = 0, b = 1)

#thr.pval1ppm = min(xtx[empval.xtx < 0.001]) 
#xtx.qval = qvalue(empval.xtx)
#thr.FDR10pcent = min(xtx[xtx.qval$qvalue<0.1])

## Plot XtX with empirically calculate threshold
#plot(core.snp.xtx$M_XtX)
#abline(h=thr.FDR10pcent,lty=2, col="red")

# compute the 1% threshold
xtx.std.thresh = quantile(xtx.std, probs=c(0.999))

# Plot XtX with empirically calculate threshold
plot(xtx.std)
abline(h=xtx.std.thresh,lty=2, col="red")

# Which SNPS were signigicant?
#snps.inf <- read.table(file = "results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aphidpool.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30.snpdet", h=F)
#
#snps.inf[which(core.snp.xtx$M_XtX > thr.FDR10pcent), ] # positive selection?
#snps.inf[which(core.snp.xtx$M_XtX < 18), ] # balanced selection?
#
#write.table(snps.inf[which(core.snp.xtx$M_XtX > thr.FDR10pcent), c(1:2)],
#            file = "baypass_outliers.txt",
#            col.names = F, row.names = F, quote = F)

## CALIBRATE THE XtX STATISTICS WITH PODs
# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution
pi.beta.coef = read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30_summary_beta_params.out",
                        header = T)$Mean

# upload the original data to obtain total allele count
pooled.data<-geno2YN("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/aphidpool.PoolSeq.PoolSNP.01.5.28Jun2021.reduced30.genobaypass")

# Create the POD
#setwd("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/")
simu.pooled <- simulate.baypass(omega.mat=omega.matrix, nsnp=100000, sample.size = 10, coverage = pooled.data$NN, 
                                beta.pi=pi.beta.coef, pi.maf=0, suffix="core_pods")
#setwd("/fs/scratch/PAS1715/aphidpool")

# RUN BAYPASS WITH POD DATA
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
# wd=/fs/scratch/PAS1715/aphidpool
# sbatch ${wd}/run_baypass_pods.sh

#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################

# get estimate of omega from the POD analysis
pod.omega = as.matrix(read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.pods_mat_omega.out"))
plot(pod.omega, omega.matrix) ; abline(a=0,b=1)
fmd.dist(pod.omega,omega.matrix)

# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.pods_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef, pi.beta.coef)  ; abline(a=0,b=1)

#######################################################
#XtX calibration
#######################################################
#get the pod XtX
pod.xtx = read.table("results/aggregated_data/minmaxcov_4_99/baypass_poolsnp/core.pods_summary_pi_xtx.out",h=T)$M_XtX

# compute the 1% threshold
pod.thresh = quantile(pod.xtx, probs=c(0.999))

# Compute the p-values of XtX calibrated with PODS XtX (From Gautier 2017 presentation)
#empval.xtx.pods = empPvals(stat = xtx, stat0 = pod.xtx) #With qvalue package
#qqplot(x=empval.xtx.pods, y=10**(-1*core.snp.xtx$log10.1.pval.))
#abline(a = 0, b = 1)

#add the thresh to the actual XtX plot
plot(core.snp.xtx$M_XtX)
abline(h=pod.thresh,lty=2, col="red")



