##########################################################
#                   SelEstim package                     # 
#              DEMOGRAPHY AND GENOME SCAN                #
##########################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

setwd("/fs/project/PAS1554/aphidpool")
renv::restore()
#renv::snapshot()
## Remove last features
rm(list=ls())
ls()

library(gplots)
library(coda)
library(lattice)
#library(corrplot)
library(ggplot2)
#library(cowplot)
library(qvalue)
library(fitdistrplus)
library(invgamma)
library(tidyverse)
library(gridExtra)
library(scales)
library(ggpubr) # cannot install
library(ggsignif)
library(xtable)
#require(plyr)
library(psych)
library(viridis)
library(VennDiagram)

### Recover R-renv environment
setwd("/fs/project/PAS1554/aphidpool")
### Import auxiliary R functions
source("/users/PAS1554/vitorpavinato/local/SelEstim/1.1.7/R/SelEstim.R")
source("AglyPoolseq/analyses/aggregated_data/aux_func.R")

### Recover R-renv environment
#setwd("/fs/scratch/PAS1715/aphidpool")

### Working dir
workingdir <- "results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/"

poolnames <- c("MN-Av.1", "MN-V.1",
               "ND-Av.1", "ND-V.1", 
               "NW-Av.1", "NW-Av.2", "NW-V.1", "NW-V.2", "NW-V.3", 
               "PA-Av.1", "PA-V.1", "PA-V.2", 
               "WI-Av.1", "WI-Av.2", "WI-V.1", "WI-V.2", "WI-V.3", 
               "WO-Av.1", "WO-V.1", "WO-V.2", "WO-V.3")

poolbiotypes <- c("Av", "V",
                  "Av", "V", 
                  "Av", "Av", "V", "V", "V", 
                  "Av", "V", "V", 
                  "Av", "Av", "V", "V", "V", 
                  "Av", "V", "V", "V")

####
####
#### ---- RANDOMIZE REFERENCE ALLELE ----
####
####
#
#### SelEstim.R::randomize.reference.allele
#randomize.reference.allele(infile = paste0(workingdir,"aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.genoselestim"),
#                           outfile = paste0(workingdir,"aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021_rand.genoselestim"),
#                           pool = TRUE)
#
####
####
#### ---- CHECK MCMC MIXIMING WITH A SMALL PILOT RUN  ----
####
####
#
#burnin.length = 50000
#mcmc.length = 100000
#thin = 40
#
#burnin.length/thin #1250
#mcmc.length/thin #2500
#
#### FROM SelEstim Run, load the tracer files
#### 1 - Visually inspect mixing
#### 2 - Estimate autocorrelation
#
#### TRACER FILE 1: LAMBDA PARAMETER
#lambda_tracer <- read.table(file=paste0(workingdir, "output_pilotrun/trace_lambda.out"),
#                            header = T)
#
#lambda_mcmc <- mcmc(data=lambda_tracer$lambda, thin = thin)
#
###Plot trace and density for lambda
#densityplot(lambda_mcmc)
#xyplot(lambda_mcmc)
#plot(lambda_mcmc) #combined plots
#acfplot(lambda_mcmc, thin = thin)
#autocorr(lambda_mcmc)
#
#
#### TRACER FILE 2: BETA HYPER PARAMETERS
#beta_tracer <- read.table(file=paste0(workingdir, "output_pilotrun/trace_beta.out"),
#                            header = T)
#
#beta_a <- mcmc(beta_tracer[,2], thin = thin)
#beta_b <- mcmc(beta_tracer[,3], thin = thin)
#beta_mcmc <- mcmc.list(beta_a, beta_b)
#
###Plot trace and density for lambda
#densityplot(beta_mcmc)
#xyplot(beta_mcmc)
#plot(beta_mcmc) #combined plots
#acfplot(beta_mcmc, thin = thin)
#autocorr(beta_mcmc)
#
#### TRACER FILE 3: M parameter
#m_tracer <- read.table(file=paste0(workingdir, "output_pilotrun/trace_M.out"),
#                       header = T)
#
#deme_1_mcmc <- mcmc(m_tracer[,2], thin = thin)
#deme_2_mcmc <- mcmc(m_tracer[,3], thin = thin)
#deme_3_mcmc <- mcmc(m_tracer[,4], thin = thin)
#deme_4_mcmc <- mcmc(m_tracer[,5], thin = thin)
#deme_5_mcmc <- mcmc(m_tracer[,6], thin = thin)
#deme_6_mcmc <- mcmc(m_tracer[,7], thin = thin)
#deme_7_mcmc <- mcmc(m_tracer[,8], thin = thin)
#deme_8_mcmc <- mcmc(m_tracer[,9], thin = thin)
#deme_9_mcmc <- mcmc(m_tracer[,10], thin = thin)
#deme_10_mcmc <- mcmc(m_tracer[,11], thin = thin)
#deme_11_mcmc <- mcmc(m_tracer[,12], thin = thin)
#deme_12_mcmc <- mcmc(m_tracer[,13], thin = thin)
#deme_13_mcmc <- mcmc(m_tracer[,14], thin = thin)
#deme_14_mcmc <- mcmc(m_tracer[,15], thin = thin)
#deme_15_mcmc <- mcmc(m_tracer[,16], thin = thin)
#deme_16_mcmc <- mcmc(m_tracer[,17], thin = thin)
#deme_17_mcmc <- mcmc(m_tracer[,18], thin = thin)
#deme_18_mcmc <- mcmc(m_tracer[,19], thin = thin)
#deme_19_mcmc <- mcmc(m_tracer[,20], thin = thin)
#deme_20_mcmc <- mcmc(m_tracer[,21], thin = thin)
#deme_21_mcmc <- mcmc(m_tracer[,22], thin = thin)
#
#m_mcmc <- mcmc.list(deme_1_mcmc, deme_2_mcmc, deme_3_mcmc, deme_4_mcmc, deme_5_mcmc, deme_6_mcmc,
#                    deme_7_mcmc, deme_8_mcmc, deme_9_mcmc, deme_10_mcmc, deme_11_mcmc, deme_12_mcmc,
#                    deme_13_mcmc, deme_14_mcmc, deme_15_mcmc, deme_16_mcmc, deme_17_mcmc, deme_18_mcmc,
#                    deme_19_mcmc, deme_20_mcmc, deme_21_mcmc)
#
###Plot trace and density for lambda
#densityplot(m_mcmc)
#xyplot(m_mcmc)
#plot(m_mcmc) #combined plots
#acfplot(m_mcmc, thin = thin)
#autocorr(m_mcmc)


###
###
### ---- CHECK MCMC CHAINS CONVERGENCY BETWEEN THREE INDEPENDENT RUNS ----
### 09/FEB/2022
###

burnin.length = 120000
mcmc.length = 50000
thin = 50

burnin.length/thin #2400
mcmc.length/thin #1000

### FROM SelEstim Run, load the tracer files to:
### 1 - Visually inspect mixing
### 2 - Estimate autocorrelation

### TRACER FILE 1: LAMBDA PARAMETER
lambda_files <-  paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "summary_lambda.out")
lambda <- lapply(lambda_files, function(x){read.table(file= x, header=T, na.strings = "NA")})

trace_lamda_files <- paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "trace_lambda.out")
trace_lamda <- lapply(trace_lamda_files, function(x){read.table(file= x, header=T, na.strings = "NA")})
  
run_lambdas <- cbind(trace_lamda[[1]]$lambda,
                     trace_lamda[[2]]$lambda,
                     trace_lamda[[3]]$lambda)

colnames(run_lambdas) <- c("lambda_0", "lambda_1", "lambda_2")

## Expected value for every run
trace_lambda_medians <- apply(run_lambdas, 2, function(x){median(x)})

## Expected value for the median value distribution
median_trace_lambda <- apply(run_lambdas, 1, function(x){median(x)})
mean(median_trace_lambda)

## Plot \lambda of every run
plot(trace_lamda[[1]]$iteration, median_trace_lambda, type="l")
lines(trace_lamda[[1]]$iteration, trace_lamda[[1]]$lambda, type="l", col="red")
lines(trace_lamda[[2]]$iteration, trace_lamda[[2]]$lambda, type="l", col="green")
lines(trace_lamda[[3]]$iteration, trace_lamda[[3]]$lambda, type="l", col="blue")

## Plot \lambda prior and hyperprior
plot(density(rinvgamma(1000, shape=3, rate = 2, scale = 1/rate)))
lines(density(run_lambdas[,1]))
lines(density(run_lambdas[,2]))
lines(density(run_lambdas[,3]))

## MCMC convergency sanity check
mcmc_lambda <- lapply(trace_lamda, function(x){mcmc(data=x, thin=thin)})   
mcmc_lambda_list <- mcmc.list(mcmc_lambda[[1]], mcmc_lambda[[2]], mcmc_lambda[[3]])

## Plot trace and density for lambda
densityplot(mcmc_lambda_list, thin = thin)
xyplot(mcmc_lambda_list)
plot(mcmc_lambda_list) #combined plots
acfplot(mcmc_lambda_list, thin = thin)
autocorr(mcmc_lambda_list)

### Chains did not mixture. 
### They also have high autocorrelation, but they don't seem to show a decrease trend
### CONCLUSION: Analyse separately, and consider selected loci when they appear in all independent analysis. 
### Average 2Ns and global selection parameters \delta of the 3 runs
### I cannot run with fixed value of \lambda

### TRACER FILE 2: BETA HYPER PARAMETERS
beta_files <-  paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "summary_beta.out")
beta <- lapply(beta_files, function(x){read.table(file= x, header=T, na.strings = "NA")})

trace_beta_files <- paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "trace_beta.out")
trace_beta <- lapply(trace_beta_files, function(x){read.table(file= x, header=T, na.strings = "NA")})

## BETA A
run_beta_a <- cbind(trace_beta[[1]]$beta_a,
                    trace_beta[[2]]$beta_a,
                    trace_beta[[3]]$beta_a)

colnames(run_beta_a) <- c("betaA_0", "betaA_1", "betaA_2")

# Expected value for every run
trace_betaA_medians <- apply(run_beta_a, 2, function(x){median(x)})

## Expected value for the median value distribution
median_trace_betaA <- apply(run_beta_a, 1, function(x){median(x)})
mean(median_trace_betaA)

## Plot \beta A of every run
plot(trace_beta[[1]]$iteration, median_trace_betaA, type="l")
lines(trace_beta[[1]]$iteration, trace_beta[[1]]$beta_a, type="l", col="red")
lines(trace_beta[[2]]$iteration, trace_beta[[2]]$beta_a, type="l", col="green")
lines(trace_beta[[3]]$iteration, trace_beta[[3]]$beta_a, type="l", col="blue")

## MCMC convergency sanity check
mcmc_betaA <- lapply(trace_beta, function(x){mcmc(data=x$beta_a, thin=thin)})   
mcmc_betaA_list <- mcmc.list(mcmc_betaA[[1]], mcmc_betaA[[2]], mcmc_betaA[[3]])

## Plot trace and density for \beta A
densityplot(mcmc_betaA_list, thin = thin)
xyplot(mcmc_betaA_list)
plot(mcmc_betaA_list) #combined plots
acfplot(mcmc_betaA_list, thin = thin)
autocorr(mcmc_betaA_list)

## BETA B
run_beta_b <- cbind(trace_beta[[1]]$beta_b,
                    trace_beta[[2]]$beta_b,
                    trace_beta[[3]]$beta_b)

colnames(run_beta_b) <- c("betaB_0", "betaB_1", "betaB_2")

# Expected value for every run
trace_betaB_medians <- apply(run_beta_b, 2, function(x){median(x)})

## Expected value for the median value distribution
median_trace_betaB <- apply(run_beta_b, 1, function(x){median(x)})
mean(median_trace_betaB)

## Plot \beta B of every run
plot(trace_beta[[1]]$iteration, median_trace_betaB, type="l")
lines(trace_beta[[1]]$iteration, trace_beta[[1]]$beta_b, type="l", col="red")
lines(trace_beta[[2]]$iteration, trace_beta[[2]]$beta_b, type="l", col="green")
lines(trace_beta[[3]]$iteration, trace_beta[[3]]$beta_b, type="l", col="blue")

## MCMC convergency sanity check
mcmc_betaB <- lapply(trace_beta, function(x){mcmc(data=x$beta_b, thin=thin)})   
mcmc_betaB_list <- mcmc.list(mcmc_betaB[[1]], mcmc_betaB[[2]], mcmc_betaB[[3]])

## Plot trace and density for lambda
densityplot(mcmc_betaB_list, thin = thin)
xyplot(mcmc_betaB_list)
plot(mcmc_betaB_list) #combined plots
acfplot(mcmc_betaB_list, thin = thin)
autocorr(mcmc_betaB_list)

## Combined List
mcmc_beta_list <- mcmc.list(mcmc_betaA[[1]], mcmc_betaA[[2]], mcmc_betaA[[3]],
                            mcmc_betaB[[1]], mcmc_betaB[[2]], mcmc_betaB[[3]])

##Plot trace and density for lambda
densityplot(mcmc_beta_list, thin = thin)
xyplot(mcmc_beta_list)
plot(mcmc_beta_list) #combined plots
acfplot(mcmc_beta_list, thin = thin)
autocorr(mcmc_beta_list)

## Plot \beta prior
plot(density(rbeta(1000, shape1 = beta[[1]][1,1], shape2 = beta[[1]][2,1])))
lines(density(rbeta(1000, shape1 = beta[[2]][1,1], shape2 = beta[[2]][2,1])))
lines(density(rbeta(1000, shape1 = beta[[3]][1,1], shape2 = beta[[3]][2,1])))

### Chains did mixture.
### Posterior estimates of \beta A and B were similar across runs.
### CONCLUSION: If run, use fixed values of shape 1 and 2 (beta a and b).

### TRACER FILE 3: M parameter
trace_M_files <- paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "trace_M.out")
trace_M <- lapply(trace_M_files, function(x){read.table(file= x, header=T, na.strings = "NA")})

## Look for convergence of M among 3 runs
mcmc_M_deme1 <- lapply(trace_M, function(x){mcmc(data=x$deme_1, thin=thin)})   
mcmc_M_deme1_list <- mcmc.list(mcmc_M_deme1[[1]], mcmc_M_deme1[[2]], mcmc_M_deme1[[3]])

mcmc_M_deme5 <- lapply(trace_M, function(x){mcmc(data=x$deme_5, thin=thin)})   
mcmc_M_deme5_list <- mcmc.list(mcmc_M_deme5[[1]], mcmc_M_deme5[[2]], mcmc_M_deme5[[3]])

mcmc_M_deme10 <- lapply(trace_M, function(x){mcmc(data=x$deme_10, thin=thin)})   
mcmc_M_deme10_list <- mcmc.list(mcmc_M_deme10[[1]], mcmc_M_deme10[[2]], mcmc_M_deme10[[3]])

##Plot trace and density for lambda
densityplot(mcmc_M_deme10_list, thin = thin)
xyplot(mcmc_M_deme10_list)
plot(mcmc_M_deme10_list) #combined plots
acfplot(mcmc_M_deme10_list, thin = thin)
autocorr(mcmc_M_deme10_list)

## Savepoint_0
##-------------
#save.image("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
#load("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
############## END THIS PART

###
###
### ---- KDL ANALYSIS OF THREE INDEPENDENT RUNS ----
### 09/FEB/2022
###

### Plot locus-selection coefficients delta
pdf("results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis//delta_runs.pdf", # File name
    width = 11, height = 8.50,                                                                # Width and height in inches
    bg = "white",                                                                             # Background color
    colormodel = "cmyk",                                                                      # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(3, 1))
plot.delta(file = paste0(workingdir,"output_0/","summary_delta.out"))
plot.delta(file = paste0(workingdir,"output_1/","summary_delta.out"))
plot.delta(file = paste0(workingdir,"output_2/","summary_delta.out"))
dev.off()

### Plot Kld
pdf("results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis//kld_runs.pdf", # File name
    width = 11, height = 8.50,                                                                # Width and height in inches
    bg = "white",                                                                             # Background color
    colormodel = "cmyk",                                                                      # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1, mfrow=c(3, 1))
plot.kld(file = paste0(workingdir,"output_0/","summary_delta.out"), 
         calibration_file = paste0(workingdir,"output_0/","calibration/","summary_delta.out"), limit = 0.001)
plot.kld(file = paste0(workingdir,"output_1/","summary_delta.out"),
         calibration_file = paste0(workingdir,"output_1/","calibration/","summary_delta.out"), limit = 0.001)
plot.kld(file = paste0(workingdir,"output_2/","summary_delta.out"),
         calibration_file = paste0(workingdir,"output_2/","calibration/","summary_delta.out"), limit = 0.001)
dev.off()


###
### -- Identify significant selected SNPs --
###

## Load summary_delta.out for three independent SELEstim runs
summary_delta_files <- paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "summary_delta.out")
summary_delta <- lapply(summary_delta_files, function(x){read.table(file= x, header=T, na.strings = "NA")})

## Load the information of SNP loci
loci_info  <- read.table(file = paste0(workingdir,"aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.snpdet")) 
names(loci_info) <- c("CHROM","POS", "REF", "ALT")

## Combine the delta_summary.out with the SNP info
loci_deltas <- lapply(summary_delta, function(x){cbind(loci_info, x)})

colnames(loci_deltas[[1]]) <- c("CHROM","POS", "REF", "ALT","locus","mean","std","KLD")
colnames(loci_deltas[[2]]) <- c("CHROM","POS", "REF", "ALT","locus","mean","std","KLD")
colnames(loci_deltas[[3]]) <- c("CHROM","POS", "REF", "ALT","locus","mean","std","KLD")

### Calculate FDR thresholds from the data for each independent run
## Posterior mean of KLD
#kld_0 = loci_deltas[[1]]$KLD
#kld_1 = loci_deltas[[2]]$KLD
#kld_2 = loci_deltas[[3]]$KLD
#
## Posterior mean of the KLD NORMALIZED
#kld_0.std = (kld_0 - mean(kld_0))/(sd(kld_0))
#kld_1.std = (kld_1 - mean(kld_1))/(sd(kld_1))
#kld_2.std = (kld_2 - mean(kld_2))/(sd(kld_2))
#
## Calculate empirical p-value
#kld_0_empval = empPvals(stat = kld_0.std, stat0 = kld_0.std[which(kld_0.std <= quantile(kld_0, 0.9999))]) #With qvalue package
## Threshold based on lowest empirical pvalue
#thr_pval1ppm_0 = min(kld_0[kld_0_empval < 0.01]) #KLD
#thr_pval1ppm_0.std = min(kld_0.std[kld_0_empval < 0.01]) #standardized KLD
#
#kld_1_empval = empPvals(stat = kld_1.std, stat0 = kld_1.std[which(kld_1.std <= quantile(kld_1, 0.9999))]) #With qvalue package
## Threshold based on lowest empirical pvalue
#thr_pval1ppm_1 = min(kld_1[kld_1_empval < 0.01]) #KLD
#thr_pval1ppm_1.std = min(kld_1.std[kld_1_empval < 0.01]) #standardized KLD
#
#kld_2_empval = empPvals(stat = kld_2.std, stat0 = kld_2.std[which(kld_2.std <= quantile(kld_2, 0.9999))]) #With qvalue package
## Threshold based on lowest empirical pvalue
#thr_pval1ppm_2 = min(kld_2[kld_2_empval < 0.01]) #KLD
#thr_pval1ppm_2.std = min(kld_2.std[kld_2_empval < 0.01]) #standardized KLD
#
#kld_thr1ppm <- c(thr_pval1ppm_0, thr_pval1ppm_1, thr_pval1ppm_2)
#kld_thr1ppm.std <- c(thr_pval1ppm_0.std, thr_pval1ppm_1.std, thr_pval1ppm_2.std)

## Select SNPs with KDL > 95% quantile KDL in each run
top_kld_snps <- lapply(loci_deltas, function(x){x[which(x[,8] > quantile(x[,8], 0.95)), ]})

## Save Top KDL SNPs names to a file
for (i in 1:3)
{
  tmp <- paste0(top_kld_snps[[i]]$CHROM,"_", top_kdl_snps[[i]]$POS)
  
  write.table(tmp,
              file = paste0(workingdir, "analysis/" ,"top5perc_empKld_snps", "_", (i-1), ".txt"),
              col.names = F, row.names = F, quote = F)
}

## HOW MANY SNPS ARE IN COMMON AMONG THE INDEPENDENT RUNS?
l0 <- paste0(top_kld_snps[[1]]$CHROM,"_", top_kld_snps[[1]]$POS) #13144
l1 <- paste0(top_kld_snps[[2]]$CHROM,"_", top_kld_snps[[2]]$POS) #13143
l2 <- paste0(top_kld_snps[[3]]$CHROM,"_", top_kld_snps[[3]]$POS) #13144

venn.obj <- venn(list(r0=l0,r1=l1,r2=l2), intersections=TRUE)
intersections <- attr(venn.obj, "intersections")

common_top_kld_snps <- intersections$`r0:r1:r2`

## Save the name of significan snps identified with KLD::SelEstim
write.table(x=common_top_kld_snps, file = "results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/lociIDs_signf_snps_kld_common.txt",
            quote=F, row.names = F, col.names = F)

## Export venn diagram to a file
venn(list(r0=l0,r1=l1,r2=l2), intersections=TRUE)

## PLOT VENN DIAGRAM
#venn.plot <- venn.diagram(list(r0=l0,r1=l1,r2=l2),
#                          category.names = c("run 1", "run 2", "run 3"),
#                          fill = c("#0073C2FF", "#EFC000FF", "#868686FF"), alpha = 0.3, filename = NULL);
#grid.draw(venn.plot)

## Savepoint_1
##-------------
#save.image("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
#load("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
############## END THIS PART

###
###
### ---- DFE ANALYSIS OF THREE INDEPENDENT RUNS ----
### 09/FEB/2022
###

##
## -- delta_j parameter: the locus-specific (global) effect of selection --
##
## Get the median delta values among runs
run_deltas <- cbind(summary_delta[[1]]$mean,
                    summary_delta[[2]]$mean,
                    summary_delta[[3]]$mean)

colnames(run_deltas) <- c("delta_0", "delta_1", "delta_2")

median_delta <- apply(run_deltas, 1, function(x){median(x)})

## Plot density for each run and median values
nsnps = 262866

df1 <- data.frame(run=c(rep("run0", nsnps), rep("run1", nsnps), 
                        rep("run2", nsnps), rep("median", nsnps)),
                 delta=c(run_deltas[,1], run_deltas[,2], 
                         run_deltas[,3], median_delta))

p1 <- ggplot(df1, aes(x=delta, fill=run))+
      geom_density(alpha=.6,) +
      labs(x = expression(italic(delta)[j]), y = "Density") + 
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 14),
                         axis.title=element_text(size=20))
p1

## DFE Goodness of fit of delta parameter

# First trim the data to remove values lower than the 1% and higher than the 99% quantile
trimmed_median_delta <- median_delta[median_delta < quantile(median_delta, 0.99) #& 
                                     #median_delta > quantile(median_delta, 0.01)
                                     ]

# Try to fit a distribution to estimate the parameters from the data
fitted_invgamma <- fitdist(trimmed_median_delta, distr = "invgamma", start=list(shape=2, rate=3))
plot(fitted_invgamma)
print(fitted_invgamma)
# Parameters:
#      estimate   Std. Error
#shape  54.42739  0.1504255
#rate  290.27009  0.8059425

## Plot histogram of median values with theoretical density
df2 <- data.frame(loci=seq_along(trimmed_median_delta), delta=trimmed_median_delta)

p2 <- ggplot(df2, aes(x=delta)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="gray", binwidth = 0.2) +
      stat_function(fun=dinvgamma, args=fitdist(df2$delta, 
                                             distr = "invgamma", start=list(shape=2, rate=3))$estimate) +
      geom_density(alpha=.2, fill="#FF6666") +
      labs(x = expression(italic(delta)[j]), y = "Density") + 
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 14),
                         axis.title=element_text(size=20))
p2

## Savepoint_2
##-------------
#save.image("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
#load("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
############## END THIS PART

##
## -- sigma_ij parameter: the locus-specific and deme-specific scale mutation effect 2Ns --
##

## Import summary_sigma.out of each run
summary_sigma_files <- paste0(workingdir, "output_", seq(from=0,to=2, by=1), "/", "summary_sigma.out")
summary_sigma <- lapply(summary_sigma_files, function(x){read.table(file= x, header=T, na.strings = "NA")})

loci_sigma <- cbind(summary_sigma[[1]][,1:3],
                    summary_sigma[[2]]$mean,
                    summary_sigma[[3]]$mean)

colnames(loci_sigma) <- c("locus","demes","sigma_0","sigma_1","sigma_2")
loci_sigma["median"] <- apply(loci_sigma[,3:5], 1, median)

##
## -- sigma_ij POOL : SCALE SELECTION COEFFICIENT FOR EACH POOL --
##

df3 <- data.frame(loci_sigma[,c(1:2,6)])

df3["names"] <- rep(poolnames, nsnps)
df3["biotype"] <- rep(poolbiotypes, nsnps)

## BOXPLOT OF \sigmas: AS A FUNCTION OF POOL
p3 <- ggplot(df3, aes(x=names, group=names, y=median, fill=biotype)) + 
      geom_boxplot() +
      scale_fill_manual(values = c("#247F00","#AB1A53"),
                        name   = "Biotypes",
                        breaks = c("Av", "V"),
                        labels = c(expression("Avirulent"), "Virulent")) +
      xlab("") +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
      
p3
#_. CONCLUSION: Something is going on with ND_V pool?

## CALCULATE THE PROPORTION OF MUTATIONS IN DIFFERENT SELECTION STRENGTH CLASSES
# Split each pool DFE in discrete intervals
# and produce a barplot with the propotions:
# 1 <= 2Ns < 10 (weak)
# 10 <= 2Ns < 100 (moderate);
# 100 <= 2Ns < 1000 (strong);

# Calculate bins of sigma
df3 <- df3 %>% 
  mutate(demes, 
         bin = cut(median, c(1,10,100,1000), include.lowest = T, right = F)
         )

# Converte to strength
df3["strength"] <- ifelse(df3$bin == "[1,10)", "weak", 
                          ifelse(df3$bin == "[10,100)", "moderate", "strong"))

## BOXPLOT OF \sigmas 2: AS A FUNCTION OF POOL AND SELECTION STRENGTH CLASS
p4 <- ggplot(df3 , aes(x=names, group=names, y=median,fill=biotype)) + 
      geom_boxplot() + 
      facet_grid(~ strength) +
      scale_fill_manual(values = c("#247F00","#AB1A53"),
                        name   = "Biotypes",
                        breaks = c("Av", "V"),
                        labels = c("Avirulent", "Virulent")) +
      xlab("") +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      #ylim(c(0, 10.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(size= 12, angle = 90, vjust = 0.5, hjust=0.5))
     
p4

## INCLUDE THE MUTATION ANNOTATION
# Upload the summary table with code of syn and non-syn mutations (generated in genetic_diversity_poolsnp.R)
summary_annotation_file <- "results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/combined_summary_statistics_annotation_table.txt"
summary_annotation <- read.table(file=summary_annotation_file, header = T, sep = "\t", na.strings = NA)

## Take only synonymous and non-synonymous mutations
# Sanity check the order of SNPs
sum(paste0(summary_annotation$CHROM, "_", summary_annotation$POS) == paste0(summary_annotation$CHROM, "_", summary_annotation$POS)) == nsnps
sum(paste0(summary_annotation$CHROM, "_", summary_annotation$POS) == paste0(loci_info$CHROM, "_", as.integer(loci_info$POS))) == nsnps

# Get the locus id for each mutation type
synonymous_snps <- which(summary_annotation$FITNES_EFFECT == "SS")
nonsynonymous_snps <- which(summary_annotation$FITNES_EFFECT == "NS")
noncoding_snps <- which(is.na(summary_annotation$FITNES_EFFECT))

# Create a vector with mutation type: SS-synonymous, NS-nonsynonymous, NC-nonconding;
df3[df3$locus %in% synonymous_snps, "mutation"] <- "SS"
df3[df3$locus %in% nonsynonymous_snps, "mutation"] <-"NS"
df3[df3$locus %in% noncoding_snps, "mutation"] <- "NC"
#table(df3$demes, df3$mutation)

# BOXPLOT OF \sigmas 3: AS A FUNCTION OF POOL, STRENGTH CLASS AND MUTATION TYPE
p5 <- ggplot(df3 , aes(x=names, group=names, y=median,fill=biotype)) + 
      geom_boxplot() + 
      facet_grid(mutation ~ strength) +
      scale_fill_manual(values = c("#247F00","#AB1A53"),
                         name   = "Biotypes",
                         breaks = c("Av", "V"),
                         labels = c("Avirulent", "Virulent")) +
      xlab("") +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      #ylim(c(0, 10.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(size= 12, angle = 90, vjust = 0.5, hjust=0.5))

p5

## BARPLOT OF THE NUMBER OF MUTATIONS OF DIFFERENT STRENGTH CLASS AND TYPE
# Create a data.frame with counts per bin/strength per pool
df3_sel_type <- as.data.frame(table(df3$names, df3$bin, df3$mutation))
colnames(df3_sel_type) <- c("pools", "bin", "type","counts")

df3_sel_type["logCounts"] <- log10(df3_sel_type$counts)
df3_sel_type$logCounts[is.infinite(df3_sel_type$logCounts)] <- 0

# BARPLOT WITH PROPORTIONS
p6 <- ggplot(df3_sel_type, aes(fill=bin, y=logCounts, x=pools)) + 
      geom_bar(position="stack", stat="identity") +
      facet_grid(type ~ .) +
      scale_fill_viridis(discrete = T) +
      labs(fill=expression(paste(sigma[ij], " bins"))) +
      xlab("") +
      ylab(expression(log[10]("counts"))) +
      #ylim(c(0, 10.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 14),
                         axis.title.y=element_text(size=16),
                         axis.title.x=element_text(size=16),
                         axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=0.5))
p6

## COMBINED PLOT 1
CP1 <- ggarrange(p5, p6,
                 ncol=2,
                 widths = c(2, 1),
                 labels=c("A","B"),
                 common.legend = F,
                 legend = "bottom",
                 font.label = list(size = 28, face = "bold"))
CP1


ggsave("results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/sigma_pools_mutations.pdf",
       CP1,
       device="pdf",
       width=15,
       height=9)


##
## -- ISOLATE CANDIDATE SNPS IN COMMON ON THE 3 INDEPENDET RUNS --
##

## Take the row number
selected_snps_idx <- row.names(summary_annotation[paste0(summary_annotation$CHROM,"_", summary_annotation$POS) %in% common_top_kld_snps, ])

## Filter the data.frame to contain only selected snps
df3.filt <- df3[df3$locus %in% selected_snps_idx, ]

# BOXPLOT OF \sigmas 4: AS A FUNCTION OF POOL, STRENGTH CLASS AND MUTATION TYPE OF SELECTED SNPS
p7 <- ggplot(df3.filt , aes(x=names, group=names, y=median, fill=biotype)) + 
      geom_boxplot() + 
      facet_grid(mutation ~ strength) +
      scale_fill_manual(values = c("#247F00","#AB1A53"),
                        name   = "Biotypes",
                        breaks = c("Av", "V"),
                        labels = c("Avirulent", "Virulent")) +
      xlab("") +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      #ylim(c(0, 10.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(size= 12, angle = 90, vjust = 0.5, hjust=0.5))
    
p7

## Export the information of selected loci
common_top_kld_snps_table <- summary_annotation[paste0(summary_annotation$CHROM,"_", summary_annotation$POS) %in% common_top_kld_snps, ]

write.table(x=common_top_kld_snps_table, file="results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/signf_snps_kld_common_annot_table.txt",
            quote = F, row.names = F, col.names = T)

## Savepoint_3
##-------------
#save.image("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
#load("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
############## END THIS PART

##
## -- AVERAGE SCALE SELECTION COEFFICIENT FOR EACH BIOTYPE --
##

## Average loci \sigmas per biotype
df4 <- df3 %>% group_by(locus, biotype, mutation)  %>%
               summarise(
                 # B1 and B4 averages
                 n = n(),
                 mean = harmonic.mean(median, na.rm = T),
                 sd = sd(median, na.rm = T),
                 se =sd(median, na.rm = T)/sqrt(sum(n()))
               )

# Calculate bins and add annotation for strength
df4["bins"] <- cut(df4$mean,  c(1,10,100,1000), include.lowest = T, right = F)
df4["strength"] <- ifelse(df4$bins == "[1,10)", "weak", 
                          ifelse(df4$bins == "[10,100)", "moderate", "strong"))


## Average loci \sigma across all pools (Global dataset)
df5 <- df3 %>% group_by(locus, mutation) %>%
               summarise(
                 # GLOBAL averaged
                 n = n(),
                 mean = harmonic.mean(median, na.rm = T),
                 sd = sd(median, na.rm = T),
                 se =sd(median, na.rm = T)/sqrt(sum(n()))
               )
# Calculate bins and add annotation for strength
df5["bins"] <- cut(df5$mean,  c(1,10,100,1000), include.lowest = T, right = F)
df5["strength"] <- ifelse(df5$bins == "[1,10)", "weak", 
                          ifelse(df5$bins == "[10,100)", "moderate", "strong"))

## Combined biotypes with global dataset
df6 <- rbind(df4,df5)
df6$biotype[is.na(df6$biotype)] <- "Global"

df6["alphaCode"] <- ifelse(df6$biotype == "Av", "A", 
                           ifelse(df6$biotype == "V", "B", "C")) 

## BOXPLOT OF \sigmas: AS A FUNCTION OF BIOTYPE
p8 <- ggplot(df6, aes(x=alphaCode, y=mean, fill=alphaCode)) + 
      geom_boxplot() +
      scale_fill_manual(values = c("#247F00","#AB1A53", "#636363"),
                        name   = "Biotypes",
                        breaks = c("A", "B", "C"),
                        labels = c("Avirulent", "Virulent", "Global")) +
      xlab("") +
      scale_x_discrete(labels = c("", "", "")) +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
     
p8

## BOXPLOT OF \sigmas 2: AS A FUNCTION OF BIOTYPE AND SELECTION STRENGTH CLASS
p9 <- ggplot(df6, aes(x=alphaCode, y=mean,fill=alphaCode)) + 
      geom_boxplot() + 
      facet_grid(~ strength) +
      scale_fill_manual(values = c("#247F00","#AB1A53", "#636363"),
                        name   = "Biotypes",
                        breaks = c("A", "B", "C"),
                        labels = c("Avirulent", "Virulent", "Global")) +
      xlab("") +
      scale_x_discrete(labels = c("", "", "")) +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      #ylim(c(0, 10.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(size= 12, angle = 90, vjust = 0.5, hjust=0.5))
      
p9

# BOXPLOT OF \sigmas 3: AS A FUNCTION OF BIOTYPE, STRENGTH CLASS AND MUTATION TYPE
p10 <- ggplot(df6 , aes(x=alphaCode, y=mean,fill=alphaCode)) + 
      geom_boxplot() + 
      facet_grid(mutation ~ strength) +
      scale_fill_manual(values = c("#247F00","#AB1A53", "#636363"),
                        name   = "Biotypes",
                        breaks = c("A", "B", "C"),
                        labels = c("Avirulent", "Virulent", "Global")) +
      xlab("") +
      scale_x_discrete(labels = c("", "", "")) +
      ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
      #ylim(c(0, 10.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                         axis.line = element_line(colour = "black"),
                         axis.text = element_text(size = 12),
                         axis.title=element_text(size=20),
                         axis.text.x = element_text(size= 12, angle = 90, vjust = 0.5, hjust=0.5))
      
p10

## BARPLOT OF THE NUMBER OF MUTATIONS OF DIFFERENT STRENGTH CLASS AND TYPE
# Create a data.frame with counts per bin/strength per pool
df6_sel_type <- as.data.frame(table(df6$alphaCode, df6$bins, df6$mutation))
colnames(df6_sel_type) <- c("pools", "bin", "type","counts")
df6_sel_type["logCounts"] <- log10(df6_sel_type$counts)
df6_sel_type$logCounts[is.infinite(df6_sel_type$logCounts)] <- 0

# BARPLOT WITH PROPORTIONS
p11 <- ggplot(df6_sel_type, aes(fill=bin, y=logCounts, x=pools)) + 
       geom_bar(position="stack", stat="identity") +
       facet_grid(type ~ .) +
       scale_fill_viridis(discrete = T) +
       labs(fill=expression(paste(sigma[ij], " bins"))) +
       xlab("") +
       scale_x_discrete(labels = c("Avirulent", "Virulent", "Global")) +
       ylab(expression(log[10]("counts"))) +
       #ylim(c(0, 10.5)) +
       theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                          axis.line = element_line(colour = "black"),
                          axis.text = element_text(size = 14),
                          axis.title.y=element_text(size=16),
                          axis.title.x=element_text(size=16),
                          axis.text.x = element_text(size=12, angle = 0, vjust = 0.5, hjust=0.5))
p11

CP2 <- ggarrange(p10, p11,
                 ncol=2,
                 widths = c(2, 1),
                 labels=c("A","B"),
                 common.legend = F,
                 legend = "bottom",
                 font.label = list(size = 28, face = "bold"))
CP2


ggsave("results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/sigma_biotypes_mutations.pdf",
       CP2,
       device="pdf",
       width=15,
       height=9)


## Filter the data.frame to contain only selected snps
df6.filt <- df6[df6$locus %in% selected_snps, ]

# BOXPLOT OF \sigmas 3: AS A FUNCTION OF BIOTYPE, STRENGTH CLASS AND MUTATION TYPE
p12 <- ggplot(df6.filt , aes(x=alphaCode, y=mean,fill=alphaCode)) + 
       geom_boxplot() + 
       facet_grid(mutation ~ strength) +
       scale_fill_manual(values = c("#247F00","#AB1A53", "#636363"),
                         name   = "Biotypes",
                         breaks = c("A", "B", "C"),
                         labels = c("Avirulent", "Virulent", "Global")) +
       xlab("") +
       scale_x_discrete(labels = c("", "", "")) +
       ylab(expression(paste(italic(sigma)[ij], " (2",italic(Ns)[ij],")"))) +
       #ylim(c(0, 10.5)) +
       theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                          axis.line = element_line(colour = "black"),
                          axis.text = element_text(size = 12),
                          axis.title=element_text(size=20),
                          axis.text.x = element_text(size= 12, angle = 90, vjust = 0.5, hjust=0.5))
     
p12

## Savepoint_4
##-------------
#save.image("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
#load("/fs/project/PAS1554/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp/analysis/selestim.poolsnp.workspace09Feb22.RData")
############## END THIS PART

## FINDIND THE BEST FIT DISTRIBUTION OF DELTA_j VALUES
#dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
#pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
#qgumbel <- function(p, a, b) a-b*log(-log(p))