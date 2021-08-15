##########################################################
#       GENOME-WIDE AND SCAFFOLD GENETIC DIVERSITY       # 
#                   FROM POPOOLATION                     #
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
library(tidyverse)
library(gridExtra)
library(scales)
library(ggpubr)
library(ggplot2)

ALPHA=0.75

header4pi <- c("chr", "pos", "nsnps", "frac_cov", "pi")
header4theta <- c("chr", "pos", "nsnps", "frac_cov", "theta")

## LOAD ONE DATASET FOR TESTING
pi <- read.table(file = "results/aggregated_data/minmaxcov_4_99/popoolation_sumstats/maf001_mac1_biallelic/MN_BIO1_S1/MN_BIO1_S1.window.pi",
                 header = F, sep = "\t", col.names = header4pi)

pi[pi$nsnps==0, "pi"] <- NA
mean(pi$pi, na.rm = T)

mean(pi$pi[pi$pi < 0.1], na.rm = T)

plot(pi ~ pos, data = pi[pi$chr == "scaffold_14", ])


pi_2 <- read.table(file = "popoolation_sumstats/MN_BIO1_S1/MN_BIO1_S1.window.pi",
                 header = F, sep = "\t", col.names = header4pi)

pi_2[pi_2$nsnps==0, "pi"] <- NA
mean(pi_2$pi, na.rm = T)

hist(pi_2$pi[pi_2$pi < 0.1])
max(pi_2$pi[pi_2$pi < 0.1], na.rm = T)


mean(pi_2$pi[pi_2$pi < 0.1], na.rm = T)

par(mfrow=c(2, 1))
plot(pi ~ pos, data = pi[pi$chr == "scaffold_14", ], type = "l")
plot(pi ~ pos, data = pi_2[pi_2$chr == "scaffold_14", ], type = "l")


theta_2 <- read.table(file = "popoolation_sumstats/MN_BIO1_S1/MN_BIO1_S1.window.theta",
                   header = F, sep = "\t", col.names = header4theta)

theta_2[theta_2$nsnps==0, "pi"] <- NA
mean(theta_2$theta, na.rm = T)

mean(theta_2$theta[theta_2$theta < 0.1], na.rm = T)

par(mfrow=c(2, 1))
plot(pi ~ pos, data = pi[pi$chr == "scaffold_14", ], type = "l")
plot(theta ~ pos, data = theta_2[theta_2$chr == "scaffold_14", ], type = "l")


exon_pi <- read.table(file = "popoolation_sumstats/MN_BIO1_S1/MN_BIO1_S1.exons.pi",
                      header = F, sep = "\t", col.names = header4pi[-c(2)])

exon_pi[exon_pi$nsnps==0, "pi"] <- NA
mean(exon_pi$pi, na.rm = T)

mean(theta_2$theta[theta_2$theta < 0.1], na.rm = T)

par(mfrow=c(2, 1))
plot(pi ~ pos, data = pi[pi$chr == "scaffold_14", ], type = "l")
plot(theta ~ pos, data = theta_2[theta_2$chr == "scaffold_14", ], type = "l")

