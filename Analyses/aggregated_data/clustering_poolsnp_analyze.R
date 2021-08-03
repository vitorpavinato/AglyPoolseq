##########################################################
#                 POPULATION STRUCTURE                   # 
#       PRNCIPAL COMPONENT ANALYSIS AND CLUSTERING       #
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
library(ggpubr)
library(missMDA)
library(ade4)
library(factoextra)
library(zoo)
library(FactoMineR)
library(Hmisc)

# Import auxiliary R functions  
source("DEST-AglyPoolseq/Analyses/aggregated_data/aux_func.R")

ALPHA=0.75

## POOL'S INFORMATION
# Pool's names
poolnames <- c("MN-Av.1", "MN-V.1",
               "ND-Av.1", "ND-V.1", 
               "NW-Av.1", "NW-Av.2", "NW-V.1", "NW-V.2", "NW-V.3", 
               "PA-Av.1", "PA-V.1", "PA-V.2", 
               "WI-Av.1", "WI-Av.2", "WI-V.1", "WI-V.2", "WI-V.3", 
               "WO-Av.1", "WO-V.1", "WO-V.2", "WO-V.3")

# Pool's/biotypes colors - for PCA
biotype.col <- c("#247F00","#AB1A53",
                 "#247F00","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53",
                 "#247F00","#247F00","#AB1A53","#AB1A53","#AB1A53",
                 "#247F00","#AB1A53","#AB1A53","#AB1A53")

biotype.sym <- c(15,17,
                 15,17,
                 15,15,17,17,17,
                 15,17,17,
                 15,15,17,17,17,
                 15,17,17,17)

## LOAD THE LIST OF COMMON SNP TO USE IN THE PRINCIPAL COMPONENT ANALYSIS
list_common_snps_names <- readLines("results/dataset_comparison/repl-aggr_aggr_common_snps_names.txt")

## LOAD ALLELE FREQUENCY DATA
# DATASET 2: ALLELE FREQUENCY DATA OF TECHNICAL REPLICATES AGGREGATED FROM READ ALIGNMENT
dt.1.imputedRefMLFreq.dataFrame <- read.table(file="results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/imputedML_REFAlleleFrequencies_pools.txt", 
                                              sep = "\t", header = T)
dim(dt.1.imputedRefMLFreq.dataFrame)

## SUBSET EACH DATASET TO KEEP ONLY COMMON SNPS
dt.1.imputedRefMLFreq.subset <-dt.1.imputedRefMLFreq.dataFrame[paste0(dt.1.imputedRefMLFreq.dataFrame$Chromosome, "_",
                                                                      dt.1.imputedRefMLFreq.dataFrame$Position) %in% 
                                                                       list_common_snps_names, ]
dim(dt.1.imputedRefMLFreq.subset)

###
###
### ---- PRINCIPAL COMPONENT ANALYSIS  ----
###
###

## DATASET 1 - WITH ALL POOLS:
W_1 <- scale(t(dt.1.imputedRefMLFreq.subset[,-c(1:4)]), scale=TRUE) #centering
W_1[1:10,1:10]
W_1[is.na(W_1)]<-0
cov.W_1<-cov(t(W_1))
eig.result_1<-eigen(cov.W_1)
eig.vec_1<-eig.result_1$vectors
lambda_1<-eig.result_1$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_1/sum(lambda_1),ylab="Fraction of total variance", ylim=c(0,0.1), type='b', cex=1.1,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_1/sum(lambda_1), col="red")

l1_1 <- 100*lambda_1[1]/sum(lambda_1) 
l2_1 <- 100*lambda_1[2]/sum(lambda_1)

# PCA plot
pdf("results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/pca_imputed_allelefreq_repl_common_snps_pools.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_1[,1],eig.vec_1[,2], col=biotype.col,
     xlim = c(-0.55, 0.45), ylim = c(-0.65, 0.45),
     xlab=paste0("PC1 (", round(l1_1,2), "%)"), ylab=paste0("PC1 (", round(l2_1,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_1[,1],eig.vec_1[,2], 
     poolnames, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.1)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## DATASET 2 - WITHOUT OUTLIERS POOLS:
W_2 <- scale(t(dt.1.imputedRefMLFreq.subset[,-c(1:4, 9,25)]), scale=TRUE) #centering
W_2[1:10,1:10]
W_2[is.na(W_2)]<-0
cov.W_2<-cov(t(W_2))
eig.result_2<-eigen(cov.W_2)
eig.vec_2<-eig.result_2$vectors
lambda_2<-eig.result_2$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_2/sum(lambda_2),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_2/sum(lambda_2), col="red")

l1_2 <- 100*lambda_2[1]/sum(lambda_2) 
l2_2 <- 100*lambda_2[2]/sum(lambda_2)

# PCA plot
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_2[,1],eig.vec_2[,2], col=biotype.col[-c(5,21)],
     xlim = c(-0.7, 0.45), ylim = c(-0.55, 0.65),
     xlab=paste0("PC1 (", round(l1_2,2), "%)"), ylab=paste0("PC1 (", round(l2_2,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_2[,1],eig.vec_2[,2], 
     poolnames[-c(5,21)], pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)


## DATASET 3 - ONE POOL PER BIOTYPE/LOCALITY (LOW % OF MISSING DATA)

# KEEP Biotype pools with less missing data in each locality
biotype.pools2keep <- c(1, 2, 3, 4, 6, 7,
                        10,11,14,16,18,20)

biotype.pools2remove <- setdiff(x=seq(1,21,1), y=biotype.pools2keep)

W_3 <- scale(t(dt.1.imputedRefMLFreq.subset[,-c(1:4,( biotype.pools2remove+4))]), scale=TRUE) #centering
W_3[1:10,1:10]
W_3[is.na(W_3)]<-0
cov.W_3<-cov(t(W_3))
eig.result_3<-eigen(cov.W_3)
eig.vec_3<-eig.result_3$vectors
lambda_3<-eig.result_3$values

# PVE
par(mar=c(5,5,4,1)+.1)
plot(lambda_3/sum(lambda_3),ylab="Fraction of total variance", ylim=c(0,0.2), type='b', cex=1.2,
     cex.lab=1.6, pch=19, col="black")
lines(lambda_3/sum(lambda_3), col="red")

l1_3 <- 100*lambda_3[1]/sum(lambda_3) 
l2_3 <- 100*lambda_3[2]/sum(lambda_3)

# PCA plot
pdf("results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/pca_imputed_allelefreq_repl_common_snps_biotypes.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
par(mar=c(5,5,4,1)+.1)
plot(eig.vec_3[,1],eig.vec_3[,2], col=biotype.col[-c(biotype.pools2remove)],
     xlim = c(-0.7, 0.55), ylim = c(-0.65, 0.55),
     xlab=paste0("PC1 (", round(l1_3,2), "%)"), ylab=paste0("PC1 (", round(l2_3,2), "%)"), 
     cex=1.5, pch=19, cex.lab=1.6)
text(eig.vec_3[,1],eig.vec_3[,2], 
     poolnames[-c(biotype.pools2remove)], pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

###
###
### ---- PRINCIPAL COMPONENT ANALYSIS WITH FactorMineR ----
###
###

## DATASET 1 - WITH ALL POOLS:

# Convert NA cells into loci means
dt.1.dataset1.NAimputed = na.aggregate(t(dt.1.imputedRefMLFreq.subset[,-c(1:4)]))

# PCA
pca.dt.1.dataset1 <- PCA(dt.1.dataset1.NAimputed, scale.unit = F, graph = F)

barplot(pca.dt.1.dataset1$eig[,1],main="Eigenvalues",names.arg=1:nrow(pca.dt.1.dataset1$eig))
summary(pca.dt.1.dataset1)

## Run cluster analysis
# 1. Loading and preparing data
pca.dt.1.dataset1.ind_coord <- pca.dt.1.dataset1$ind$coord

cluster_discovery_dataset1 = fviz_nbclust(pca.dt.1.dataset1.ind_coord, kmeans, method = "gap_stat")
#ggsave("cluster_discovery_dataset1.pdf", cluster_discovery_dataset1,  width =4, height = 4)

# 2. Compute k-means at K=2, K=3 and at K=4
set.seed(123)
kmeans_k2_dataset1 <- kmeans(pca.dt.1.dataset1.ind_coord, 2, nstart = 20)
kmeans_k3_dataset1 <- kmeans(pca.dt.1.dataset1.ind_coord, 3, nstart = 20)
kmeans_k4_dataset1 <- kmeans(pca.dt.1.dataset1.ind_coord, 4, nstart = 20)
kmeans_k5_dataset1 <- kmeans(pca.dt.1.dataset1.ind_coord, 5, nstart = 20)

# Combine Pool ID with cluster assigment
data.frame(k2_clusters = kmeans_k2_dataset1$cluster,
           k3_clusters = kmeans_k3_dataset1$cluster,
           k4_clusters = kmeans_k4_dataset1$cluster,
           k5_clusters = kmeans_k5_dataset1$cluster) %>% mutate(sampleId = poolnames) -> pca.dt.1.dataset1_cluster_data

### COMBINED PCA PLOTS
pdf("results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/pca_imputed_allelefreq_repl_common_snps_factorminer_clusters_pools.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)

par(mar=c(5,5,4,1)+.1, mfrow=c(2, 2))
# PLOT PCA COLORED BY BIOTYPE
plot(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], col=biotype.col,
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=1",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[2,2],2), "%)"), 
     cex=1.7, pch=19)
text(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], 
     pca.dt.1.dataset1_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K2
plot(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], col=pca.dt.1.dataset1_cluster_data$k2_clusters,
     xlim = c(-2.5, 2), ylim = c(-3, 2),  main="K=2",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym)
text(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], 
     pca.dt.1.dataset1_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2", "Avirulent", "Virulent"), col = c(1,2, "grey", "grey"), 
       pch = c(19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K3
plot(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], col=pca.dt.1.dataset1_cluster_data$k3_clusters,
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=3",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym)
text(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], 
     pca.dt.1.dataset1_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2","K3", "Avirulent", "Virulent"), col = c(1,2,3, "grey", "grey"), 
       pch = c(19,19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K4
plot(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], col=pca.dt.1.dataset1_cluster_data$k4_clusters,
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=4",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym)
text(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], 
     pca.dt.1.dataset1_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2","K3","K4", "Avirulent", "Virulent"), col = c(1,2,3,4, "grey", "grey"), 
       pch = c(19,19,19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K5
#plot(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], col=pca.dt.1.dataset1_cluster_data$k5_clusters,
#     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=5",
#     xlab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[1,2],2), "%)"), 
#     ylab=paste0("PC1 (", round(pca.dt.1.dataset1$eig[2,2],2), "%)"), 
#     cex=1.7, pch=biotype.sym)
#text(pca.dt.1.dataset1$svd$U[,1], pca.dt.1.dataset1$svd$U[,2], 
#     pca.dt.1.dataset1_cluster_data$sampleId, pos=2 , cex = 0.6)
#legend("bottomleft", 
#       legend = c("K1", "K2","K3","K4","K5", "Avirulent", "Virulent"), col = c(1,2,3,4,5, "grey", "grey"), 
#       pch = c(19,19,19,19,19,15,17), bty = "n", cex = 1.2)
#abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Savepoint_1
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/clustering.poolsnp.workspace29Jul21.RData.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/clustering.poolsnp.workspace29Jul21.RData.RData")
############## END THIS PART


## DATASET 2 - WITHOUT OUTLIERS POOLS:

# Convert NA cells into loci means
dt.1.dataset2.NAimputed = na.aggregate(t(dt.1.imputedRefMLFreq.subset[,-c(1:4,9,25)]))

# PCA
pca.dt.1.dataset2 <- PCA(dt.1.dataset2.NAimputed, scale.unit = F, graph = F)

barplot(pca.dt.1.dataset2$eig[,1],main="Eigenvalues",names.arg=1:nrow(pca.dt.1.dataset2$eig))
summary(pca.dt.1.dataset2)

## Run cluster analysis
# 1. Loading and preparing data
pca.dt.1.dataset2.ind_coord <- pca.dt.1.dataset2$ind$coord

cluster_discovery_dataset2 = fviz_nbclust(pca.dt.1.dataset2.ind_coord, kmeans, method = "gap_stat")
#ggsave("cluster_discovery_dataset2.pdf", cluster_discovery_dataset2,  width =4, height = 4)

# 2. Compute k-means at K=2, K=3 and at K=4
set.seed(123)
kmeans_k2_dataset2 <- kmeans(pca.dt.1.dataset2.ind_coord, 2, nstart = 20)
kmeans_k3_dataset2 <- kmeans(pca.dt.1.dataset2.ind_coord, 3, nstart = 20)
kmeans_k4_dataset2 <- kmeans(pca.dt.1.dataset2.ind_coord, 4, nstart = 20)
kmeans_k5_dataset2 <- kmeans(pca.dt.1.dataset2.ind_coord, 5, nstart = 20)

# Combine Pool ID with cluster assigment
data.frame(k2_clusters = kmeans_k2_dataset2$cluster,
           k3_clusters = kmeans_k3_dataset2$cluster,
           k4_clusters = kmeans_k4_dataset2$cluster,
           k5_clusters = kmeans_k5_dataset2$cluster) %>% mutate(sampleId = poolnames[-c(5,21)]) -> pca.dt.1.dataset2_cluster_data

### COMBINED PCA PLOTS
pdf("results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/pca_imputed_allelefreq_repl_common_snps_factorminer_clusters_noutliers.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)

par(mar=c(5,5,4,1)+.1, mfrow=c(2, 2))
# PLOT PCA COLORED BY BIOTYPE
plot(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], col=biotype.col[-c(5,21)],
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=1",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[2,2],2), "%)"), 
     cex=1.7, pch=19)
text(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], 
     pca.dt.1.dataset2_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K2
plot(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], col=pca.dt.1.dataset2_cluster_data$k2_clusters,
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=2",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], 
     pca.dt.1.dataset2_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2", "Avirulent", "Virulent"), col = c(1,2, "grey", "grey"), 
       pch = c(19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K3
plot(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], col=pca.dt.1.dataset2_cluster_data$k3_clusters,
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=3",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], 
     pca.dt.1.dataset2_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2","K3", "Avirulent", "Virulent"), col = c(1,2,3, "grey", "grey"), 
       pch = c(19,19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K4
plot(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], col=pca.dt.1.dataset2_cluster_data$k4_clusters,
     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=4",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(5,21)])
text(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], 
     pca.dt.1.dataset2_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2","K3","K4", "Avirulent", "Virulent"), col = c(1,2,3,4, "grey", "grey"), 
       pch = c(19,19,19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K5
#plot(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], col=pca.dt.1.dataset2_cluster_data$k5_clusters,
#     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=5",
#     xlab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[1,2],2), "%)"), 
#     ylab=paste0("PC1 (", round(pca.dt.1.dataset2$eig[2,2],2), "%)"), 
#     cex=1.7, pch=biotype.sym[-c(5,21)])
#text(pca.dt.1.dataset2$svd$U[,1], pca.dt.1.dataset2$svd$U[,2], 
#     pca.dt.1.dataset2_cluster_data$sampleId, pos=2 , cex = 0.6)
#legend("bottomleft", 
#       legend = c("K1", "K2","K3","K4","K5", "Avirulent", "Virulent"), col = c(1,2,3,4,5, "grey", "grey"), 
#       pch = c(19,19,19,19,19,15,17), bty = "n", cex = 1.2)
#abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Savepoint_2
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/clustering.poolsnp.workspace29Jul21.RData.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/clustering.poolsnp.workspace29Jul21.RData.RData")
############## END THIS PART

## DATASET 3 - ONE POOL PER BIOTYPE/LOCALITY (LOW % OF MISSING DATA):

# Convert NA cells into loci means
dt.1.dataset3.NAimputed = na.aggregate(t(dt.1.imputedRefMLFreq.subset[,-c(1:4,(biotype.pools2remove+4))]))

# PCA
pca.dt.1.dataset3 <- PCA(dt.1.dataset3.NAimputed, scale.unit = F, graph = F)

barplot(pca.dt.1.dataset3$eig[,1],main="Eigenvalues",names.arg=1:nrow(pca.dt.1.dataset3$eig))
summary(pca.dt.1.dataset3)

## Run cluster analysis
# 1. Loading and preparing data
pca.dt.1.dataset3.ind_coord <- pca.dt.1.dataset3$ind$coord

cluster_discovery_dataset3 = fviz_nbclust(pca.dt.1.dataset3.ind_coord, kmeans, method = "gap_stat")
#ggsave("cluster_discovery_dataset3.pdf", cluster_discovery_dataset3,  width =4, height = 4)

# 2. Compute k-means at K=2, K=3 and at K=4
set.seed(123)
kmeans_k2_dataset3 <- kmeans(pca.dt.1.dataset3.ind_coord, 2, nstart = 20)
kmeans_k3_dataset3 <- kmeans(pca.dt.1.dataset3.ind_coord, 3, nstart = 20)
kmeans_k4_dataset3 <- kmeans(pca.dt.1.dataset3.ind_coord, 4, nstart = 20)
kmeans_k5_dataset3 <- kmeans(pca.dt.1.dataset3.ind_coord, 5, nstart = 20)

# Combine Pool ID with cluster assigment
data.frame(k2_clusters = kmeans_k2_dataset3$cluster,
           k3_clusters = kmeans_k3_dataset3$cluster,
           k4_clusters = kmeans_k4_dataset3$cluster,
           k5_clusters = kmeans_k5_dataset3$cluster) %>% mutate(sampleId = poolnames[-c(biotype.pools2remove)]) -> pca.dt.1.dataset3_cluster_data

### COMBINED PCA PLOTS
pdf("results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/pca_imputed_allelefreq_repl_common_snps_factorminer_clusters_biotypes.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)

par(mar=c(5,5,4,1)+.1, mfrow=c(2, 2))
# PLOT PCA COLORED BY BIOTYPE
plot(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], col=biotype.col[-c(biotype.pools2remove)],
     xlim = c(-2.5, 2), ylim = c(-2.5, 1.2), main="K=1",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[2,2],2), "%)"), 
     cex=1.7, pch=19, cex.lab=1.6, cex.axis=1.4)
text(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], 
     pca.dt.1.dataset3_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K2
plot(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], col=pca.dt.1.dataset3_cluster_data$k2_clusters,
     xlim = c(-2.5, 2), ylim = c(-2.5, 1.2), main="K=2",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(biotype.pools2remove)], cex.lab=1.6, cex.axis=1.4)
text(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], 
     pca.dt.1.dataset3_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("K1", "K2", "Avirulent", "Virulent"), col = c(1,2, "grey", "grey"), 
       pch = c(19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K3
plot(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], col=pca.dt.1.dataset3_cluster_data$k3_clusters,
     xlim = c(-2.5, 2), ylim = c(-2.5, 1.2), main="K=3",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(biotype.pools2remove)], cex.lab=1.6, cex.axis=1.4)
text(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], 
     pca.dt.1.dataset3_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("K1", "K2","K3", "Avirulent", "Virulent"), col = c(1,2,3, "grey", "grey"), 
       pch = c(19,19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K4
plot(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], col=pca.dt.1.dataset3_cluster_data$k4_clusters,
     xlim = c(-2.5, 2), ylim = c(-2.5, 1.2), main="K=4",
     xlab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[2,2],2), "%)"), 
     cex=1.7, pch=biotype.sym[-c(biotype.pools2remove)], cex.lab=1.6, cex.axis=1.4)
text(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], 
     pca.dt.1.dataset3_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("topleft", 
       legend = c("K1", "K2","K3","K4", "Avirulent", "Virulent"), col = c(1,2,3,4, "grey", "grey"), 
       pch = c(19,19,19,19,15,17), bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)

# PLOT PCA COLORED BY K5
#plot(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], col=pca.dt.1.dataset3_cluster_data$k5_clusters,
#     xlim = c(-2.5, 2), ylim = c(-3, 2), main="K=5",
#     xlab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[1,2],2), "%)"), 
#     ylab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[2,2],2), "%)"), 
#     cex=1.7, pch=biotype.sym[-c(biotype.pools2remove)])
#text(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], 
#     pca.dt.1.dataset3_cluster_data$sampleId, pos=2 , cex = 0.6)
#legend("bottomleft", 
#       legend = c("K1", "K2","K3","K4","K5", "Avirulent", "Virulent"), col = c(1,2,3,4,5, "grey", "grey"), 
#       pch = c(19,19,19,19,19,15,17), bty = "n", cex = 1.2)
#abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Savepoint_3
##-------------
#save.image("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/clustering.poolsnp.workspace29Jul21.RData.RData")
#load("/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/clustering.poolsnp.workspace29Jul21.RData.RData")
############## END THIS PART

### PCA AND CLUSTER DISCOVERY GAP STATISTICS
### ONE POOL PER BIOTYPE/LOCALITY (LOW % OF MISSING DATA): 

### COMBINED PCA PLOTS
pdf("results/aggregated_data/minmaxcov_4_99/clustering_poolsnp/pca_imputed_allelefreq_repl_common_snps_factorminer_clusters_discov_biotypes.pdf",         # File name
    width = 16, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)

par(mar=c(5,5,4,1)+.1, mfrow=c(1, 2))
# PLOT CLUSTER DISCOVERY GAP STATISTICS
errbar(seq(1,10), 
       cluster_discovery_dataset3$data$gap,
       cluster_discovery_dataset3$data$ymax,
       cluster_discovery_dataset3$data$ymin,
       type='b', pch=19, ylim=c(-0.05, 0.60),
       ylab="Gap statistics (k)", xlab="Number of clusters k", cex.lab=1.6, cex.axis=1.4)
abline(v=1, col="#b2182b", lty=3, lwd=2)
text(x=2.8, y=0.6, "Optimal number of clusters", col="#b2182b",cex=1.2)

# PLOT PCA COLORED BY BIOTYPE
plot(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], col=biotype.col[-c(biotype.pools2remove)],
     xlim = c(-2.5, 2), ylim = c(-2.5, 1.2),
     xlab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[1,2],2), "%)"), 
     ylab=paste0("PC1 (", round(pca.dt.1.dataset3$eig[2,2],2), "%)"), 
     cex=1.7, pch=19, cex.lab=1.6, cex.axis=1.4)
text(pca.dt.1.dataset3$svd$U[,1], pca.dt.1.dataset3$svd$U[,2], 
     pca.dt.1.dataset3_cluster_data$sampleId, pos=2 , cex = 0.6)
legend("bottomleft", 
       legend = c("Avirulent", "Virulent"), col = c("#247F00","#AB1A53"), 
       pch = 19, bty = "n", cex = 1.2)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

