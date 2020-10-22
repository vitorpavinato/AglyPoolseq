library(adegenet)
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(vcfR)
library(patchwork)
library(reshape2)
library(zoo)
library(matrixStats)
library(data.table)
library(SeqArray)
library(LEA)
library(gdsfmt)
library(SNPRelate)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)

source("./ThinLDinR_SNPtable.R")

##################################
# Part 1 -- Prepare Files
##################################

### open GDS file
  genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.ann.gds")

### get target populations
  samps <- fread("/scratch/yey2sn/DEST/samps.csv")

  samps <- rbind(samps[set=="DrosRTEC"],
                samps[set=="DrosEU"],
                samps[set=="dgn"]
                )
 
### get subsample of data to work on
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=samps$sampleId)

  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

## choose number of alleles
 snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
  seqSetFilter(genofile, sample.id=samps$sampleId,
              snps.dt[chr%in%c("2L", "2R", "3L", "3R")][missing<.05][af>.2]$variant.id)

### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")

  dat <- ad$data/dp
  dim(dat)  
  
## Add metadata
    colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
  
  rownames(dat) <- seqGetData(genofile, "sample.id")

# Remove unwated samples -- Because 

samples_to_remove = c(
                      "AT_gr_12_fall",
                      "AT_gr_12_spring",
                      "ES_ba_12_fall",
                      "ES_ba_12_spring",
                      "SIM",
                      "B",
                      "T"
                      )

dat_filt = dat[-which(rownames(dat) %in% samples_to_remove),]

left_join(data.frame(sampleId=rownames(dat_filt)), as.data.frame(samps)) -> DEST_DGN_metadata

# Remove Physical Linkage

snp_info = colnames(dat_filt) 

data.frame(t(dat_filt)) -> dat_filt_t

dat_filt_t %<>% mutate(SNPid=rownames(.)) %>% separate(SNPid, into = c("chr","pos","id"))

rownames(dat_filt_t) = snp_info

dat_filt_t$chr = as.character(dat_filt_t$chr) 
dat_filt_t$pos = as.numeric(dat_filt_t$pos)

picksnps_500<- pickSNPs(dat_filt_t,dist=500)

dat_filt_t_LD500 = dat_filt_t[picksnps_500,] 

dat_filt_t_LD500[,-which(names(dat_filt_t_LD500) %in% c("chr","pos","id") )] %>% t() -> dat_filt_LD500

#Convert NA cells into loci means
dat_filt_maf_LD500_naimp = na.aggregate(dat_filt_LD500)

save(dat_filt_maf_LD500_naimp,DEST_DGN_metadata, file="./DEST_DGN_AllSNPs_Metadata.Rdata")

##################################
# Part 2 -- PCA Analysis
##################################

load("./DEST_DGN_AllSNPs_Metadata.Rdata")

dat_filt_maf_LD500_naimp %>% PCA(scale.unit = F, graph = F, ncp = 50) -> LD500_naimp_PCA_50PCs_object

LD500_naimp_PCA_50PCs_object$ind$coord %>% as.data.frame() %>% mutate(sampleId = rownames(.)) %>% left_join(., DEST_DGN_metadata) -> PCA_coords_metadata

PC1_PVE=round(LD500_naimp_PCA_50PCs_object$eig[1,2],2)
PC2_PVE=round(LD500_naimp_PCA_50PCs_object$eig[2,2],2)
PC3_PVE=round(LD500_naimp_PCA_50PCs_object$eig[3,2],2)

#Graph
label.df_PCA <- PCA_coords_metadata %>% 
  group_by(country) %>% 
  summarize(Dim.1 = mean(Dim.1), Dim.2 = mean(Dim.2), Dim.3 = mean(Dim.3), Dim.4 = mean(Dim.4) ) 

ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + geom_point(data =PCA_coords_metadata ,  aes(x=Dim.1, y=Dim.2, fill = as.factor(country), shape = as.factor(Continental_clusters)),size = 3,  alpha = 0.8)  +  ggrepel::geom_label_repel(data = label.df_PCA, size = 1.7 ,  alpha = 0.7, aes(x=Dim.1, y=Dim.2, ,label = country))  + theme_classic() + theme(legend.position = "none") + xlab("PC 1 (24% VE)") + ylab("PC 2 (9% VE)") + scale_shape_manual(values = c(21,22,23,24)) -> DEST_Woldwide_PCA12

ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + geom_point(data =PCA_coords_metadata ,  aes(x=Dim.1, y=Dim.3, fill = as.factor(country), shape = as.factor(Continental_clusters)),size = 3,  alpha = 0.8)  +  ggrepel::geom_label_repel(data = label.df_PCA, size = 1.7 ,  alpha = 0.7, aes(x=Dim.1, y=Dim.3, ,label = country))  + theme_classic() + theme(legend.position = "none") + xlab("PC 1 (24% VE)") + ylab("PC 3 (3.9% VE)") + scale_shape_manual(values = c(21,22,23,24))  -> DEST_Woldwide_PCA13


#Graph maps
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-125.15, 45.00), ylim = c(-43.00, 65.00), expand = FALSE) + theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + geom_point(data =PCA_coords_metadata , aes(x=as.numeric(long), y=as.numeric(lat), fill = (Dim.1), shape = as.factor(Continental_clusters)),size = 2.5,  alpha = 0.7)  +
  scale_fill_gradient(low = "steelblue", high = "firebrick", na.value = NA) + xlab("Lon") + ylab("Lat") + ggtitle("PC 1 projections onto World Map") + theme(legend.position = "none") + scale_shape_manual(values = c(21,22,23,24)) -> Wolrd_PC1

#ggsave("Wolrd_PC1.pdf",(Wolrd_PC1),  width = 6, height = 4)

####PC 2

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim =  c(-12, 41.00), ylim = c(32.00, 63.00), expand = FALSE) + theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + geom_point(data =PCA_coords_metadata , aes(x=as.numeric(long), y=as.numeric(lat), fill = (Dim.3), shape = as.factor(Continental_clusters)),size = 2.5,  alpha = 0.7)  +
  scale_fill_gradient2(low = "steelblue", high = "darkblue", na.value = NA) + xlab("Lon") + ylab("Lat") + ggtitle("PC 3 projections onto World Map") + theme(legend.position = "none") + scale_shape_manual(values = c(21,22,23,24)) -> Wolrd_PC3

#ggsave("Wolrd_PC2.pdf",(Wolrd_PC2),  width = 6, height = 4)

ggsave("PCA_figure.pdf",( (DEST_Woldwide_PCA12+DEST_Woldwide_PCA13)/(Wolrd_PC1+Wolrd_PC3) ),  width = 8, height = 6)

################################## ###################### ######################
# Part 3 -- Train DAPC to discover demography informative eigenvectors 
################################## ###################### ######################

#load("./DEST_DGN_AllSNPs_Metadata.Rdata")

################################## ###################### ######################
# Train Continental Cluster Model 
################################## ###################### ######################

left_join(data.frame(sampleId=rownames(dat_filt_maf_LD500_naimp)),
          DEST_DGN_metadata[,c("sampleId","Continental_clusters")]
      ) -> dat_filt_maf_LD500_naimp_labels

rownames(dat_filt_maf_LD500_naimp) == dat_filt_maf_LD500_naimp_labels$sampleId

#Checkpoint
save(
	dat_filt_maf_LD500_naimp,
	dat_filt_maf_LD500_naimp_labels, 
	file = "./data_to_train_Cluster_DAPC.Rdata"
	)

# Train model on continental cluster
xval_cluster <- xvalDapc(dat_filt_maf_LD500_naimp, as.factor(dat_filt_maf_LD500_naimp_labels$Intracontinent_cluster), n.pca.max = 300, training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = FALSE)

#save(xval_cluster, file = "./xval_DAPC_optimization.Rdata")
#load("./xval_DAPC_optimization.Rdata")


xval$`Cross-Validation Results` %>% ggplot(aes(x=n.pca, y=success))  + geom_point(shape = 21, size = 3, fill = "grey") + geom_vline(xintercept = 40, linetype = "dashed") + xlab("Number of PCs retained") + ylab("Probability of Successul Prediction") + theme_classic() -> cross_val_results
ggsave("cross_val_results.pdf",cross_val_results,  width =4, height = 4)

################################## ###################### ######################
# Train model on country/state model
################################## ###################### ######################

# Make some corrections to the prior label
left_join(data.frame(sampleId=rownames(dat_filt_maf_LD500_naimp)),
          samps[,c("sampleId","locality")]
      ) -> dat_filt_maf_LD500_naimp_States


dat_filt_maf_LD500_naimp_States[c(48:51,53),"locality"] = "ETHI"
dat_filt_maf_LD500_naimp_States[c(89),"locality"] = "GABO"
dat_filt_maf_LD500_naimp_States[c(97),"locality"] = "NY"
dat_filt_maf_LD500_naimp_States[c(144),"locality"] = "NC"
dat_filt_maf_LD500_naimp_States[c(159,161,165),"locality"] = "SOAF"
dat_filt_maf_LD500_naimp_States[c(239),"locality"] = "UGA"
dat_filt_maf_LD500_naimp_States[c(255),"locality"] = "CA"

dat_filt_maf_LD500_naimp_States %<>% separate(locality, into = c("State","City"), sep = "_")

dat_filt_maf_LD500_naimp_States$State %>% table %>% as.data.frame() -> State_Count
names(State_Count)[1]="State"

# We have to reduce the data set to popualtions with N >= 3 to train the model.
State_Count[which(as.numeric(State_Count$Freq) >= 3 ),] %>% .$State %>% as.character -> StatesToKeep

dat_filt_maf_LD500_naimp_States %>% .[which(.$State %in% StatesToKeep), "sampleId"] -> Samples_to_Keep

dat_filt_maf_LD500_naimp %>% .[which(rownames(.) %in%  Samples_to_Keep),] -> dat_filt_maf_LD500_naimp_forStateTrain

dat_filt_maf_LD500_naimp_States %>% .[which(.$sampleId %in%  Samples_to_Keep),] -> dat_filt_maf_LD500_naimp_Labels_StatesTrain

rownames(dat_filt_maf_LD500_naimp_forStateTrain) == dat_filt_maf_LD500_naimp_Labels_StatesTrain$sampleId

#Checkpoint
save(
	dat_filt_maf_LD500_naimp_forStateTrain,
	dat_filt_maf_LD500_naimp_Labels_StatesTrain, 
	file = "./data_to_train_State_DAPC.Rdata"
	)

xval_state <- xvalDapc(dat_filt_maf_LD500_naimp_forStateTrain, as.factor(dat_filt_maf_LD500_naimp_Labels_StatesTrain$State), n.pca.max = 300, training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = FALSE)

