ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem
module load gcc/7.1.0  
module load openmpi/3.1.4
module load gdal
module load proj
module load goolf R/4.0.0
module load bcftools
R

#rm(list = ls())
#args = commandArgs(trailingOnly=TRUE)
#set <- as.numeric(args[1])
#caller <- args[2]
caller <- c("PoolSNP", "SNAPE")

#### libraries
  library(data.table)
  library(sp)
  library(foreach)
  library(tidyverse)
  library(magrittr)
  library(patchwork)
  library(viridis)

### load maf filters
  if(caller=="SNAPE") {
    maf_thresh <- 0.05
  } else if(caller=="PoolSNP") {
    maf_thresh <- 0.001
  }

### load sample names
  pops <- names(fread("bcftools view -S /scratch/yey2sn/DEST_endemism/good.samps.txt /project/berglandlab/DEST/vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | head -n 40", nrows=1, skip="#CHR"))[-(1:9)]

##### if using pre Oct 2020 data, do this. Otherwise, don't
#    pops <- gsub("AT_gr_12", "newSpain", pops)
#    pops <- gsub("ES_ba_12", "newAust", pops)
#    pops <- gsub("newSpain", "ES_ba_12", pops)
#    pops <- gsub("newAust", "AT_gr_12", pops)

### samps
  samps <- fread("/scratch/yey2sn/DEST/populationInfo/samps.csv")
  setkey(samps, sampleId)

#add clusters for AT/ES Pops
samps$Continental_clusters[grep("ES_ba_12", samps$sampleId)] = "1.Europe_W"
samps$Continental_clusters[grep("AT_gr_12", samps$sampleId)] = "3.Europe_E"


### get distances
  pw.dist <- spDists(x=as.matrix(samps[set%in%c("DrosRTEC", "DrosEU"), c("long", "lat"), with=F]), longlat=T)
  diag(pw.dist) <- NA
  pw.dist[lower.tri(pw.dist)] <- NA

### private
### caller[1] is poolsnp -- unless overwriten by user
### caller[2] is SANPE -- unless overwriten by user
  fn <- list.files("/scratch/yey2sn/DEST_endemism/", paste("geo_endemic.", caller[2], ".goodSamps.delim", sep=""), full.name=T)


  priv.dt <- foreach(fn.i=fn)%do%{
    message(fn.i)
    fread(fn.i)
  }
  
  priv.dt <- rbindlist(priv.dt)

######
# SNAPE portion
## Add names to the prov object

## build metadata object
meata_data = left_join(data.frame(sampleId=pops), samps) %>% separate(Continental_clusters, into = c("Cluster_code","Cluster_name"), sep = "\\.")

## Add names to columns
names(priv.dt) = c(
					"AF_bin_aob",
					"chr",
					"pos",
					"nPop",
					"missingPop",
					"ref",
					"alt",
					"Pop_codes",
					"afs_zero"
					)

## Restrict any further analyses to the 4 main chromosomes				
priv.dt = priv.dt[which(priv.dt$chr %in% c("2L","2R","3L","3R")),]					

#Parse the allele phylogenetic concrodance data
 unlist_vect = lapply( priv.dt[,8] , unlist )
 split_vect = lapply( unlist_vect , strsplit, split="\\;" )
 remove_1stobj = lapply( split_vect$Pop_codes, function(x) x[-1])
 num_ids = lapply( remove_1stobj , as.numeric )
 extract_popcodes = lapply( num_ids , function(x) meata_data$Cluster_code[x]  )
 num_clusters = lapply( extract_popcodes , as.numeric )
##It turns out that one can use "function(x) abs(ceiling(var(x))-1)" to estimate the probability of concrodance
 phylo_concordance = lapply( num_clusters , function(x) abs(ceiling(var(x))-1) )
 Prob_CC = as.data.frame(do.call("rbind", phylo_concordance))


# Add the allele frequency bin vector to the priv object
	priv.dt %<>% 
	mutate(Prob_CC = Prob_CC$V1) %>% 
	mutate(Prob_CC_data = paste(chr,nPop,Prob_CC, sep = "_"))

as.data.frame(table(priv.dt$Prob_CC_data  , dnn = list("Prob_CC")), responseName = "Freq") %>% 
separate(Prob_CC, into = c("chr","nPop","P_cc"), sep = "_") -> o

o[,c("nPop","P_cc","Freq")] = sapply(o[,c("nPop","P_cc","Freq")], as.numeric)

o_dcast =  o %>% .[which(.$nPop > 1),]  %>% dcast(chr+nPop~P_cc, value.var = "Freq")

o_dcast[is.na(o_dcast)] <- 0

o_dcast %<>% mutate(caller = caller[2], Total_obs = `0`+`1`) %>% mutate(Probality_obs = `1`/Total_obs)

# Generate Expected Distribution
expected_o = list()

N=1000
for(j in 1:N){
sim_out = data.frame()
for(i in 2:75){


sim_out[i,"nPop"] = i
sim_out[i,"sim"] = j
sim_out[i,"success"] = abs(ceiling(var(sample(1:3, i, replace = T)))-1)

expected_o[[j]] = sim_out

} #i
} #j

expected_o = as.data.frame(do.call("rbind", expected_o))
expected_o %<>% group_by(nPop) %>% summarize(Freq = sum(success, na.rm = T)) %>% mutate(Probality_exp = Freq/N ) 

## Save
save(expected_o, o_dcast, file="./SNAPE.Probability_CC.Rdata")


### Plot
load("./SNAPE.Probability_CC.Rdata")

rbind(
melt(id = "nPop", o_dcast[,c("nPop","Probality_obs")]),
melt(id = "nPop", expected_o[,c("nPop","Probality_exp")])
	) %>% separate(variable, into = c("stat","type"), sep = "_" ) %>%
ggplot(aes(x=nPop, y=value, color=type)) + 
geom_smooth(span = 1/10) +
xlim(0,75) +
theme_bw() ->
prob_graph

ggsave("prob_SNAPE_graph.pdf",prob_graph, width = 4, height = 5)
