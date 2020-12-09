ijob -A jcbnunez -c 2 --mem=120G  --partition=largemem
module load gcc/7.1.0  
module load openmpi/3.1.4
module load gdal
module load proj
module load goolf R/4.0.0
module load bcftools
R

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

  priv.dt[,V8:=paste(V9, paste(V6, V7, V8, sep=""), sep=";")]

## Add names to the prov object
names(priv.dt) = c(
					"AF_bin_aob",
					"chr",
					"pos",
					"nPop",
					"missingPop",
					"ref",
					"alt",
					"afs_mut_pop",
					"afs_zero"
					)

## Restrict any further analyses to the 4 main chromosomes				
priv.dt = priv.dt[which(priv.dt$chr %in% c("2L","2R","3L","3R")),]					
									
#Parse the allele frequency data 
 unlist_vect = lapply( priv.dt[,9] , unlist )
 split_vect = lapply( unlist_vect , strsplit, split="\\+" )
 remove_1stobj = lapply( split_vect$afs_zero, function(x) x[-1])
 num_vect = lapply( remove_1stobj , as.numeric )
 mean_AF = lapply( num_vect , mean )
 AF_vector = as.data.frame(do.call("rbind", mean_AF))
 
# Add the allele frequency bin vector to the priv object
	priv.dt %<>% 
	mutate(AF = round(AF_vector$V1,2) ) %>% 
	.[which(.$AF >= maf_thresh & .$AF <= (1-maf_thresh) ),] %>% 
	 mutate(Raster_data = paste(chr,nPop,AF, sep = "_"))

as.data.frame(table(priv.dt$Raster_data  , dnn = list("Raster_data")), responseName = "Freq") %>% 
separate(Raster_data, into = c("chr","nPop","AF"), sep = "_") -> o


o[,c("nPop","AF","Freq")] = sapply(o[,c("nPop","AF","Freq")], as.numeric)

o %<>% mutate(caller = caller[2])

save(o, file="./SNAPE.MAFRASTER.Rdata")


### Plot
load("./SNAPE.MAFRASTER.Rdata")

o %>% 
mutate(AF_fold = ifelse(.$AF > 0.5, (1-.$AF), .$AF)) %>% 
ggplot(aes(x=nPop, y=AF_fold, color = log10(Freq))) + 
geom_point(alpha = 0.8, size = 0.9, shape = 15) + 
theme_bw() + 
scale_color_viridis() -> 
rasterMAF

o %>% 
group_by(nPop,chr) %>% 
summarize(Nloc = sum(Freq)) %>% 
ggplot(aes(x=nPop, y=log10(Nloc), color = chr)) + 
theme_bw() + 
geom_line() -> function_count

ggsave("SNAPE_fig.pdf",function_count/rasterMAF, width = 4, height = 5)
