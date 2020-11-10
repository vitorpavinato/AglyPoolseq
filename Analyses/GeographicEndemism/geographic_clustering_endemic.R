#module load bcftools gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


args = commandArgs(trailingOnly=TRUE)
set <- as.numeric(args[1])
caller <- args[2]

#set <- 132; caller <- "SNAPE"

### libraries

#### now R
#### libraries
  library(data.table)
  library(sp)
  library(foreach)

### load maf filters
  #maf <- unique(fread("/scratch/aob2x/dest/geo_endemic/jobs.txt", header=F)$V1)[set%%10+1]
  if(caller=="SNAPE") {
    maf <- 0.05
  } else if(caller=="PoolSNP") {
    maf <- 0.001
  }

### load sample names
  #pops <- names(fread("zcat /scratch/aob2x/dest/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | head -n40", nrows=1, skip="#CHR"))[-(1:9)]

  pops <- names(fread("bcftools view -S /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/goodSamps.delim /scratch/aob2x/dest/dest.PoolSeq.PoolSNP.001.50.ann.vcf.gz | head -n 40",
        nrows=1, skip="#CHR"))[-(1:9)]


  #### if using pre Oct 2020 data, do this. Otherwise, don't
    pops <- gsub("AT_gr_12", "newSpain", pops)
    pops <- gsub("ES_ba_12", "newAust", pops)
    pops <- gsub("newSpain", "ES_ba_12", pops)
    pops <- gsub("newAust", "AT_gr_12", pops)
  
### samps
  samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")
  setkey(samps, sampleId)

### get distances
  pw.dist <- spDists(x=as.matrix(samps[set%in%c("DrosRTEC", "DrosEU"), c("long", "lat"), with=F]), longlat=T)
  diag(pw.dist) <- NA
  pw.dist[lower.tri(pw.dist)] <- NA

### private
  #fn <- list.files("/scratch/aob2x/dest/geo_endemic", paste("geo_endemic.", caller, ".", maf, sep=""), full.name=T)
  fn <- list.files("/scratch/aob2x/dest/geo_endemic", paste("geo_endemic.", caller, ".goodSamps", sep=""), full.name=T)

  priv.dt <- foreach(fn.i=fn)%do%{
    message(fn.i)
    fread(fn.i)
  }
  priv.dt <- rbindlist(priv.dt)

  #priv.dt <- priv.dt[V3>1]

  #priv.dt[,list(.N), list(V4)]


  priv.dt[,V8:=paste(V9, paste(V6, V7, V8, sep=""), sep=";")]


  priv.dt.small <- priv.dt[,list(pops=sample(V8, 100, replace=T), n=.N), list(chr=V2, nPop=V4)]
  setkey(priv.dt.small, chr)

  priv.dt.small <- priv.dt.small[J(c("2L", "2R", "3L", "3R", "X"))]
  rm(priv.dt)
  priv.dt.small <- na.omit(priv.dt.small)
  dim(priv.dt.small)[1]

  getmean <- function(x) {
    # x <- priv.dt.small.tmp$pops
    af_string <- tstrsplit(x, ";")[[1]]
    mean(as.numeric(unlist(tstrsplit(af_string, "\\+")[-1])))
  }

  priv.dt.small[,id:=c(1:dim(priv.dt.small)[1])]

  o <- foreach(i=priv.dt.small$id[1:10000])%dopar%{
    #i<-201
    if(i%%2==0) message(paste(i, dim(priv.dt.small)[1], sep=" / "))

    ### get average allele frequency
      priv.dt.small.tmp <- priv.dt.small[id==i]
      priv.dt.small.tmp[,af:=getmean(pops)]

    ### obs
      set.obs <- samps[J(pops[as.numeric(unlist(tstrsplit(priv.dt.small.tmp$pops, ";")[-c(1,2)]))])]

      pw.dist.obs <- spDists(x=as.matrix(set.obs[,c("long", "lat"), with=F]), longlat=T)

    ### exp
      set.exp <- samps[J(sample(pops, nrow(set.obs)))]

      pw.dist.exp <- spDists(x=as.matrix(set.exp[,c("long", "lat"), with=F]), longlat=T)


      o <- cbind(rbind(priv.dt.small.tmp, data.table(chr=priv.dt.small.tmp$chr, nPop=priv.dt.small.tmp$nPop, pops=NA, n=NA, id=NA), fill=T),
          data.table(set=c("obs", "exp"),
                      meanDist=c(mean(pw.dist.obs[lower.tri(pw.dist.obs)]),
                                mean(pw.dist.exp[lower.tri(pw.dist.exp)])),

                      minDist=c(min(pw.dist.obs[lower.tri(pw.dist.obs)]),
                                min(pw.dist.exp[lower.tri(pw.dist.exp)])),

                      maxDist=c(max(pw.dist.obs[lower.tri(pw.dist.obs)]),
                                  max(pw.dist.exp[lower.tri(pw.dist.exp)])),

                      cc_equal=c(length(unique(set.obs$Continental_clusters))==1,
                                length(unique(set.exp$Continental_clusters))==1)))
      o[,mt:=tstrsplit(pops, ";")[[2]]]
      o[,caller:=caller]
      o[,maf:=maf]
      return(o)
  }


  o <- rbindlist(o)

  save(o, file=paste("/scratch/aob2x/dest/geo_endemic/goodSamps/summarySet.", set, ".", caller, ".", maf, ".Rdata", sep=""))
