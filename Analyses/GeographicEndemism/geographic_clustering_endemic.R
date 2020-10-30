#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


args = commandArgs(trailingOnly=TRUE)
set <- as.numeric(args[1])
#set <- 1

### libraries

#### now R
#### libraries
  library(data.table)
  library(sp)
  library(foreach)

### load sample names
  pops <- names(fread("zcat /scratch/aob2x/dest/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | head -n40", nrows=1, skip="#CHR"))[-(1:9)]

### samps
  samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")
  setkey(samps, sampleId)

### get distances
  pw.dist <- spDists(x=as.matrix(samps[set%in%c("DrosRTEC", "DrosEU"), c("long", "lat"), with=F]), longlat=T)
  diag(pw.dist) <- NA
  pw.dist[lower.tri(pw.dist)] <- NA

### private
  priv.dt <- fread("/scratch/aob2x/dest/geo_endemic/geo_endemic.delim")
  priv.dt <- priv.dt[V3>1]
  priv.dt[,list(.N), list(V3)]

  priv.dt.small <- priv.dt[,list(pops=sample(V7, 20, replace=T), n=.N), list(chr=V1, V3)]
  setkey(priv.dt.small, chr)

  priv.dt.small <- priv.dt.small[J(c("2L", "2R", "3L", "3R", "X"))]

  rm(priv.dt)

  dim(priv.dt.small)[1]

  o <- foreach(i=1:dim(priv.dt.small)[1])%dopar%{
    #i<-1
    if(i%%2==0) message(paste(i, dim(priv.dt.small)[1], sep=" / "))

    ### obs
      set.obs <- samps[J(pops[as.numeric(unlist(tstrsplit(priv.dt.small$pops[i], ";"))[-1])])]

      pw.dist.obs <- spDists(x=as.matrix(set.obs[,c("long", "lat"), with=F]), longlat=T)

    ### exp
      set.exp <- samps[J(sample(pops, length(set.obs)))]

      pw.dist.exp <- spDists(x=as.matrix(set.exp[,c("long", "lat"), with=F]), longlat=T)

      cbind(priv.dt.small[i],
          data.table(set=c("obs", "exp"),
                      meanDist=c(mean(pw.dist.obs[lower.tri(pw.dist.obs)]),
                                mean(pw.dist.exp[lower.tri(pw.dist.exp)])),
                      minDist=c(min(pw.dist.obs[lower.tri(pw.dist.obs)]),
                                min(pw.dist.exp[lower.tri(pw.dist.exp)])),
                      maxDist=c(max(pw.dist.obs[lower.tri(pw.dist.obs)]),
                                  max(pw.dist.exp[lower.tri(pw.dist.exp)])),
                      cc_equal=c(length(set.obs$Continental_clusters)==1,
                                length(set.exp$Continental_clusters)==1)))

  }


  o <- rbindlist(o)

  save(o, file=paste("/scratch/aob2x/dest/geo_endemic/summarySet", set, ".Rdata", sep=""))
