### libraries
  library(data.table)
  library(foreach)

### load summary data for good samps SNAPE
  #fn <- list.files("SummarySet.goodSamps.endemism.Rdata", full.names=T)
  fn <- "/scratch/aob2x/dest/geo_endemic/goodSamps/SummarySet.goodSamps.endemism.Rdata"

  o.ag.snape <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    return(o.ag[caller=="SNAPE"])
  }

  o2.ag.snape <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    #maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    #o2.ag[,maf:=maf]
    return(o2.ag[caller=="SNAPE"])
  }

  o3.ag.snape <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    #maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    #o4.ag[,maf:=maf]
    return(o3.ag[caller=="SNAPE"])
  }


  o4.ag.snape <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    #maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    #o4.ag[,maf:=maf]
    return(o4.ag[caller=="SNAPE"])
  }



### load summary data for PoolSeq PoolSNP
  #fn <- list.files("SummarySet.goodSamps.endemism.Rdata", full.names=T)
  fn <- "/scratch/aob2x/dest/geo_endemic/goodSamps/SummarySet.PoolSeq.PoolSNP.endemism.Rdata"

  o.ag.PoolSeq <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    return(o.ag)
  }

  o2.ag.PoolSeq <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    #maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    #o2.ag[,maf:=maf]
    return(o2.ag)
  }

  o3.ag.PoolSeq <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    #maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    #o4.ag[,maf:=maf]
    return(o3.ag)
  }


  o4.ag.PoolSeq <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    #maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    #o4.ag[,maf:=maf]
    return(o4.ag)
  }


### combine
  o.ag.all <- rbind(o.ag.snape, o.ag.PoolSeq)
  o2.ag.all <- rbind(o2.ag.snape, o2.ag.PoolSeq)
  o4.ag.all <- rbind(o4.ag.snape, o4.ag.PoolSeq)
  o3.ag.all <- rbind(o3.ag.snape, o3.ag.PoolSeq)

  save(o.ag.all, o2.ag.all, o3.ag.all, o4.ag.all, dgrp.ag, file="~/allSummarySet_endemism.bothCallers.Rdata")
  save(o.ag.all, o2.ag.all, o3.ag.all, o4.ag.all, dgrp.ag, file="/project/berglandlab/alan/allSummarySet_endemism.bothCallers.Rdata")
