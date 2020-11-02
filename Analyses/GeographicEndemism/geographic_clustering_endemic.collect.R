#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### gather
library(data.table)
library(foreach)

fn <- list.files("/scratch/aob2x/dest/geo_endemic", "summarySet", full.names=T)

o <- foreach(fn.i=fn[1:20])%do%{
  message(fn.i)
  #fn.i <- fn[1]
  load(fn.i)
  return(o)
}
o <- rbindlist(o, fill=T)


### geogrpahic distances
  o.ag <- o[,list(mean.dist=mean(meanDist), sd.dist=sd(meanDist),
                  med.dist=median(meanDist), lci.dist=quantile(meanDist, .025), uci.dist=quantile(meanDist, .975),
                  mean.cc_equal=mean(cc_equal)),
             list(set, nPop, caller)]

### mutation classes
  o2.ag <- o[,list(.N), list(mt, nPop, caller)]
  o2.ag <- na.omit(o2.ag)[,list(freq=N/sum(N), mt), list(nPop, caller)]

### allele freq. binning
  maf <- function(x) min(x, 1-x)

  o[,af.bin:=round(af, 2)]
  o[,maf:=sapply(o$af, maf)]
  o[,maf.bin:=round(log10(maf), 1)]

  o3.ag <- o[,list(delta.dist=mean(meanDist[set=="obs"] - meanDist[set=="exp"], na.rm=T)),
             list(nPop, maf.bin=maf.bin, caller)]

### how many
  o4.ag <- o[set=="obs",list(nSites=mean(n)), list(nPop, caller)]

  ### get DGRP mt
    dgrp <- fread("/scratch/aob2x/dest/dgrp2.bim")
    dgrp[,mt:=paste(V5, V6)]

    dgrp <- dgrp[V1<=5][nchar(V5)==1 & nchar(V6)==1]
    dgrp[,mt:=paste(V5, V6, sep="")]
    dgrp <- dgrp[!grepl("0", mt)]
    dgrp.ag <- dgrp[,list(.N), mt]
    dgrp.ag <- dgrp.ag[,list(freq=N/sum(N), mt)]




save(o.ag, o2.ag, o3.ag, o4.ag, dgrp.ag, file="~/allSummarySet_endemism.Rdata")
