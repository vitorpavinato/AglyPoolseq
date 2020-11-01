#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### gather
library(data.table)
library(foreach)

fn <- list.files("/scratch/aob2x/dest/geo_endemic", "summarySet", full.names=T)

o <- foreach(fn.i=fn)%do%{
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
             list(set, nPop=V3)]

### mutation classes
  o2.ag <- o[,list(.N), list(mt, nPop=V3)]
  o2.ag <- na.omit(o2.ag)[,list(freq=N/sum(N), mt), list(nPop)]

### allele freq. binning
  maf <- function(x) min(x, 1-x)

  o[,af.bin:=round(af, 2)]
  o[,maf:=sapply(o$af, maf)]
  o[,maf.bin:=round(log10(maf), 1)]

  o3.ag <- o[,list(delta.dist=mean(meanDist[set=="obs"] - meanDist[set=="exp"], na.rm=T)),
             list(nPop=V3, maf.bin=maf.bin)]

### how many
  o4.ag <- o[set=="obs",list(nSites=mean(n)), list(nPops=V3)]

  ### get DGRP mt
    dgrp <- fread("/scratch/aob2x/dest/dgrp2.bim")
    dgrp[,mt:=paste(V5, V6)]

    dgrp <- dgrp[V1<=5][nchar(V5)==1 & nchar(V6)==1]
    dgrp[,mt:=paste(V5, V6, sep="")]
    dgrp <- dgrp[!grepl("0", mt)]
    dgrp.ag <- dgrp[,list(.N), mt]
    dgrp.ag <- dgrp.ag[,list(freq=N/sum(N), mt)]




save(o.ag, o2.ag, o3.ag, o4.ag, dgrp.ag, file="~/allSummarySet_endemism.Rdata")



#### plot
  scp aob2x@rivanna.hpc.virginia.edu:~/allSummarySet_endemism.Rdata ~/.
  #scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/dest/geo_endemic/summarySet1.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(viridis)

  load("~/allSummarySet_endemism.Rdata")
  #load("~/summarySet1.Rdata")

### make expectation MA table
  exp.mt <- data.table(mt=c("AT", "AC", "AG", "TA", "TC", "TG", "CA", "CT", "CG", "GA", "GT", "GC"),
              rate=c(6, 3, 5.5, 6, 5.5, 3., 9, 20, 6.5, 20, 9, 6.5)/100)


### plots

    mutation.plot <- ggplot(data=o2.ag, aes(x=nPop, y=freq, group=mt, color=mt)) + geom_line() + cowplot::theme_cowplot() +
    xlab("Number of polymorphic populations") + ylab("Frequency") +
    geom_hline(data=exp.mt, aes(yintercept=rate, color=mt), linetype="dashed") +
    geom_hline(data=dgrp.ag, aes(yintercept=freq, color=mt), linetype="solid") +
    facet_wrap(~mt)


    dist.plot <- ggplot(data=o.ag, aes(x=nPop, y=mean.dist, group=set, color=set)) +
    geom_ribbon(aes(x=nPop, ymin=mean.dist-sd.dist, ymax=mean.dist+sd.dist, group=set, fill=set), alpha=.5) +
    geom_line(size=1) +
    xlab("Number of polymorphic populations") +
    ylab("Average geographic distance between populations (km") +
    cowplot::theme_cowplot()

    phyloConcord.plot <- ggplot(data=o.ag, aes(x=nPop, y=mean.cc_equal, group=set, color=set)) +
    geom_line() +
    xlab("Number of polymorphic populations") +
    ylab("Probability that all populations \nare in same phylogeographic cluster") +
    cowplot::theme_cowplot()

    n.plot <- ggplot(o4.ag, aes(x=nPops, y=log10(nSites))) +
    geom_line() +
    cowplot::theme_cowplot()

    ggplot(na.omit(o3.ag), aes(x=(nPop), y=maf.bin, fill=delta.dist)) + geom_tile(interpolate=T) +
    scale_fill_viridis()


(dist.plot | phyloConcord.plot | n.plot) / mutation.plot








o[V3==2][set=="obs"]

  ggplot() +
  geom_boxplot(data=o, aes(x=as.factor(V3), y=meanDist, fill=set) )

  +

  xlab("Number of populations polymorphic for SNP")

  pw.dist <- spDists(x=as.matrix(samps[set%in%c("DrosRTEC", "DrosEU"), c("long", "lat"), with=F]), longlat=T)
  diag(pw.dist) <- NA
  oe <- foreach(i=1:dim(priv.dt)[1])%dopar%{

    if(i%%500==0) message(paste(i, dim(priv.dt)[1], sep=" / "))

    median(pw.dist[sample(246,2), sample(246,2)])

  }

  oe <- do.call("c", oe)

  t.test(o$meanDist, oe)

  save(o, oe, file="~/oeo.Rdata")


  load("~/oeo.Rdata")
  library(data.table)
  library(ggplot2)

  dt <- data.table(var=c(rep("obs", dim(o)[1]), rep("exp", length(oe))),
                  dist=c(o$meanDist, oe))


  ggplot(data=dt, aes(x=var, y=dist)) + geom_boxplot()
