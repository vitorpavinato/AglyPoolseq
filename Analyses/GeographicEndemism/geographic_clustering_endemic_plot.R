### libraries
  library(data.table)
  library(foreach)

### load summary data
  fn <- list.files("/scratch/aob2x/dest/geo_endemic/maf/", "SNAPE.SummarySet", full.names=T)

  o.ag.all <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    return(o.ag)
  }

  o2.ag.all <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    o2.ag[,maf:=maf]
    return(o2.ag)
  }


  o4.ag.all <- foreach(fn.i=fn, .combine="rbind")%do%{
    load(fn.i)
    maf <- gsub(".endemism.Rdata", "", gsub("/scratch/aob2x/dest/geo_endemic/maf//SNAPE.SummarySet.", "", fn.i))
    o4.ag[,maf:=maf]
    return(o4.ag)
  }

  save(o.ag.all, o2.ag.all, o4.ag.all, dgrp.ag, file="~/allSummarySet_endemism.SNAPE.maf.Rdata")


#### plot
  scp aob2x@rivanna.hpc.virginia.edu:~/allSummarySet_endemism.SNAPE.maf.Rdata ~/.
  #scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/dest/geo_endemic/summarySet1.Rdata ~/.



##### compares MAF thresholds for SNAPE
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(viridis)
  library(ggridges)

  load("~/allSummarySet_endemism.SNAPE.maf.Rdata")
  o.ag.ag <- o.ag.all[,list(mean.dist=mean(mean.dist), mean.cc_equal=mean(mean.cc_equal), sd.dist=mean(sd.dist)), list(set, caller, nPop, maf)]
  o2.ag.ag <- o2.ag.all[,list(freq=mean(freq)), list(nPop, caller, mt, maf)]
  o4.ag.ag <- o4.ag.all[,list(nSites=sum(nSites)), list(nPop, caller, maf)]

  setkey(o2.ag.ag, nPop, caller, maf)
  setkey(o4.ag.ag, nPop, caller, maf)
  o2.ag.ag <- merge(o2.ag.ag, o4.ag.ag)

  exp.mt <- data.table(mt=c("AT", "AC", "AG", "TA", "TC", "TG", "CA", "CT", "CG", "GA", "GT", "GC"),
              rate=c(6, 3, 5.5, 6, 5.5, 3., 9, 20, 6.5, 20, 9, 6.5)/100)


### n plpt
  n.plot <- ggplot(o4.ag.ag, aes(x=(nPop), y=log10(nSites), group=interaction(caller, maf), color=as.factor(maf), linetype=caller)) +
  geom_line() +
  cowplot::theme_cowplot()


  dist.plot <- ggplot(data=o.ag.ag[nPop<50], aes(x=nPop, y=mean.dist, group=interaction(set, caller), color=set, linetype=caller)) +
  geom_ribbon(aes(x=nPop, ymin=mean.dist-sd.dist, ymax=mean.dist+sd.dist, group=set, fill=set), alpha=.5) +
  geom_line(size=1) +
  facet_wrap(.~maf) +
  xlab("Number of polymorphic populations") +
  ylab("Ave. dist. (km)") +
  cowplot::theme_cowplot()

  mutation.plot <- ggplot(data=o2.ag.ag[nSites>500], aes(x=nPop, y=freq, group=mt, color=mt, linetype=caller)) + geom_line() + cowplot::theme_cowplot() +
  xlab("Number of polymorphic populations") + ylab("Frequency") +
  geom_hline(data=exp.mt, aes(yintercept=rate, color=mt), linetype="dashed") +
  geom_hline(data=dgrp.ag, aes(yintercept=freq, color=mt), linetype="solid") +
  facet_grid(maf~mt) +
  theme(legend.position = "none")

  phyloConcord.plot <- ggplot(data=o.ag.ag[nPop<50],
    aes(x=nPop, y=(mean.cc_equal), group=interaction(set, caller), color=set, linetype=caller)) +
  geom_line() +
  facet_wrap(~maf+caller) +
  xlab("Number of polymorphic populations") +
  ylab("Prob same phylo. clust.") +
  cowplot::theme_cowplot()

design <- "
AAAAAAA
AAAAAAA
BBBBBBB
BBBBBBB
CCCC###
CCCC###
"



o.plot <- (n.plot + dist.plot + phyloConcord.plot + phyloConcord.plot)  + plot_layout(design = design)
#o.plot

ggsave(o.plot, file="~/geo_endemic.SNAPE.maf.pdf", h=8, w=11)


  ggsave(mutation.plot, file="~/mutation_plot.SNAPE.maf.pdf", h=8, w=11)



##### compares SNAPE & PoolSNP (I need to remake the input file)
### libraries
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(viridis)

  load("~/allSummarySet_endemism.Rdata")
  o.ag.ag <- o.ag[,list(mean.dist=mean(mean.dist), mean.cc_equal=mean(mean.cc_equal), sd.dist=mean(sd.dist)), list(set, caller, nPop)]
  o2.ag.ag <- o2.ag[,list(freq=mean(freq)), list(nPop, caller, mt)]
  o4.ag.ag <- o4.ag[,list(nSites=sum(nSites)), list(nPop, caller)]


  setkey(o2.ag.ag, nPop, caller)
  setkey(o4.ag.ag, nPop, caller)
  o2.ag.ag <- merge(o2.ag.ag, o4.ag.ag)
  #load("~/summarySet1.Rdata")

### make expectation MA table and do a bit of rearragment of MT data
  exp.mt <- data.table(mt=c("AT", "AC", "AG", "TA", "TC", "TG", "CA", "CT", "CG", "GA", "GT", "GC"),
              rate=c(6, 3, 5.5, 6, 5.5, 3., 9, 20, 6.5, 20, 9, 6.5)/100)


### mutation plot
    mutation.plot <- ggplot(data=o2.ag.ag[nSites>500], aes(x=nPop, y=freq, group=mt, color=mt)) + geom_line() + cowplot::theme_cowplot() +
    xlab("Number of polymorphic populations") + ylab("Frequency") +
    geom_hline(data=exp.mt, aes(yintercept=rate, color=mt), linetype="dashed") +
    geom_hline(data=dgrp.ag, aes(yintercept=freq, color=mt), linetype="solid") +
    facet_grid(caller~mt) +
    theme(legend.position = "none")

### total number plot
    n.plot <- ggplot(o4.ag, aes(x=(nPop), y=log10(nSites), fill=caller)) +
    geom_col(position = "dodge", size=2) +
    cowplot::theme_cowplot()

### distance plot
    dist.plot <- ggplot(data=o.ag.ag[nPop<50], aes(x=nPop, y=mean.dist, group=interaction(set), color=set)) +
    geom_ribbon(aes(x=nPop, ymin=mean.dist-sd.dist, ymax=mean.dist+sd.dist, group=set, fill=set), alpha=.5) +
    geom_line(size=1) +
    facet_grid(.~caller) +
    xlab("Number of polymorphic populations") +
    ylab("Ave. dist. (km)") +
    cowplot::theme_cowplot()

### phylo concord plot
    phyloConcord.plot <- ggplot(data=o.ag.ag[nPop<50],
      aes(x=nPop, y=(mean.cc_equal), group=interaction(set, caller), color=set, linetype=caller)) +
    geom_line() +
    #facet_grid(~caller) +
    xlab("Number of polymorphic populations") +
    ylab("Prob same phylo. clust.") +
    cowplot::theme_cowplot()



  design <- "
  AAAAAA
  AAAAAA
  BBBBBB
  BBBBBB
  CCC###
  CCC###
  "



  o.plot <- (n.plot + dist.plot + phyloConcord.plot + phyloConcord.plot)  + plot_layout(design = design)
  #o.plot

  ggsave(o.plot, file="~/geo_endemic.pdf", h=8, w=11)
  ggsave(mutation.plot, file="~/mutation_plot.pdf", h=8, w=11)






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
