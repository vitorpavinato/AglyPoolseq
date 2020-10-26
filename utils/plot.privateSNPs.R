# module load intel/18.0 intelmpi/18.0 R/3.6.3; R



### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### fl
  fl <- list.files("/scratch/aob2x/dest/privateSNPs", "pnps.Rdata", full.name=T)

### load
  priv.dt <- foreach(fl.i=fl)%dopar%{
    #fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    m[is.na(af), af:=ac/dp]

    m.ag <- m[,list(n=length(ac),
                    mean_ac=mean(ac, na.rm=T),
                    mean_dp=mean(dp, na.rm=T),
                    mean_af=mean(af, na.rm=T)),
                list(class, ann, pop)]
    return(m.ag)
  }
  priv.dt <- rbindlist(priv.dt)

### save
  save(priv.dt, file="/scratch/aob2x/dest/priv.dt")


## scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/dest/priv.dt ~/.

###
  library(data.table)
  library(ggplot2)
  library(patchwork)

### load

  load("~/priv.dt")
  setnames(priv.dt, "pop", "sampleId")



### pnps
  priv.dt.ag <- priv.dt[,list(pn=sum(n[grepl("missense_variant", ann)]),
                ps=sum(n[grepl("synonymous_variant", ann)]),
                n=sum(n),
                af=median(mean_af),
                ac=median(mean_ac)),
           list(class, sampleId)]
  setkey(priv.dt.ag, sampleId)

### load mps object
  load("/Users/alanbergland/Documents/GitHub/DEST/populationInfo/mps.Rdata")
  setkey(mps, sampleId)
    priv.dt.ag <-  merge(priv.dt.ag, mps)
    priv.dt.ag[class=="common", class:="control_private"]


### pn/ps of private SNPs in SNAPE verus matched controls
  ac.plot <- ggplot() +
  geom_line(data=priv.dt.ag, aes(x=class, y=ac, group=sampleId, color=set))

  af.plot <- ggplot() +
  geom_line(data=priv.dt.ag, aes(x=class, y=af, group=sampleId, color=set))

  priv.control.plot <- ggplot() +
  geom_line(data=priv.dt.ag, aes(x=class, y=pn/ps, group=sampleId, color=set))

  pnps.n.plot <- ggplot(data=priv.dt.ag[class=="private"][auto==T], aes(x=n, y=pn/ps, color=set, shape=class)) + geom_point() +
  ggtitle("Private SNPs only")

  pcr.plot <- ggplot(data=priv.dt.ag[auto==T], aes(x=pcrDup, y=pn/ps, color=set, shape=class)) + geom_point()

  sim.plot <- ggplot(data=priv.dt.ag[auto==T][class=="private"], aes(x=propSimNorm, y=pn/ps, color=set, shape=class)) + geom_point() +
  ggtitle("Private SNPs only")



### merge
(ac.plot | af.plot) / (priv.control.plot | pnps.n.plot) / (pcr.plot | sim.plot) +
plot_annotation(tag_levels = 'A')


   +
  geom_boxplot(data=priv.dt.ag, aes(x=class, y=pn/ps, group=class))
