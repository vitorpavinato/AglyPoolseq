#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

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



o.ag <- o[,list(mu=mean(meanDist), sd=sd(meanDist),
                med=median(meanDist), lci=quantile(meanDist, .025), uc=quantile(meanDist, .975)),
           list(set, nPop=V3)]

o2.ag <- o[,list(.N), list(mt, nPop=V3)]
o2.ag.freq <- na.omit(o2.ag)[,list(freq=N/sum(N), mt), list(nPop)]

save(o.ag, o2.ag, file="~/allSummarySet_endemism.Rdata")


  scp aob2x@rivanna.hpc.virginia.edu:~/allSummarySet_endemism.Rdata ~/.
  library(ggplot2)
  library(data.table)

  load("~/allSummarySet_endemism.Rdata")




    ggplot(data=o.ag, aes(x=nPop, y=med, group=set, color=set)) +
    geom_line()


    ggplot(data=o.ag, aes(x=nPop, y=med, group=set, color=set)) +
    geom_ribbon(aes(x=nPop, ymin=lci, ymax=uc, group=set, fill=set), alpha=.5) +
    geom_line()


    mutationplot <- ggplot(data=o2.ag.freq, aes(x=nPop, y=freq, group=mt, color=mt)) + geom_line() + cowplot::theme_cowplot() +
    xlab("Number of polymorphic populations") + ylab("Frequency")





    distplot <- ggplot(data=o.ag, aes(x=nPop, y=mu, group=set, color=set)) +
    geom_ribbon(aes(x=nPop, ymin=mu-sd, ymax=mu+sd, group=set, fill=set), alpha=.5) +
    geom_line(size=1) +
    xlab("Number of polymorphic populations") +
    ylab("Average distance between populations") + cowplot::theme_cowplot()


distplot | mutationplot




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
