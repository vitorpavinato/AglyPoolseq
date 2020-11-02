
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

### make expectation MA table and do a bit of rearragment of MT data
  exp.mt <- data.table(mt=c("AT", "AC", "AG", "TA", "TC", "TG", "CA", "CT", "CG", "GA", "GT", "GC"),
              rate=c(6, 3, 5.5, 6, 5.5, 3., 9, 20, 6.5, 20, 9, 6.5)/100)


### plots

    mutation.plot <- ggplot(data=o2.ag, aes(x=nPop, y=freq, group=mt, color=mt)) + geom_line() + cowplot::theme_cowplot() +
    xlab("Number of polymorphic populations") + ylab("Frequency") +
    geom_hline(data=exp.mt, aes(yintercept=rate, color=mt), linetype="dashed") +
    geom_hline(data=dgrp.ag, aes(yintercept=freq, color=mt), linetype="solid") +
    facet_grid(caller~mt) +
    scale_linetype_manual(values = c("Assaf" = 1, "DGRP" = 2))


    dist.plot <- ggplot(data=o.ag, aes(x=nPop, y=mean.dist, group=set, color=set)) +
    geom_ribbon(aes(x=nPop, ymin=mean.dist-sd.dist, ymax=mean.dist+sd.dist, group=set, fill=set), alpha=.5) +
    geom_line(size=1) +
    facet_grid(~caller) +
    xlab("Number of polymorphic populations") +
    ylab("Average geographic distance \nbetween populations (km") +
    cowplot::theme_cowplot()



    phyloConcord.plot <- ggplot(data=o.ag, aes(x=nPop, y=log10(mean.cc_equal), group=set, color=set)) +
    geom_line() +
    facet_grid(~caller) +
    xlab("Number of polymorphic populations") +
    ylab("Probability that all populations \nare in same phylogeographic cluster") +
    cowplot::theme_cowplot()

    n.plot <- ggplot(o4.ag, aes(x=(nPop), y=log10(nSites), fill=caller)) +
    geom_col(position = "dodge", size=2) +
    cowplot::theme_cowplot()

    ggplot(na.omit(o3.ag), aes(x=(nPop), y=maf.bin, fill=delta.dist)) + geom_tile(interpolate=T) +
    scale_fill_viridis()

design <- "
AAABBBCCC
AAABBBCCC
AAABBBCCC
"

design <- "
AAAA
AAAA
BBCC
BBCC
"

design <- "
AAAAAA
AAAAAA
BBBCCC
BBBCCC
"



o.plot <- (n.plot + dist.plot + phyloConcord.plot)  + plot_layout(design = design)
o.plot

ggsave(o.plot, file="~/geo_endemic.pdf")







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
