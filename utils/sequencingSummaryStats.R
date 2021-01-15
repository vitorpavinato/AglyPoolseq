### module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(Rsamtools)

### functions
  getRD <- function(gds.fn, q) {
    ### open GDS file
      genofile <- seqOpen(gds.fn, allow.duplicate=T)

    ### get SNP index data.table
      snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                            pos=seqGetData(genofile, "position"),
                            variant.id=seqGetData(genofile, "variant.id"),
                            nAlleles=seqNumAllele(genofile))

    ### get read depths
      nSNPs <- 10000
      seqSetFilter(genofile, variant.id=as.numeric(sample(as.character(snps.dt[nAlleles==2]$variant.id), nSNPs)))
      #tmp.ad <- seqGetData(genofile, "annotation/format/AD")
      #tmp.rd <- seqGetData(genofile, "annotation/format/RD")
      tmp.dp <- seqGetData(genofile, "annotation/format/DP")

      dat <- data.table(dp=expand.grid(tmp.dp$data)$Var1,
                        sampleId=rep(seqGetData(genofile, "sample.id"), dim(tmp.dp$data)[2]),
                        variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(tmp.dp$data)[1]))
      dat[,variant.id:=as.numeric(variant.id)]

      setkey(snps.dt, variant.id)
      setkey(dat, variant.id)

      dat <- merge(dat, snps.dt)
      setkey(dat, chr)
      dat.ag <- dat[J(c("2L", "2R", "3L", "3R", "X")),
                      list(mu=mean(dp, na.rm=T), nmissing=sum(is.na(dp)), q=q), list(sampleId, auto=!(chr=="X"))]

      setnames(dat.ag, c("mu", "nmissing"), paste(c("mu", "nmissing"), q, sep="."))
      setkey(dat.ag, sampleId, auto)

      seqClose(genofile)
      return(dat.ag)
  }

### open GDS file
  dat <- getRD(gds.fn="/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.10Nov2020.ann.gds", q=25)

### incorporate metadata
  samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")

  m <- merge(dat, samps, by="sampleId", all.y=T)


### pull in PCR dup rate
  fns <- system("ls /project/berglandlab/DEST/dest_mapped/*/*/*duplicates_report.txt", intern=T)

  pcrs <- foreach(fn=fns)%do%{
    #fn <- fns[1]
    data.table(pcrDup=fread(fn, skip="LIBRARY", nrows=1)$PERCENT_DUPLICATION, sampleId=tstrsplit(fn, "/")[[7]])
  }
  pcrs <- rbindlist(pcrs)
  setkey(pcrs, sampleId)

  mp <- merge(m, pcrs, all.x=T)

### simulans contamination rate
  simContam <- foreach(samp.i=m[set!="dgn"][auto==T]$sampleId)%do%{
    #samp.i=m[set!="dgn"][auto==T]$sampleId[7]
    message(samp.i)
    simBam <- gsub("mark_duplicates_report.txt", "sim.bam", fns[grepl(samp.i, fns)])
    simIdx <- paste(simBam, "bai", sep=".")

    melBam <- gsub("mark_duplicates_report.txt", "mel.bam", fns[grepl(samp.i, fns)])
    melIdx <- paste(melBam, "bai", sep=".")

    if(!file.exists(melIdx)) {
      message("indexing Mel")
      indexBam(melBam)
    }

    if(!file.exists(simIdx)) {
      message("indexing Sim")

      indexBam(simBam)
    }

    simidx.out <- as.data.table(idxstatsBam(file=simBam, index=simIdx))[grepl("2L|2R|3L|3R|X|^4$|Y", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]
    melidx.out <- as.data.table(idxstatsBam(file=melBam, index=melIdx))[grepl("2L|2R|3L|3R|X|^4$|Y", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]

    idx.out <- merge(melidx.out, simidx.out, by="seqnames")

    idx.out[,nReads:=mapped.x + mapped.y]
    idx.out[,simChr:=grepl("sim", seqnames)]
    idx.out[,chr:=gsub("sim_", "", seqnames)]
    idx.out[,nReadsNorm:=nReads/seqlength.x]


    idx.out.ag <- idx.out[,list(propSim=nReads[simChr==T]/sum(nReads),
                                propSimNorm=nReadsNorm[simChr==T]/sum(nReadsNorm),
                                nMelReads=nReads[simChr==F],
                                melChrLen=seqlength.x[simChr==F],
                                sampleId=samp.i),
             list(chr)]
   idx.out.ag[,mappingEffort:=nMelReads/melChrLen]


    idx.out.ag
  }
  simContam <- rbindlist(simContam)
  simContam.ag <- simContam[,list(propSimNorm=mean(propSimNorm, na.rm=T)), sampleId]


  setkey(mp, sampleId)
  setkey(simContam.ag, sampleId)

### diversity statistics (PoolSNP)
  fn <- list.files("/scratch/aob2x/dest/PoolGenOut_PoolSNP/", full.names=T)
  div <- foreach(fn.i=fn)%dopar%{
    message(paste(which(fn.i==fn), length(fn), sep=" / "))
    #fn.i <- fn[1]
    div <- fread(fn.i)
    setnames(div, names(div), c("chr", "mid", "est", "V4"))
    div[,sampleId:=tstrsplit(last(tstrsplit(fn.i, "/")), "\\.")[[1]]]
    div[,stat:=last(tstrsplit(fn.i, "\\."))]
    div
  }
  div <- rbindlist(div)
  div.ag <- div[,list(div_estimate=median(est, na.rm=T)), list(chr, stat, sampleId)]

  mpd <- merge(mp, div.ag, all=T, by="sampleId", allow.cartesian=T)
  save(mpd, div, file="~/mpd.Rdata")

  scp aob2x@rivanna.hpc.virginia.edu:~/mpd.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/mpd.Rdata")

  ggplot(data=mpd[!is.na(stat)][!is.na(Continental_clusters)], aes(x=chr, y=div_estimate, color=Continental_clusters)) + geom_boxplot() + facet_wrap(~stat, scales="free")


  ggplot(data=mpd[!is.na(stat)][!is.na(Continental_clusters)], aes(x=lat, y=div_estimate, color=Continental_clusters)) + geom_point() + facet_wrap(~stat, scales="free")


  summary(lm(div_estimate~lat*Continental_clusters, mpd[!is.na(stat)][!is.na(Continental_clusters)]))
  mpd[,dist_ch:=sqrt((lat-51.2763)^2 +  (long-30.2219)^2)]

  ggplot(data=mpd[!is.na(stat)][!is.na(Continental_clusters)], aes(x=dist_ch, y=div_estimate, color=Continental_clusters)) + geom_point() + facet_wrap(~stat, scales="free")
  

#### new vs. odl
  library(gdata)


  div.new.old <- foreach(sheet.i=c("pi", "Theta", "TajimaD", "pi_40x", "Theta_40x", "TajimaD_40x"))%dopar%{
    message(sheet.i)
    div.old <- read.xls("/scratch/aob2x/dest/DEST/utils/DrosEU2014_PopGen.xlsx", sheet=sheet.i)
    div.old.long <- as.data.table(melt(div.old, id.vars=c("Chrom", "Pos")))
    setnames(div.old.long, c("Chrom", "Pos", "variable"), c("chr", "mid", "sampleId"))

    if(grepl("pi", sheet.i)) {
      div.old.long[,stat:="pi"]
      div.old.long[,stat_old:=sheet.i]

    } else if(grepl("Theta", sheet.i)) {
      div.old.long[,stat:="th"]
      div.old.long[,stat_old:=sheet.i]

    } else if(grepl("Tajima", sheet.i)) {
      div.old.long[,stat:="D"]
      div.old.long[,stat_old:=sheet.i]

    }

    setkey(div.old.long, chr, mid, stat, sampleId)
    setkey(div, chr, mid, stat, sampleId)
    tmp <- merge(div, div.old.long)
    tmp
  }
  div.new.old <- rbindlist(div.new.old)
  save(div.new.old, file="~/div_new_old.Rdata")


scp aob2x@rivanna.hpc.virginia.edu:~/div_new_old.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/div_new_old.Rdata")
  div.ag <- div.new.old[,list(new=median(est, na.rm=T), old=mean(value, na.rm=T)), list(chr, stat, sampleId, downsample=paste("40x:", grepl("40x", stat_old), sep=" "))]
  div.ag2 <- div.new.old[,list(new=mean(est, na.rm=T), old=mean(value, na.rm=T)), list(chr, mid, stat, sampleId, downsample=paste("40x:", grepl("40x", stat_old), sep=" "))]
  save(div.ag2, file="div_new_old.Rdata")

  ggplot(div.ag, aes(x=old, y=new, color=chr)) + geom_point() + facet_wrap(downsample~stat, scales="free")
  ggplot(div.ag2, aes(x=old, y=new, color=chr)) + geom_point() + facet_wrap(downsample~stat, scales="free") + geom_abline(intercept=0,slope=1)









  mps <- merge(mp, simContam.ag, all.x=T)

  save(mps, simContam, file="~/mps.Rdata")
  save(mps, simContam, file="/scratch/aob2x/dest/DEST/populationInfo/mps.Rdata")

###
  # scp aob2x@rivanna.hpc.virginia.edu:~/mps.Rdata ~/.

### plots
  library(ggplot2)
  library(data.table)
  library(cowplot)

### load adata
  load("~/mps.Rdata")
  setnames(mps, "mu.15", "AveReadDepth.15")
  setnames(mps, "mu.25", "AveReadDepth.25")

### a few small fixes
  mps[continent=="North_America", continent:="NorthAmerica"]

  mps[,effRD.15:=(AveReadDepth.15 * 2*nFlies) / (AveReadDepth.15 + 2*nFlies)]
  mps[,effRD.25:=(AveReadDepth.25 * 2*nFlies) / (AveReadDepth.25 + 2*nFlies)]

### rank x-axis to mean read depth
  mps <- mps[auto==T]
  mps[,x:=rank(AveReadDepth, ties.method="first")]
  mps[,x.id:=factor(sampleId, levels=mps$sampleId[mps$AveReadDepth])]

  mps[,propMissing:=nmissing/10000]

### wide to long
  mpsl <- melt(mps,
              id.vars=c("x", "sampleId", "continent", "auto", "set"),
              measure.vars=c("AveReadDepth", "propMissing", "nFlies", "effRD", "pcrDup", "propSimNorm"))

  mpsl[,xf:=as.factor(x)]

### summary plot
  summaryStat.plot <- ggplot(data=mpsl, aes(x=xf, y=value, color=continent, fill=continent)) +
          geom_point(pch=21, alpha=.5, size=2) +
          facet_grid(variable~set, scales="free", space="free_x") +
          geom_point(data=mpsl[sampleId%in%mps[propMissing>.1]$sampleId], aes(x=xf, y=value), color="black", size=.5) +
          theme(axis.text.x = element_blank(), panel.spacing = unit(1, "lines"))


### mapping rates per chr
  mappingRat.ag <- simContam[,list(sumMelReads=sum(nMelReads), sumMelGenome=sum(melChrLen)), list(sampleId)]
  setkey(mappingRat.ag, sampleId)
  mappingRate <- merge(simContam, mappingRat.ag, by="sampleId")
  mappingRate[,mr:=(nMelReads/sumMelReads)/(melChrLen/sumMelGenome)]

  mrmps <- merge(mps, mappingRate, by="sampleId")
  mrmps[,chr:=factor(chr, levels=c("2L", "2R", "3L", "3R", "X", "Y", "4"))]

  mpsl.ag <- mpsl[,list(xf=mean(x)), list(sampleId)]

  mrmps <- merge(mrmps, mpsl.ag, by="sampleId")
  mrmps[,xf:=as.factor(x)]
  #mrmps[,continent:=factor(continent, levels=unique(mps$continent))]


  mappingRate.plot <- ggplot(data=mrmps, aes(x=xf, y=mr, color=continent, fill=continent)) +
  geom_point(pch=21, alpha=.5, size=2) +
  facet_grid(chr~set, scales="free_x", space="free_x") +
  geom_point(data=mrmps[propMissing>.1], aes(x=xf, y=mr), color="black", size=.5) +
  theme(axis.text.x = element_blank(), panel.spacing = unit(1, "lines"))

  summaryPlot <- summaryStat.plot + mappingRate.plot

  ggsave(summaryPlot, file="~/summaryPlot.Aug9_2020.001.50.pdf")




  ggsave(summaryStat.plot, file="~/summaryStat_plot.pdf", h=10, w=8)


### any collinearity with propMissing? Not really
  cor(mps[, c("AveReadDepth", "propMissing", "nFlies", "effRD", "pcrDup", "propSim"), with=F], use="complete")

  summary(lm(propMissing ~ propSim * pcrDup * , mps))
  summary(lm(propSim ~ pcrDup  , mps))

### labeled samples
    missing <- ggplot(data=mps[auto==T], aes(x=as.factor(mu.r), y=nmissing/10000, color=continent, fill=continent)) +
    geom_point(pch=21, alpha=.5, size=2) +
    facet_grid(auto~set, scales="free_x", space="free_x") +
    theme(axis.text.x = element_blank()) +
    geom_text(aes(label=ifelse((nmissing/10000)>.1, as.character(sampleId),'')), hjust=0,vjust=0)
