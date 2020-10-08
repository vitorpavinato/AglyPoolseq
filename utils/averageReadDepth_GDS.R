### module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(Rsamtools)

### open GDS file
  gds.fn <- "/scratch/aob2x/dest/dest.July6_2020.001.10.ann.gds"
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
                  list(mu=mean(dp, na.rm=T), nmissing=sum(is.na(dp))), list(sampleId, auto=!(chr=="X"))]

### incorporate metadata
  samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")

  m <- merge(dat.ag, samps)


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
    #samp.i=m[set!="dgn"][auto==T]$sampleId[1]
    message(samp.i)
    simBam <- gsub("mark_duplicates_report.txt", "sim.bam", fns[grepl(samp.i, fns)])
    simIdx <- paste(simBam, "bai", sep=".")

    melBam <- gsub("mark_duplicates_report.txt", "mel.bam", fns[grepl(samp.i, fns)])
    melIdx <- paste(melBam, "bai", sep=".")

    simidx.out <- as.data.table(idxstatsBam(file=simBam, index=simIdx))[grepl("2L|2R|3L|3R|X", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]
    melidx.out <- as.data.table(idxstatsBam(file=melBam, index=melIdx))[grepl("2L|2R|3L|3R|X", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]

    idx.out <- merge(melidx.out, simidx.out, by="seqnames")

    data.table(propSim=sum(idx.out$mapped.y)/(sum(idx.out$mapped.y) + sum(idx.out$mapped.x)),
              sampleId=samp.i)

  }
  simContam <- rbindlist(simContam)
  setkey(mp, sampleId)
  setkey(simContam, sampleId)
  mps <- merge(mp, simContam, all.x=T)

  save(mps, file="~/mps.Rdata")
  save(mps, file="/scratch/aob2x/dest/DEST/populationInfo/mps.Rdata")

###
  # scp aob2x@rivanna.hpc.virginia.edu:~/mps.Rdata ~/.

### plots
  library(ggplot2)
  library(data.table)
  library(cowplot)

### load adata
  load("~/mps.Rdata")

### a few small fixes
  mps[continent=="North_America", continent:="NorthAmerica"]

  mps[,effRD:=(mu* 2*nFlies) / (mu + 2*nFlies)]
  setnames(mps, "mu.15", "AveReadDepth.15")
  setnames(mps, "mu.25", "AveReadDepth.25")
  mps[,effRD.15:=(AveReadDepth.15 * 2*nFlies) / (AveReadDepth.15 + 2*nFlies)]
  mps[,effRD.25:=(AveReadDepth.25 * 2*nFlies) / (AveReadDepth.25 + 2*nFlies)]

### rank x-axis to mean read depth
  mps <- mps[auto==T]
  mps[,x:=rank(AveReadDepth.25, ties.method="first")]
  mps[,x.id:=factor(sampleId, levels=mps$sampleId[mps$mu.r])]

  mps[,propMissing.25:=nmissing.25/10000]

### wide to long
  mpsl <- melt(mps,
              id.vars=c("x", "sampleId", "continent", "auto", "set"),
              measure.vars=c("AveReadDepth.25", "propMissing.25", "nFlies", "effRD.25", "pcrDup", "propSim"))

  mpsl[,xf:=as.factor(x)]

### plot
  summaryStat.plot <- ggplot(data=mpsl, aes(x=xf, y=value, color=continent, fill=continent)) +
          geom_point(pch=21, alpha=.5, size=2) +
          facet_grid(variable~set, scales="free", space="free_x") +
          theme(axis.text.x = element_blank(), panel.spacing = unit(1, "lines")) +
          geom_hline(data=data.table(variable="propMissing", value=.1), aes(yintercept=value))

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

###
  library(reticulate)
  # py_install("pandas")
  source_python("/scratch/aob2x/dest/DEST/misc_obsolete/pickle_reader.py")
  cov <- read_pickle_file("C:/tsa/dataset.pickle")





drosEU.pcr <- ggplot(data=mps[set=="DrosEU"], aes(x=year, y=pcrDup)) +
geom_point() +
geom_text(data=mps[set=="DrosEU"][pcrDup>.25], aes(x=year, y=pcrDup, label=sampleId), check_overlap=TRUE, hjust = 0, nudge_x = 0.05) +
ylim(0,1)



drosRTEC.pcr <- ggplot(data=mps[set=="DrosRTEC"], aes(x=year, y=pcrDup)) +
geom_point() +
geom_text(data=mps[set=="DrosRTEC"][pcrDup>.25], aes(x=year, y=pcrDup, label=sampleId), check_overlap=TRUE, hjust = 1, nudge_x = 0.1) +
ylim(0,1)



#library(patchwork)
drosEU.pcr | drosRTEC.pcr
