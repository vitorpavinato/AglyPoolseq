### module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(SeqArray)
  library(data.table)

### open GDS file
  gds.fn <- "/scratch/aob2x/dest/dest.July6_2020.001.10.ann.gds"
  genofile <- seqOpen(gds.fn, allow.duplicate=T)


### get SNP index data.table
  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile))

### get depths
  seqSetFilter(genofile, variant.id=snps.dt[nAlleles==2]$variant.id[1:1000])
  #tmp.ad <- seqGetData(genofile, "annotation/format/AD")
  #tmp.rd <- seqGetData(genofile, "annotation/format/RD")
  tmp.dp <- seqGetData(genofile, "annotation/format/DP")

  dat <- data.table(dp=expand.grid(tmp.dp$data)$Var1,
                    population=rep(seqGetData(genofile, "sample.id"), dim(tmp.dp$data)[2]),
                    variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(tmp.dp$data)[1]))

  dat.ag <- dat[,list(mu=mean(dp, na.rm=T), nmissing=sum(is.na(dp))), list(population)]
