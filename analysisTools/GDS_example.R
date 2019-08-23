### libraries
  library(SeqArray)
  library(data.table)
  library(doMC)
  registerDoMC(20)
  library(ggplot2)

### open GDS file
  genofile <- seqOpen("/mnt/ssd/DrosEU.DrosRTEC.DPGP.09042019.ref.snpEFF.vcf.reheader.gds", allow.duplicate=T)

### get summary stats per sample for random SNPs

  ### get SNP index data.table
    snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                          pos=seqGetData(genofile, "position"),
                          variant.id=seqGetData(genofile, "variant.id"))

  ### random sample of 1000 SNPs per chromosome
    snp.samp <- snps.dt[,list(variant.id=sample(variant.id, 1000)), list(chr)]

  ### extract out average read depth; seqArray reports the results of this seqGetData command a bit oddly so it takes a little work to get it into the right shape.
    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.samp$variant.id)
    tmp <- seqGetData(genofile, "annotation/format/AD")


    tmp.names <- data.table(col=cumsum(unlist(sapply(tmp$length, function(x) rep(1, x)))),
                            var=unlist(sapply(tmp$length, function(x) list(c("refDP", "altDP")[1:x]))),
                            variant.id=unlist(apply(cbind(snp.samp$variant.id, tmp$length), 1, function(x) rep(x[1], each=x[2]))))

    tmp.ag <- tmp.names[,list(n=length(col)), list(variant.id)]

    tmp.dt <- data.table(refDP = expand.grid(tmp$data[,tmp.names$var=="refDP"])[,1],
                         sample.id = rep(seqGetData(genofile, "sample.id"), dim(snp.samp)[1]),
                         variant.id = rep(snp.samp$variant.id, each=length(seqGetData(genofile, "sample.id"))))



    setkey(tmp.dt, variant.id)
    setkey(snp.samp, variant.id)

     tmp.dt <- merge(tmp.dt, snp.samp)
     tmp.dt.ag <- tmp.dt[,list(medRD=mean(refDP, na.rm=T)), list(chr, sample.id)]

 ### plot
     ggplot(data=tmp.dt.ag[sort(medRD)], aes(x=chr, y=log10(medRD))) + geom_boxplot()







### needs some work
tmp.dt <- data.table(refDP = expand.grid(tmp$data[,tmp.names[tmp.names$variant.id%in%tmp.ag[n==2]$variant.id]$var=="refDP"])[,1],
                     altDP = expand.grid(tmp$data[,tmp.names[tmp.names$variant.id%in%tmp.ag[n==2]$variant.id]$var=="altDP"])[,1],
                     sample.id = rep(seqGetData(genofile, "sample.id"), dim(snp.samp)[1]),
                     variant.id = rep(snp.samp$variant.id, each=length(seqGetData(genofile, "sample.id"))))
