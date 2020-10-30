


#### this is all done in awk
zcat dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | grep -v "#" | awk '{
npoly=0
nmissing=0
pops=""
if(length($5)==1) {
  for(i=10; i<=NF; i++) {
    split($i, sp, ":")
    if(sp[1]=="0/1") npoly++
    if(sp[1]=="0/1") pops=pops";"i-9
    if(sp[1]=="./.") nmissing++
    }
    if(npoly<25) print $1"\t"$2"\t"npoly"\t"nmissing"\t"$4"\t"$5"\t"pops
  }
}' > /scratch/aob2x/dest/polym.delim





#### now R
#### libraries
  library(data.table)


### load sample names
  pops <- names(fread("zcat dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | head -n40", nrows=1, skip="#CHR"))[-(1:9)]

### samps
  samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")
  setkey(samps, sampleId)

### private
  priv.dt <- fread("/scratch/aob2x/dest/polym.delim")
  priv.dt <- priv.dt[V3==2]

  dim(priv.dt)[1]
  o <- foreach(i=1:dim(priv.dt)[1])%dopar%{

    if(i%%500==0) message(paste(i, dim(priv.dt)[1], sep=" / "))
    set <- samps[J(pops[as.numeric(unlist(tstrsplit(priv.dt$V7[i], ";"))[-1])])]
    pw.dist <- spDists(x=as.matrix(set[,c("long", "lat"), with=F]), longlat=T)

    cbind(priv.dt[i], data.table(meanDist=mean(pw.dist[lower.tri(pw.dist)])))

  }
  o <- rbindlist(o)

  pw.dist <- spDists(x=as.matrix(samps[set%in%c("DrosRTEC", "DrosEU"), c("long", "lat"), with=F]), longlat=T)
  diag(pw.dist) <- NA
  oe <- foreach(i=1:dim(priv.dt)[1])%dopar%{

    if(i%%500==0) message(paste(i, dim(priv.dt)[1], sep=" / "))

    pw.dist[sample(246,1), sample(246,1)]

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
