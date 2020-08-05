module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(Rsamtools)

### per mapping rate

  fns <- system("ls /project/berglandlab/DEST/dest_mapped/*/*/*mel.bam", intern=T)



mappingRate <- foreach(fn.i=fns)%do%{
  #fn.i=fns[1]
  print(fn.i)
  melidx.out <- as.data.table(idxstatsBam(file=fn.i, index=paste(fn.i, ".bai", sep="")))[grepl("2L|2R|3L|3R|X", seqnames)][!grepl("Het|het|Sac|Sca", seqnames)]

  melidx.out[,samp:=last(tstrsplit(fn.i, "/"))]

  melidx.out
}

mappingRate <- rbindlist(mappingRate)
mappingRate[,rate:=mapped/seqlength]

save(mappingRate, file="~/mappingRate.Rdata")


# scp aob2x@rivanna.hpc.virginia.edu:~/mappingRate.Rdata ~/.

library(ggplot2)
load("~/mappingRate.Rdata")

mappingRate.ag <- mappingRate[,list(sum=sum(mapped)), list(samp)]
mappingRate <- merge(mappingRate, mappingRate.ag, by="samp")
mappingRate[,rate:=mapped/sum]
ggplot(data=mappingRate[!grepl("sim", seqnames)], aes(x=samp, y=rate)) +geom_point() + facet_grid(~seqnames)
