#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### capture input argument: 1-246

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### open GDS for common SNPs (PoolSNP)
  genofile.PoolSeq <- seqOpen("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.ann.gds", allow.duplicate=T)
  length(seqGetData(genofile.PoolSeq, "variant.id"))
  seqClose(genofile.PoolSeq)

  genofile.all <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.ann.gds", allow.duplicate=T)
  length(seqGetData(genofile.all, "variant.id"))
  seqClose(genofile.all)

### get summaries of private SNps
  fn <- list.files("/scratch/aob2x/dest/privateSNPs", ".summary.delim", full.names=T)

  o <- foreach(fn.i=fn)%do%{
    #fn.i=fn[1]
    fread(fn.i)
  }
  o <- rbindlist(o, fill=T)
  o[is.na(set), set:="SNAPE"]

  o.ag <- o[,list(nPrivate=sum(nPrivate)), list(set)]
