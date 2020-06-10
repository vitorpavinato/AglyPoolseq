#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

args = commandArgs(trailingOnly=TRUE)
job=args[1]
tmpdir=args[2]

#job="2L_41301_55067"; tmpdir="/scratch/aob2x/test"
jobId=gsub(",", "_", job)

### libraries
  library(data.table)
  library(foreach)

### get input files
  files <- list.files(tmpdir, pattern=jobId)

### import
  o <- foreach(files.i=files)%do%{
    #files.i=files[199]
    tmp <- fread(files.i)
    tmp[,pop:=gsub("_$", "", gsub(jobId, "", files.i))]
    tmp
  }
  o <- rbindlist(o)

### long to wide
  ow <- dcast(o, V1+V2~pop, value.var="V4")

## get reference
  ow.ref <- o[pop=="AT_gr_12_fall", c("V1", "V2", "V3"), with=F]

  setkey(ow, V1, V2)
  setkey(ow.ref, V1, V2)

  owr <- merge(ow.ref, ow)

### output
  write.table(owr, quote=F, row.names=F, col.names=F, sep="\t", file=paste(tmpdir, "/allpops.sites", sep=""))
  write.table(names(owr)[-c(1,2,3)], quote=F, row.names=F, col.names=F, sep="\t", file=paste(tmpdir, "/allpops.names", sep=""))
