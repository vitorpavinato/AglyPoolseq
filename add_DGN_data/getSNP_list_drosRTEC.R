### libraries
  library(data.table)
  library(foreach)

### functions
  loadData <- function(obj) {
   ### status

           print(obj)

   ### load data
           load(obj)

   ### twist around frequency object
           freq$chr <- info$X.CHROM
           freq$pos <- info$POS

           freql <- melt(freq, id.vars=c("chr", "pos"))
           freql <- as.data.table(freql)
           setnames(freql, c("variable", "value"), c("pop", "af"))
           setkey(freql, chr, pos, pop)

   ### twist around frequency object
           dp$chr <- info$X.CHROM
           dp$pos <- info$POS
           dp <- as.data.table(dp)
           setnames(dp, names(dp), gsub("\\.1", "", names(dp)))

           dpl <- melt(dp, id.vars=c("chr", "pos"))
           setnames(dpl, c("variable", "value"), c("pop", "dp"))
           setkey(dpl, chr, pos, pop)

   ### merge
           merge(freql, dpl)
         }

### data
  fileList <- system("ls /scratch/aob2x/dest/drosRTEC/mel_freqdp*.Rdata", intern=T)
  dat <- foreach(i=fileList)%do%loadData(i)
  dat <- rbindlist(dat)

  setkey(dat, chr, pos)

### extract out sites
  dat <- dat[!duplicated(dat)]

### tack in chr
  dat[,chr:=gsub("chr", "", chr)]
  dat[,chr:=paste("chr", chr, sep="")]
  dat[,start:=pos]
  dat[,stop:=pos+1]

### disable scientific notation
  options(scipen = 999)

### export
  write.table(dat[,c("chr", "start", "stop"), with=F], file="/scratch/aob2x/dest/drosRTEC/drosRTEC_sites.dm3.bed",
                          quote=F, row.names=F, col.names=F, sep="\t")
