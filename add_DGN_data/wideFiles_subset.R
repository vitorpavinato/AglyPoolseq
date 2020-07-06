### libraries
  library(data.table)

### set wd
  setwd("/scratch/aob2x/dest")

### load DGN samps to use
  dgn.use <- fread("DEST/populationInfo/dpgp.ind.use.csv")

### load DGN available
  dgn.available <- fread("dgn/dgn_wideFiles.delim", sep="\t", sep2="_")
  dgn.available[,sampleId:=tstrsplit(V2, "_")[[1]]]
  dgn.available[,Stock.ID:=tstrsplit(V2, "_")[[2]]]

### merge
  setkey(dgn.use, "sampleId", "Stock.ID")
  setkey(dgn.available, "sampleId", "Stock.ID")

  dgn.final <- merge(dgn.use, dgn.available)


  setkey(dgn.final, V2)

  dgn.final <- dgn.final[!duplicated(V2)]
  dgn.final <- dgn.final[,V1:=c(1:dim(dgn.final)[1])][,c("V1", "V2"), with=F]

  write.table(dgn.final, file="dgn/dgn_wideFiles.final.delim", quote=F, row.names=F, sep="\t", col.names=F)
