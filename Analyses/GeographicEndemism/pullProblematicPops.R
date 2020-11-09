#
### libraries
	library(data.table)
  library(googlesheets4)

### get data from shared google sheet
sheet.id <- "1PJEQE6Pn8H_vwA0cgn8e4nECrWrSGzKj6kgg7b8a3Ms"

pops <- as.data.table(read_sheet(ss=sheet.id))
setnames(pops, "problematic runs", "prob")

write.table(pops[is.na(prob)]$Sample, row.names=F, col.names=F, quote=F, file="/Users/alanbergland/Documents/GitHub/DEST/Analyses/GeographicEndemism/goodSamps.delim")
