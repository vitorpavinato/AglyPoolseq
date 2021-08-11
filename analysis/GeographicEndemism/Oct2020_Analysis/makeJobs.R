
library(data.table)
tmp <- data.table(maf=rep(seq(from=1, to=7, by=2), each=3)/c(1000, 100, 10), chr=rep(c("2L", "2R", "3L", "3R", "X"), each=12))
tmp <- tmp[maf<=.3]
write.table(tmp, file="/scratch/aob2x/dest/geo_endemic/jobs.txt", col.names=F, row.names=F, sep="\t")
