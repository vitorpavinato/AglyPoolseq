#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

library(SeqArray)
library(data.table)

seqVCF2GDS("~/test.all.vcf", "~/tmp.gds", storage.option="ZIP_RA")
genofile <- seqOpen("~/tmp.gds")
samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")
vSamps <- seqGetData(genofile, "sample.id")


setkey(samps, sampleId)

samps[J(vSamps)]
