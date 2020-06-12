#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

library(SeqArray)
library(data.table)

seqVCF2GDS("~/test.all.vcf", "~/tmp.gds", storage.option="ZIP_RA")
genofile <- seqOpen("~/tmp.gds")
samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")
vSamps <- seqGetData(genofile, "sample.id")


setkey(samps, sampleId)

samps[J(vSamps)]

samps[!samps$sampleId%in%vSamps]$sampleId



snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"), pos=seqGetData(genofile, "position"), nAlleles=seqGetData(genofile, "$num_allele"), id=seqGetData(genofile, "variant.id"))

seqSetFilter(genofile, variant.id=1200)

seqGetData(genofile, "annotation/format/FREQ")
seqGetData(genofile, "annotation/format/AD")
seqGetData(genofile, "annotation/format/RD")
seqGetData(genofile, "annotation/format/DP")

seqGetData(genofile, "$dosage")
seqGetData(genofile, "@genotype")
seqGetData(genofile, "allele")
