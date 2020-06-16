#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

library(SeqArray)
library(data.table)
library(foreach)
library(ggplot2)

#seqVCF2GDS("/scratch/aob2x/dest/dest.June14_2020.ann.vcf", "/scratch/aob2x/dest.June14_2020.ann.gds", storage.option="ZIP_RA")

args = commandArgs(trailingOnly=TRUE)
vcf.fn=args[[1]]
gds.fn=gsub(".vcf", ".gds", vcf.fn)

seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA")


#genofile <- seqOpen("/scratch/aob2x/dest.June14_2020.maf001.gds")
#samps <- fread("/scratch/aob2x/dest/DEST/populationInfo/samps.csv")
#vSamps <- seqGetData(genofile, "sample.id")
#
## are all populations there?
#setkey(samps, sampleId)
#samps[J(vSamps)]
#samps[!samps$sampleId%in%vSamps]$sampleId
#
## tri-allelic sites
#
### how many are there?
#  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
#                        pos=seqGetData(genofile, "position"),
#                        nAlleles=seqGetData(genofile, "$num_allele"),
#                        id=seqGetData(genofile, "variant.id"))
#  table(snp.dt$chr)
#  prop.table(table(snp.dt$nAlleles))

## Do both alts at tri allelic sites pass filtering?
#  seqSetFilter(genofile, variant.id=snp.dt[nAlleles==3]$id)

#  parseSeqData <- function(x) {
#    #x<- seqGetData(genofile, "annotation/format/AD")

#    x.cs <- cumsum(x$length)
#    x.dt <- data.table(locus=1:length(x.cs), start=x.cs-x$length+1, stop=x.cs)
#    x.dt.sites <- x.dt[stop==(start+1), list(i=start:(stop)), list(locus)]
#    ret <- x$data[,x.dt.sites$i]
#    return(ret)
#  }


#  ad <- parseSeqData(seqGetData(genofile, "annotation/format/AD"))
#  ad1 <- ad[,c(1:dim(ad)[2])[seq(1,dim(ad)[2],by=2)]]
#  ad2 <- ad[,c(2:dim(ad)[2])[seq(1,dim(ad)[2],by=2)]]

#  freq <- parseSeqData(seqGetData(genofile, "annotation/format/FREQ"))
#  freq1 <- freq[,c(1:dim(freq)[2])[seq(1,dim(freq)[2],by=2)]]
#  freq2 <- freq[,c(2:dim(freq)[2])[seq(1,dim(freq)[2],by=2)]]

#  tri <- data.table(locus=c(1:dim(ad1)[2], 1:dim(ad1)[2]),
#                        allele=rep(c(1,2), each=dim(ad1)[2]),
#                       npop=c(apply(ad1, 2, function(x) sum(!is.na(x))),
#                              apply(ad2, 2, function(x) sum(!is.na(x)))),
#                       mac=c(apply(ad1, 2, function(x) sum(x, na.rm=T)),
#                             apply(ad2, 2, function(x) sum(x, na.rm=T))),
#                       af=c(apply(freq1, 2, function(x) mean(x, na.rm=T)),
#                             apply(freq2, 2, function(x) mean(x, na.rm=T))))
#  save(tri, file="~/tri.Rdata")

#  # load("~/tri.Rdata")

#  ggplot(data=tri, aes(x=log10(npop), y=log10(mac), color=as.factor(allele))) + geom_point()
#  ggplot(data=tri, aes(x=af, y=log10(mac), color=as.factor(allele))) + geom_point() + facet_wrap(~allele)
#  ggplot(data=tri, aes(x=af, y=log10(npop), color=as.factor(allele))) + geom_point() + facet_wrap(~allele)


#  ggplot(data=tri, aes(x=af, y=log10(npop), color=as.factor(allele))) + geom_point() + facet_wrap(~allele)



#  mean(tri[,list(bad=any(mac<=4)), locus]$bad)
