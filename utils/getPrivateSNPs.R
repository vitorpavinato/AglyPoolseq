
# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### capture input argument: 1-246

args = commandArgs(trailingOnly=TRUE)
pop.n <- as.numeric(args[1])
set <- args[2]

# pop.n <- 1; set<-"PoolSNP"

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### function
  getPrivateSNPs <- function(pop.i, write.out=T) {
    message(pop.i)
    #pop.i<-"FL_ho_10_fall"
    #pop.i <- pops[2]
    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt$id, sample.id=pop.i)
    tmp <- seqGetData(genofile, "annotation/format/AD")[[2]]

    pop.dt <- snp.dt
    pop.dt[,pop.ac:=t(tmp)[,1]]
    pop.dt[,pop:=pop.i]
    pop.dt[,set:=set]


    if(write.out==T) {
      out.fn <- paste("/scratch/aob2x/dest/privateSNPs/", pop.i, ".", set, ".csv", sep="")

      write.csv(pop.dt[ac==pop.ac], row.names=F, quote=F, file=out.fn)

      pop.ag <- pop.dt[,list(nPrivate=sum(ac==pop.ac, na.rm=T),
                          nPoly=sum(pop.ac>0, na.rm=T),
                          nInform=sum(!is.na(pop.ac)),
                          n=length(pop.ac)),
                      list(chr, pop, set)]

      out.fn <- gsub(".csv", ".summary.delim", out.fn)
      write.table(pop.ag, out.fn, sep="\t", quote=F, row.names=F)
    }

    return(pop.dt[,list(nPrivate=sum(ac==pop.ac, na.rm=T),
                        nPoly=sum(pop.ac>0, na.rm=T),
                        nInform=sum(!is.na(pop.ac)),
                        n=length(pop.ac)),
                    list(chr, pop, set)])

  }

### open GDS file
  if(set=="PoolSNP") genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.ann.gds", allow.duplicate=T)
  if(set=="SNAPE") genofile <- seqOpen("/scratch/aob2x/dest/dest.PoolSeq.SNAPE.001.50.ann.gds", allow.duplicate=T)
  if(set=="all") genofile <- seqOpen("/scratch/aob2x/dest/dest.all.PoolSNP.001.50.ann.gds", allow.duplicate=T)

### make SNP table
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        id=seqGetData(genofile, "variant.id"))
  snp.dt <- snp.dt[nAlleles==2]
  seqSetFilter(genofile, snp.dt$id)

  snp.dt[,adp:=seqGetData(genofile, "annotation/info/ADP")]
  snp.dt[,ac:=seqGetData(genofile, "annotation/info/AC")[[2]]]

### iterate through populations
  seqResetFilter(genofile)
  pops <- seqGetData(genofile, "sample.id")
  private.dt <- foreach(i=pops[pop.n])%do%getPrivateSNPs(pop.i=i)
  private.dt

  seqClose(genofile)
