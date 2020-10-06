#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### function
  getPrivateSNPs <- function(pop.i) {
    message(pop.i)
    #pop.i<-"FL_ho_10_fall"
    #pop.i <- pops[2]
    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt$id, sample.id=pop.i)
    tmp <- seqGetData(genofile, "annotation/format/AD")[[2]]

    pop.dt <- snp.dt
    pop.dt[,pop.ac:=t(tmp)[,1]]
    pop.dt[,pop:=pop.i]

    return(pop.dt[,list(nPrivate=sum(ac==pop.ac, na.rm=T),
                        nPoly=sum(pop.ac>0, na.rm=T),
                        nInform=sum(!is.na(pop.ac)),
                        n=length(pop.ac)),
                    list(chr, pop)])

  }

### open GDS file
  genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.ann.gds", allow.duplicate=T)

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
  private.dt <- foreach(i=pops[1:50])%dopar%getPrivateSNPs(pop.i=i)
  private.dt <- rbindlist(private.dt)
