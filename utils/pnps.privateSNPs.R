# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### capture input argument: 1-246

args = commandArgs(trailingOnly=TRUE)
pop.n <- as.numeric(args[1])
set <- args[2]

# pop.n <- 1

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### open GDS for common SNPs (PoolSNP)
  genofile.common <- seqOpen("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.ann.gds", allow.duplicate=T)

### common SNP.dt
  snp.dt.common <- data.table(chr=seqGetData(genofile.common, "chromosome"),
                        pos=seqGetData(genofile.common, "position"),
                        nAlleles=seqGetData(genofile.common, "$num_allele"),
                        id=seqGetData(genofile.common, "variant.id"))
  snp.dt.common <- snp.dt.common[nAlleles==2]
  seqSetFilter(genofile.common, snp.dt.common$id)

  snp.dt.common[,af:=seqGetData(genofile.common, "annotation/info/AF")[[2]]]

### open GDS for private SNPs (SNAPE)
  genofile.priv <- seqOpen("/scratch/aob2x/dest/dest.PoolSeq.SNAPE.001.50.ann.gds", allow.duplicate=T)

### file list
  priv.fn <- list.files("/scratch/aob2x/dest/privateSNPs/", pattern="SNAPE.csv", full.names=T)

### pnps function


  pnps.func <- function(pop.n) {

    ### load list of private SNPs,
      priv.pop <- fread(priv.fn[pop.n])
      pop <- gsub(".SNAPE.csv", "", last(tstrsplit(priv.fn[pop.n], "/")))

      message(pop)

    ### tabulate basic characteristics of MAF, AC, etc. for private SNPs
      message("get private")
      # priv.pop <- priv.pop[1:500]

      seqSetFilter(genofile.priv, variant.id=priv.pop$id, sample.id=pop)
      tmp.dp <- seqGetData(genofile.priv, "annotation/format/DP")
      priv.pop[,dp:=t(tmp.dp$data)[,1]]

      tmp <- seqGetData(genofile.priv, "annotation/info/ANN")
      len1 <- tmp$length
      len2 <- tmp$data

      snp.dt1 <- data.table(len=rep(len1, times=len1),
                            ann=len2,
                            id=rep(priv.pop$id, times=len1))

      # Extracting data between the 2nd and third | symbol
        snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]

      # Collapsing additional annotations to original SNP vector length
        snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                              list(variant.id=id)]

        snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
        priv.pop[,ann:=snp.dt1.an$col]

  ### get common matched control SNPs
    message("get control")
    seqResetFilter(genofile.common)
    seqSetFilter(genofile.common, sample.id=pop, variant.id=snp.dt.common$id)

    common.pop <- snp.dt.common
    tmp.dp <- seqGetData(genofile.common, "annotation/format/DP")
    common.pop[,dp:=t(tmp.dp$data)[,1]]

    tmp.ac <- seqGetData(genofile.common, "annotation/format/AD")
    common.pop[,ac:=t(tmp.ac$data)[,1]]

  ### get control SNPs; this matches on CHR, total read depth and alternate count
    control.SNPs <- merge(priv.pop[,c("chr", "dp", "ac"),with=F],
                          common.pop,
                          by=c("chr", "dp", "ac"), allow.cartesian=T)


    setkey(control.SNPs, id)
    control.SNPs <- control.SNPs[!duplicated(control.SNPs)]

    seqSetFilter(genofile.common, sample.id=pop, variant.id=control.SNPs$id)
    tmp <- seqGetData(genofile.common, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(control.SNPs$id, times=len1))

    # Extracting data between the 2nd and third | symbol
      snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]

    # Collapsing additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                            list(variant.id=id)]

      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
      control.SNPs[,ann:=snp.dt1.an$col]
      control.SNPs[,pop:=pop]

  ### tabulate pnps
    message("save")
    priv.pop[,class:="private"]
    control.SNPs[,class:="common"]
    m <- rbind(priv.pop, control.SNPs, fill=T)


    out.fn <- gsub("SNAPE.csv", "pnps.Rdata", priv.fn[pop.n])
    save(m, file=out.fn)

    return(m)
  }

  m <-  pnps.func(pop.n=pop.n)

  fisher.test(table(m[grepl("synonymous_variant|missense_variant", ann)][!grepl("splice", ann)]$ann,
                    m[grepl("synonymous_variant|missense_variant", ann)][!grepl("splice", ann)]$class))


  table(m[grepl("synonymous_variant|missense_variant", ann)][!grepl("splice", ann)]$ann,
                    m[grepl("synonymous_variant|missense_variant", ann)][!grepl("splice", ann)]$class)
