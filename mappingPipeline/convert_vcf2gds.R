### libraries
  library(SeqArray)

### modify VCF header to make read faster; do this only once.
  vcf_fn <- "/mnt/ssd/DrosEU.DrosRTEC.DPGP.09042019.ref.snpEFF.vcf"
  gds_fn <- "/mnt/ssd/DrosEU.DrosRTEC.DPGP.09042019.ref.snpEFF.vcf.reheader.gds"

  h <- seqVCF_Header(vcf_fn)
  h$format
  h$format[1, 3] <- "Integer"

  seqVCF2GDS(vcf_fn, gds_fn, header=h)

### done
