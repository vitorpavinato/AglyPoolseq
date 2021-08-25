### Remove last features
rm(list=ls())
ls()

snpsANN          = 'vcf/aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.vcf.gz.info.tsv'
fstsnpsBED       = 'results/aphid_poolseq_B4_v2_outliers/signf_multifst.sorted.bed'

snpsANN_header <- c('CHROM', 'POS', 'REF_ALT', 'ADP', 'DP', 'NC', 'AC', 'AF', 'ANN')
fstsnpsBED_header <- c('CHROM', 'start', 'end')

## read the bed file containing the SNPs positions and annotation
snps_ann <- read.table(file = snpsANN, sep = '\t', col.names = snpsANN_header)

#snp_alleles <- do.call('rbind', strsplit(snps_ann$REF_ALT, ','))

# Look at the ANN column to learn about the annotation features
ann_features <- do.call('rbind', strsplit(snps_ann$ANN, '[|]' ))

## read the bed file containing the significant 'multi-fst snp window limits' 
fstsnps <- read.table(file = fstsnpsBED, sep = '\t', col.names = snps_header) 

dist2snp=10000
res=NULL
for(i in seq_along(fstsnps$CHROM)){
  
  fstsnps_ <- fstsnps[i, ]
  snps_ann_ <- snps_ann[snps_ann$CHROM %in% fstsnps_$CHROM, ]
  
  snps_in_window <- snps_ann_[snps_ann_$POS >= (fstsnps_$end - dist2snp) & 
                              snps_ann_$POS < (fstsnps_$end + dist2snp), ]
  
  d <- data.frame(CHROM=snps_in_window$CHROM, 
                  multfstSNP=paste0(fstsnps_$start,':', fstsnps_$end), 
                  window_size = ((fstsnps_$end + dist2snp)-(fstsnps_$end - dist2snp)), 
                  snps_in_window[,2:9])
  
  res <- rbind(res, d)
  
}

## Save the preliminar Annotation
write.table(x=res, file = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_snp_fst_ann.txt',
            quote = F, row.names = F, col.names = T)

## Extract only the first annotation of each SNP in the res table
res_annotation <- do.call('rbind', strsplit(res$ANN, '[|]' ))
res_annotation_names <- c('MUTATION', 'TYPE', 'IMPACT', 'GENE', "REGION", "FUNCTION", "NT_CHANGE", "AA_CHANGE")
res_annotation <- res_annotation[,c(1,2,3,4,6,8,10,11)]
colnames(res_annotation) <- res_annotation_names

## Add a NS or SS label for missense_variant OR synonymous_variant
fitness_effect  <- ifelse(res_annotation[,2] == "missense_variant" | res_annotation[,2] == "missense_variant&splice_region_variant", "NS", 
                                           ifelse(res_annotation[,2] == "synonymous_variant" | res_annotation[,2] == "splice_region_variant&synonymous_variant", "SS", NA))

## Combine this annotation with the rest of the res table
res_annotation <- data.frame(res[,-c(11)], res_annotation, FITNES_EFFECT=fitness_effect)

## Save the final annotation
write.table(x=res_annotation, file = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_snp_fst_final_ann.txt',
            quote = F, row.names = F, col.names = T)


## Take only the NS AND SS mutations
res_annotation_genes <- res_annotation[which(res_annotation$FITNES_EFFECT == "NS" | res_annotation$FITNES_EFFECT == "SS"), ]

## Calculat the pNpS for the genes
genes_fitness = NULL
for (i in seq_along(unique(res_annotation_genes$GENE)))
{
  
  gene_set <- res_annotation_genes[which(res_annotation_genes$GENE == unique(res_annotation_genes$GENE)[i]), ]
  ns = sum(gene_set$FITNES_EFFECT == "NS")
  ss = sum(gene_set$FITNES_EFFECT == "SS")
  pnps = ns/ss
  
  g <- data.frame(GENE=gene_set[1, 4], pNpS=pnps)
  genes_fitness <- rbind(genes_fitness, g)
  
}

genes_fitness[is.infinite(genes_fitness$pNpS), 'pNpS'] <- 1

## Save pNpS table
write.table(x=genes_fitness, file = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_snp_gene_pnps.txt',
            quote = F, row.names = F, col.names = T)









#########
## GENES CLOSE TO THE MULTILOCUS FST IDENTIFIED BY BEDTOOLS

#snpsBED          = 'results/aphid_poolseq_B4_v2_outliers/snps_ordered_sorted.bed'
#fstsnps2GenesBED = 'results/aphid_poolseq_B4_v2_outliers/outliers_genes_windows.txt'

#snps <- read.table(file = snpsBED, sep = '\t', col.names = snps_header)

## read the output from bedtools windows
#genes <- read.table(file = fstsnps2GenesBED, sep = '\t')
#genes <- genes[genes$V6 == 'gene', ]

#annotate.fst.snps <- function(fstsnps, snps, genes, dist2snp=5000, dist2gene=500)
#{
  res=NULL
  for(i in seq_along(fstsnps$CHROM)){
    
    fstsnps_ <- fstsnps[7, ]
    snps_ <- snps[snps$CHROM %in% fstsnps_$CHROM, ]
    genes_ <- genes[genes$V1 %in% fstsnps_$CHROM, ]
    
    if (dim(genes_)[1] == 0) {
      genes_ <- as.data.frame(t(rep(NA, length(genes_))))
      
    } else if (dim(genes_)[1] == 1) {
      g_ <- strsplit(genes_$V12, "=")[[1]][2]
      g_ <- strsplit(g_, ";")[[1]][1]
      genes_$V12 <- g_
      
    } else {
      g_ <- do.call('rbind', strsplit(genes_$V12, "="))
      g_ <- do.call('rbind', strsplit(g_[,2], ";"))
      genes_$V12 <- g_
    } # end of if
    
    snps_in_window <- snps_[snps_$end >= (fstsnps_$end - dist2snp) & 
                              snps_$end < (fstsnps_$end + dist2snp), ]
    
    if (all(is.na(genes_$V7))){
      close2gene <- NA
    } else if (length(genes_$V7) == 1) {
      close2gene <- ifelse(snps_in_window$end < (genes_$V7 - dist2gene) | 
                             snps_in_window$end > (genes_$V8 + dist2gene), NA, genes_$V12)
    } else {
      # remove duplicated entries
      genes_ <- genes_[!duplicated(genes_$V12), ]
      
      # reduce the list to the closes gene entry
      close2gene = NULL
      for (j in seq_along(snps_in_window$end))
      {
        c_s <- snps_in_window$end[2] - genes_$V7
        c_e <- snps_in_window$end[2] - genes_$V8
        
        w <- ifelse(c_s < 0 & c_e <0 | c_s > 0 & c_e >0, 'outside', 'within')
        
        if (any(w == 'within')){
          c_ <- genes_[which(w == 'within'), ]
        } else {
          
          p_s <- c_s[which(w == 'outside')]
          p_e <- c_e[which(w != 'outside')]
          
          if (all(p_s > 0) & all(p_e > 0))
          {
            c_ <- genes_[which.min(p_s + p_e), ]
          } else if (all(p_s < 0) & all(p_e < 0)){
            c_ <- genes_[which.max(p_s + p_e)]
          } else {
            if (length(p_s) = 2)
            {
              mid <- abs(p_s[2] + p_e[1])/2
              if (snps_in_window$end[2] - p_e[1] <= mid)
              {
                c_ <- genes_[which(p_e == p_e[1])]
              } else {
                c_ <- genes_[which(p_e == p_e[2])]
              }
                
            }
          
            p_s <- c(1000, -530)
            p_e <- c(500, -1530)
            
            
            p_s[-length(p_s)] + diff(p_s)/2
            
            p_s[which.min(test)]
            p_e[which.min(test)]
            
          }
          
        }
        
        
        #c_ <- genes_[which.min(abs(c1_ + c2_ )), ]
        #c_ <- genes_[which.min(abs(snps_in_window$end[1] - genes_$V7)), ]
        
        tmp <- ifelse(snps_in_window$end[j] < (c_$V7 - dist2gene) | 
                        snps_in_window$end[j] > (c_$V8 + dist2gene), NA, c_$V12) 
        
        close2gene <- rbind(close2gene, tmp)
      } # end of for j loop
    } # end of if
    
    d <- data.frame(CHROM=snps_in_window$CHROM, 
                    multfstSNP=paste0(fstsnps_$start,':', fstsnps_$end), 
                    window_size = ((fstsnps_$end + dist2snp)-(fstsnps_$end - dist2snp)), 
                    snps_in_window[,2:3],
                    closest_gene = close2gene)
    
    res <- rbind(res, d)
    
  } # end of for i loop
  
  return(res)
  
#} # end of a function


#annotate.fst.snps(fstsnps = fstsnps, snps = snps, genes = genes)





intralocus_he_signfst <- intralocus_he[intralocus_he$Chromosome %in% unique(high.multilocus.fst.ordered$CHRM), ]