##########################################################
#             ADD MULTILOCUS-FST VALUES TO THE           # 
#                  SNP ANNOTATION TABLE                  #
##########################################################

## Vitor Pavinato
## correapavinato.1@osu.edu
## CFAES

###
###
### ---- SETUP ANALYZES ENVIRONMENT ----
###
###

# Recover R-renv environment
setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
#renv::snapshot()

# Remove last features
rm(list=ls())
ls()

###
###
### ---- UPLOAD FILES ----
###
###

## Load table containing the SNPs positions, annotation and summary stats
snps_ann_path = 'results/aggregated_data/minmaxcov_4_99/diversity_poolsnp/combined_summary_statistics_annotation_table'

snps_ann <- read.table(file = snps_ann_path, header = T, sep = "\t")

## Load the bed file containing the significant 'multilocus-fst snp window limits' 
snp_fstBED_path      = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst.bed'
snp_fstBED_header <- c('CHROM', 'start', 'end')

fstsnps <- read.table(file = snp_fstBED_path, sep = '\t', col.names = snp_fstBED_header) 
dim(fstsnps)

## Load the file with multilocus-fst values
snp_fstVALUES_path = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_multilocus_fst_values.txt'

fstVALUES <- read.table(file = snp_fstVALUES_path, header = T) 
dim(fstVALUES)

## Combine the Multilocus-FST coordinates with values

# Sanity check if the two file are in the same order
sum(paste0(fstsnps$CHROM,"_",fstsnps$end) == paste0(fstVALUES$CHRM,"_",fstVALUES$POS))

fstsnps['multilocus_fst'] <- fstVALUES$FST

## Keep only SNPs in the SNP annotation file that felt within the X Kbp limits around the multilocus-FST
dist2snp=10000
significant_snps_window = NULL
for(i in seq_along(fstsnps$CHROM)){
  
  fstsnps_ <- fstsnps[i, ]
  snps_ann_ <- snps_ann[snps_ann$CHROM %in% fstsnps_$CHROM, ]
  
  snps_in_window <- snps_ann_[snps_ann_$POS >= (fstsnps_$end - dist2snp) & 
                              snps_ann_$POS < (fstsnps_$end + dist2snp), ]
  
  d <- data.frame(snps_in_window, 
                  multfstSNP=paste0(fstsnps_$start,':', fstsnps_$end), 
                  window_size = ((fstsnps_$end + dist2snp)-(fstsnps_$end - dist2snp)), 
                  multilocus_fst=fstsnps_$multilocus_fst)
  
  significant_snps_window <- rbind(significant_snps_window, d)
  
}

head(significant_snps_window)

## Save the preliminar Annotation
write.table(x=significant_snps_window, file = 'results/aggregated_data/minmaxcov_4_99/multilocusfst_poolsnp/multilocus_fst_signfsnp_ann_sumstats.txt',
            quote = F, row.names = F, col.names = T)

## Take only the NS AND SS mutations
#res_annotation_genes <- res_annotation[which(res_annotation$FITNES_EFFECT == "NS" | res_annotation$FITNES_EFFECT == "SS"), ]
#
### Calculat the pNpS for the genes
#genes_fitness = NULL
#for (i in seq_along(unique(res_annotation_genes$GENE)))
#{
#  
#  gene_set <- res_annotation_genes[which(res_annotation_genes$GENE == unique(res_annotation_genes$GENE)[i]), ]
#  ns = sum(gene_set$FITNES_EFFECT == "NS")
#  ss = sum(gene_set$FITNES_EFFECT == "SS")
#  pnps = ns/ss
#  
#  g <- data.frame(GENE=gene_set[1, 4], pNpS=pnps)
#  genes_fitness <- rbind(genes_fitness, g)
#  
#}
#
#genes_fitness[is.infinite(genes_fitness$pNpS), 'pNpS'] <- 1
#
### Save pNpS table
#write.table(x=genes_fitness, file = 'results/aggregated_data/minmaxcov_4_99/poolfstat_poolsnp/signf_snp_gene_pnps.txt',
#            quote = F, row.names = F, col.names = T)


##REMOVE
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