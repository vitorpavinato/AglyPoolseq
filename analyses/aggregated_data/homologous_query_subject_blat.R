### Remove last features
rm(list=ls())
ls()

### LOAD BLAT RESULTS
BLAT_RESULTS <- 'results/aphid_poolseq_B4_v2_outliers/homology_v1_Ag_bt4_genes/blat_qGeneOutliers_sB4v1Prot.txt'
blat_column_names <- c('qseqid',	'sseqid',	'pident',	'align_length',	'mismatch',	'gapopen',	
                       'qstart',	'qend',	'sstart',	'send',	'evalue',	'bitscore')

blat.results <- read.table(BLAT_RESULTS,sep = '\t', na.strings = NA, col.names = blat_column_names, header = F)
head(blat.results)

### LOAD SAMTOOLS FAIDX
SAMTOOLS_FAIDX <- 'results/aphid_poolseq_B4_v2_outliers/outliers_genes_bedtools-windows_aa.fas.fai'
faix_column_names <- c('names', 'query_length','offset', 'linebases', 'linewidth')

samtools.faix <- read.table(SAMTOOLS_FAIDX, sep = '\t', na.strings = NA, col.names = faix_column_names, header = F)
samtools.faix <- samtools.faix[, c(1,2)]
head(samtools.faix)

### COMBINE DATASETS TO INCLUDE QUERY LENGTH TO THE BLAT OUTPUT
merged_data <- merge(blat.results, samtools.faix, by.x = c(1), by.y = c(1), all = TRUE)

write.table(x=merged_data, file='results/aphid_poolseq_B4_v2_outliers/homology_v1_Ag_bt4_genes/blat_qGeneOutliers_sB4v1Prot_withQlength.txt',
            quote = F, col.names = F, row.names = F, sep = '\t')

### KEEP HIGHLY HOMOLOGOUS PROTEIN SEQUENCES
homologous_proteins <- merged_data[which(merged_data$pident > 75.00 & merged_data$bitscore >= (merged_data$query_length * 0.70)), ]
#homologous_proteins <- merged_data[which(merged_data$bitscore >= (merged_data$query_length * 0.70)), ]

write.table(x=homologous_proteins, file='results/aphid_poolseq_B4_v2_outliers/homology_v1_Ag_bt4_genes/blat_qGeneOutliers_sB4v1Prot_highlyhomologous_.txt',
            quote = F, col.names = F, row.names = F, sep = '\t')
