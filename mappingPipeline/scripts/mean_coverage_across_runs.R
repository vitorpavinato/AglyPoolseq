#!/usr/bin/env Rscript

#############################################################
# Access sequencing consistency across technical replicates #
#############################################################

# Each library, that corresponds to pooled sequecing of 5 individuals that share the same
# barcode, were sequenced 4-5x in different days. The libraries were the same, but they were 
# loaded and sequenced in different runs.

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES

## Remove last features
rm(list=ls())
ls()

## Receive external arguments
args = commandArgs(TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("ERROR: Need to supply <coverage_file_path> <sample_name> <output_dir> as positional arguments.n", call.=FALSE)
}

# point each argument to an object
input=args[1]
sample=args[2]
output=args[3]

#input="/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output"
#sample="MN_BIO4_S1"
#output="/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output"

# prepare the prefix of each replicate
prefix=c("_140711", "_140722", "_140725", "_140730", "_AddRun")


## Program
filenames_coverage <- paste0(input, "/", sample, prefix, "/",
                             sample, prefix, "_", "meandepth", ".txt")

# column names
column_names <- c("rname",  "startpos", "endpos", "numreads",
                  "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")

# read.table of each replicate
if (file_test("-f", filenames_coverage[5])){
  datalist_coverage <- lapply(filenames_coverage, 
                              function(x){read.table(file= x, header=F, col.names = column_names, na.strings = "NA")})
} else {
  datalist_coverage <- lapply(filenames_coverage[-c(5)], 
                              function(x){read.table(file= x, header=F, col.names = column_names, na.strings = "NA")})
}

# coverage
datalist_meancov <- sapply(datalist_coverage, function(x) x[[6]])
if (dim(x = datalist_meancov)[2] == 5)
{
  colnames(datalist_meancov) <- paste0("run", prefix)
  rownames(datalist_meancov) <- datalist_coverage[[1]]$rname
} else {
  colnames(datalist_meancov) <- paste0("run", prefix[-c(5)])
  rownames(datalist_meancov) <- datalist_coverage[[1]]$rname
}

# mean depth
datalist_meandepth <- sapply(datalist_coverage, function(x) x[[7]])
if (dim(x = datalist_meandepth)[2] == 5)
{
  colnames(datalist_meandepth) <- paste0("run", prefix)
  rownames(datalist_meandepth) <- datalist_coverage[[1]]$rname
} else {
  colnames(datalist_meandepth) <- paste0("run", prefix[-c(5)])
  rownames(datalist_meandepth) <- datalist_coverage[[1]]$rname
}

## prepare to output files

output_folder <- paste0(output, "/", sample)

# check if the folder exists
if (!file_test("-d", output_folder)){
  dir.create(file.path(output_folder))
}

# From the help of the pairs() function
# Print the absolute correlation between columns
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.4/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# plot to summarize the data
# Mean coverage
pdf(file = paste0(output_folder, "/", sample, "_meancov_pairwise_plot.pdf"))
pairs(datalist_meancov, lower.panel = panel.smooth, upper.panel = panel.cor, gap=0, row1attop=FALSE)
dev.off()

# Mean depth
datalist_log10meandepth <- apply(datalist_meandepth, 2, log10)
datalist_log10meandepth <- datalist_log10meandepth[is.finite(datalist_log10meandepth[,1]), ]

# log10 scale
pdf(file = paste0(output_folder, "/", sample, "_log10meandepth_pairwise_plot.pdf"))
pairs(datalist_log10meandepth, lower.panel = panel.smooth, upper.panel = panel.cor, gap=0, row1attop=FALSE)
dev.off()

# normal scale
pdf(file = paste0(output_folder, "/", sample, "_meandepth_pairwise_plot.pdf"))
pairs(datalist_meandepth, lower.panel = panel.smooth, upper.panel = panel.cor, gap=0, row1attop=FALSE)
dev.off()

# summary statistics
# from https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
MatVar <- function(x, dim = 1, ...) {
  if(dim == 1){
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  } else if (dim == 2) {
    rowSums((t(x) - colMeans(x, ...))^2, ...)/(dim(x)[1] - 1)
  } else stop("Please enter valid dimension")
}

# Mean coverage
row_means_meancov <- rowMeans(datalist_meancov)
row_std_meancov <- sqrt(MatVar(datalist_meancov, 1))
row_means_stds_meancov <- cbind(row_means_meancov, row_std_meancov)
write.table(row_means_stds_meancov, 
            file = paste0(output_folder, "/", sample, "_meancov_scaffolds_meanstd.txt"),
            quote = FALSE)
  
col_means_meancov <- colMeans(datalist_meancov)
col_std_meancov <- sqrt(MatVar(datalist_meancov, 2))
col_means_stds_meancov <- rbind(col_means_meancov, col_std_meancov)
write.table(col_means_stds_meancov, 
            file = paste0(output_folder, "/", sample, "_meancov_runs_meanstd.txt"),
            quote = FALSE)

# Mean depth
row_means_meandepth <- rowMeans(datalist_meandepth)
row_std_meandepth <- sqrt(MatVar(datalist_meandepth, 1))
row_means_stds_meandepth <- cbind(row_means_meandepth, row_std_meandepth)
write.table(row_means_stds_meandepth, 
            file = paste0(output_folder, "/", sample, "_meandepth_scaffolds_meanstd.txt"),
            quote = FALSE)

col_means_meandepth <- colMeans(datalist_meandepth)
col_std_meandepth <- sqrt(MatVar(datalist_meandepth, 2))
col_means_stds_meandepth <- rbind(col_means_meandepth, col_std_meandepth)
write.table(col_means_stds_meandepth, 
            file = paste0(output_folder, "/", sample, "_meandepth_runs_meanstd.txt"),
            quote = FALSE)


cat("\n Coverage Analysis Finished\n\n")
