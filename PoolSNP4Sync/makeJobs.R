#module load intel/18.0 intelmpi/18.0 R/3.6.0; R


### libraries
  library(data.table)
  library(foreach)

### define working directory where
  wd="/project/berglandlab/DEST"

### get paths to bgzipped sync files
  sync_files <- system(paste("ls ",
                             wd,
                             "/dest_mapped/*/*masked.sync.gz",
                             sep=""), intern=T)

### get chromosome names and lengths
  
