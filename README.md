# DEST
  > Scripts for mapping and quality control of DEST dataset

## Metadata
  > `populationInfo\`: Has supplemental data from the DrosEU, DrosRTEC, and DPGP files to make a unified meta-datafile; attaches GHCND station based on lat. and long.

## Download PoolSeq data from SRA
  > `downloadData\`: Download data from SRA using information from `populationInfo/samps.csv`

## PoolSeq mapping pipeline
  > `mappingPipeline\`: Contains dockerized mapping pipeline. Produces bam files, filter files, gSYNC files

## Incorporate DGN data
  > `add_DGN_datda\`: Downloads, formats, liftsover DGN data into gSYNC format

## SNP calling
  > `PoolSNP4Sync\`: SNP calling based on gSYNC files
  > `snape\`: SNP calling based on snape output
