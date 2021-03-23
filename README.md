# DEST-AglyPoolseq
  > Modified Scripts for mapping and quality control of DEST dataset. The pipeline and Dockerfile were modified to work with Aphis glycines genome.

## Metadata
  > `populationInfo\`: Has supplemental data from field and laboratory populations.

## PoolSeq mapping pipeline
  > `mappingPipeline\`: Contains dockerized mapping pipeline. Downloads data, produces bam files, filter files, gSYNC files

## Incorporate DGN data
  > `add_DGN_datda\`: Downloads, formats, lifts-over DGN data into gSYNC format. Not used in this project.

## SNP calling
  > `PoolSNP4Sync\`: SNP calling based on gSYNC files </br>
  > `SNAPE\`: SNP calling based on snape output

## Utility scripts
  > `utility\`: (a) download individual files; (b)
