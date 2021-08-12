## Agly-Poolseq
This repository was forked from [alanbergland/DEST](https://github.com/alanbergland/DEST). It contain modified scripts for mapping and quality control that were used on [DEST](https://www.biorxiv.org/content/10.1101/2021.02.01.428994v3) dataset (Drosophila over Space and Time). The pipeline and Dockerfile were modified to work with Aphis glycines genome.

### Metadata
  > `populationInfo\`: Has supplemental data from field populations.

### PoolSeq mapping pipeline
  > `mappingPipeline\`: Contains dockerized mapping pipeline. Downloads data, produces bam files, filter files, gSYNC files

### SNP calling
  > `snpCalling\`: SNP calling based on gSYNC files </br>
  > `SNAPE\`: SNP calling based on snape output

### Utility scripts
  > `utils\`: (a) download individual files; (b)

### Analysis scripts
> `analysis\`: script to reproduce the SNP-based analysis;
