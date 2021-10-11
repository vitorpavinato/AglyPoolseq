## Agly-Poolseq
This repository was forked from [alanbergland/DEST](https://github.com/alanbergland/DEST). It contain modified scripts for mapping and quality control that were used on [DEST](https://www.biorxiv.org/content/10.1101/2021.02.01.428994v3) dataset (Drosophila over Space and Time). The pipeline and Dockerfile were modified to work with _Aphis glycines_ [genome](https://academic.oup.com/g3journal/article/10/3/899/6026189).

### Metadata
  > `populationInfo\`: Has supplemental data from field populations.

### PoolSeq mapping pipeline
  > `mappingPipeline\`: Contains dockerized mapping pipeline. Downloads data, produces bam files, filter files, gSYNC files

### SNP calling
  > `snpCalling\`: SNP calling based on gSYNC files </br>
  > `SNAPE\`: SNP calling based on snape output

### Utility scripts
  > `utils\`: scripts for additional analyses (_e.g_ identification of population private alleles, etc ). Some scripts were used for the analyses of soybean aphid data (_e.g_  `PNPS4VCF.py` for the calculation of genome-wide $p_N$, $p_S$ and $p_N/p_S$) , some were not used. I kept these all scripts here to mirror the [alanbergland/DEST](https://github.com/alanbergland/DEST) original repository.
  

### Analysis scripts
> `analyses\`: script to reproduce the SNP-based analyses presented in the manuscript;
