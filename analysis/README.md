# Dataset used for the analysis

### The dataset were produced by running the mapping pipeline with different minimun and maximum read coverage and by running the variant calling pipeline with different _BAM_ (from each pipeline run). This was done because each instance of the mapping pipeline started with different sequencing coverage.

|                                  VCF file                                               |            Source            |  # of Samples   |   Min_cov  | Max_cov  |  MAF    |  MAC                         | Missing franction | # of SNP variants |  Associated Script               |
|  aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz                    |  Technical replicates  |        87              |        3        |       99%   |  0.05    |     5                            |  50%                     | 40,105                  |  check_replicates_pca.R     |
|  aphidpool.PoolSeq.PoolSNP.0.20.paramTest.ann.vcf.gz  | Aggregated data        |         21             |        4       |        99%   |  many   | 5, 10, 15, 20, 50, 10  |  50%                     | many values         |  exploratory_macmaf_pnps_missing.R |
|  aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz      | Aggregated data        |         21             |        4       |        99%   | 0.001    |      5                           |  50%                     |  344,524               |  poolfstat_poolsnp_mac_mac_limits.R |
|  aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz        | Aggregated data        |         21             |        4       |        99%   | 0.05      |      5                           |  50%                     |  262,866               |  poolfstat_poolsnp_analysis.R              |

These vcf files and the associated `R scripts` will become available with the manuscript submission

#### Why so many datasets? What did each dataset informe about?

- **aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz:** this dataset was used to check if there were some issues with the sequencing. We checked if each sequencing of technical replicates clustered together in the PCA of allele frequencies;
- **aphidpool.PoolSeq.PoolSNP.0.20.paramTest.ann.vcf.gz:** this collection of  `vcf` files contains a sample of annotated SNPs that were used to check the best combination of MAC/MAF for SNP filtering. We checked the relationship between MAC/MAF and the number of SNPs after filtering, the percentage of missing data, and the proportion of non-synonimous to synonimous mutations (`p_N/p_S`).
- **aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz:** this dataset was used to check what was the real minor allele frequency associated with the defined MAF. In theory, MAC and MAF are independent, but for our dataset with MAC=5 the lowest MAF was  > 0.001;
- **aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz:** this dataset was used for the population genomics analysis - `F_{ST}` scan; from this dataset two sets were obtained 1) with all samples and SNPs; and 2) a subset containing the best 12 samples (one each biotype/locality) and SNPs with proportion of missing data (across samples) < 25%.
