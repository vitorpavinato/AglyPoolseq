# Dataset used in the analyses

### Variant calling files produced by the mapping + variant calling pipeline and used in the different analyses steps

Dictionary: 
- `technical replicates`: each library from each sample (or pool) was sequenced 4-5x (in different days/runs) in MiSeq in order to produce sufficient coverage and to allow us check the sequencing quality before all samples were totally sequenced. For this reason, we had > 1 paired-end (PE) read file for each sample;
- `aggregated data`: this refers to the combination (aggregation) of sequence data (actually alignments) of each sample after initial read processing steps (trimming, alignment, duplication removal etc). After aggregation, sorting, duplication removal and indel realignment were re-do;
- minor allele frequency (MAC);
- minor allele count (MAC);
- `Poolfstat` an `R` package for the analyses of Pool-seq data;

🗂️ 👾The datasets were produced by running the mapping pipeline with different minimun and maximum read coverage and by running the variant calling pipeline with different _BAM_ (from each pipeline run). This was done because each instance of the mapping pipeline started with different sequencing coverage.


|                                  VCF file                                                                       |            Pipeline           |  # of Samples  |  Min_cov  |  Max_cov  |  SNP Calling |  MAF   |  MAC  | Missing % |  # of SNPs |  Script                                                    |
|------------------------------------------------------------------------------------|-------------------------|------------------|-------------|------------|----------------|---------|---------|------------- |------------- |---------------------------------------------|
|  aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz                                            |  Technical replicates  |        87             |        3       |       99%   |  PoolSNP      |  0.05    |     5     |  50%          | 40,105       |  check_replicates_pca.R                      |
|  aphidpool.PoolSeq.PoolSNP.0.`x`.paramTest.ann.vcf.gz                           |  Aggregated data       |         21           |        4       |        99%   |  PoolSNP      |  0-0.2  | 5-100  |  50%         |  `values`   |  exploratory_macmaf_pnps_missing.R |
|  aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz                               | Aggregated data        |         21           |        4       |        99%   |  PoolSNP      | 0.001  |      5     |  50%         |  344,524    |  poolfstat_poolsnp_mac_mac_limits.R |
|  aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz                                 | Aggregated data        |         21           |        4       |        99%   |  PoolSNP      | 0.05    |      5     |  50%         |  262,866    |  poolfstat_poolsnp_analyses.R              |
|  aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.filtered.pools.vcf         | Aggregated data        |         21           |        4       |        99%   |  PoolSNP      | 0.05    |      5     |  5%            |  15,248     |  clustering_poolsnp_analyses.R              |
|  aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.ann.filtered.biotypes.vcf   | Aggregated data        |         21           |        4       |        99%   |  PoolSNP      | 0.05    |      5     |  5%            |  15,841     |  clustering_poolsnp_analyses.R              |

These vcf files and the associated `R scripts` will become available with the manuscript submission

#### Why so many datasets? What did each dataset informe about?

- **aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz:** this dataset was used to quality control (QC) the sequenced reads. We checked if each sequencing of technical replicates clustered together in the PCA of allele frequencies;
- **aphidpool.PoolSeq.PoolSNP.0.20.paramTest.ann.vcf.gz:** this collection of  `vcf` files contains a sample of annotated SNPs that were used to QC the best combination of MAC/MAF for SNP filtering. We checked the relationship between MAC/MAF and the number of SNPs after filtering, the percentage of missing data, and the proportion of non-synonimous to synonimous mutations p<sub>_N_</sub> / p<sub>_S_</sub>.
- **aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz:** this dataset was used to check what was the real minor allele frequency associated with the defined MAF. In theory, MAC and MAF are independent, but for our dataset with MAC=5 the lowest MAF was  > 0.001;
- **aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz:** this dataset was used for the population genomics analyses - _e.g_ F<sub>_ST_</sub>  scan; from this dataset two sets were obtained 1) with all samples and SNPs; and 2) a subset containing the best 12 samples (one each biotype/locality) and SNPs with proportion of missing data (across samples) < 25%.

#### Script, dataset and analyses

##### Script `check_replicates_pca.R`; folder `technical_replicates`: 
Dataset: `aphidpool.all.PoolSNP.05.5.07Jul2021.vcf.gz`
1. Maximum likelihood estimation of reference allele counts (ML-counts) of technical replicates SNPs;
2. Calculation of SNPs  reference allele frequencies from the ML-counts;
3. Principal componenet analyses (PCA) of reference allele frequencies data to check if different technical replicate clustered together;
4. Correlation between repliate overall proportion of missing data and sequence coveraget;
5. Correlation between repliate overall proportion of missing data and base/alignment quality;
6. Technical repliate read count data were aggregated in samples, and pilot PCA analyses with clustering were performed (with all 21 pools, and with 19 pools);
7. Correlation between aggregated read coverage % missing data and aggregated coverage;

##### Script `exploratory_macmaf_pnps_missing.R`; folder `aggregated_data`:
Dataset: `aphidpool.PoolSeq.PoolSNP.0.x.paramTest.ann.vcf.gz` x from 0 to 0.2
1. Missing data as a function of MAC and MAF;
2. p<sub>_N_</sub> / p<sub>_S_</sub> as a function of MAC and MAF;
3. Identify problematic samples to remove based on % Missing and p<sub>_N_</sub> / p<sub>_S_</sub>;
4. T<sub>_s_</sub> / T<sub>_v_</sub> as a function of MAC and MAF (one dataset with all SNPs with `min_cov=5`, and `max_cov=95%`);

##### Script `poolfstat_poolsnp_mac_mac_limits.R`; folder `aggregated_data`:
Dataset: `aphidpool.PoolSeq.PoolSNP.001.5.22Jun2021.vcf.gz`
1. Investigate the actual MAF and MAC after load `VCF` data with `Poolfstat`;

##### Script `poolfstat_poolsnp_analyses.R`; folder `aggregated_data`:
Dataset: `aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021.vcf.gz`
Sub-datasets: 
- 21 samples (all) and all SNPs that passed the thresholds (`min_cov=4`, `max_cov=99%`, `MAC=5`, `MAF=0.05`, `missing_fraction < 0.5`);
- 12 samples that had the lowest % of missing data within each class (one per biotype / per locality); SNP were filtered to contain only SNPs with `missing_fraction <= 0.25`;

Analyses:
1. Maximum likelihood estimation of reference allele ML-counts of samples;
2. Calculation of SNPs  reference allele frequencies from the ML-counts;
3. Genome scan with multi-locus F<sub>_ST_</sub> scan;
4. SNP- Expected Heterozygosity H<sub>_E_</sub>;
5. SNP- average per-sample depth of basis (ADP);
6. Global and pairwise F<sub>_ST_</sub>;
7. Partial mantel test and causal modelling;
