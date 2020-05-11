# call SNPs on joined sync file based on minimum allele count (across all samples) or minimum allele frequency (across all samples).

# Note, that only positions with a proportion of missing data <= than the miss-frac threshold are retained.

# The max-cov argument is only used for the header.

```bash
python /Users/mkapun/Documents/GitHub/DEST/PoolSNP4Sync/PoolSnp.py \
--sync testMpileup2Sync/joined_masked.sync.gz \
--min-cov 4 \
--max-cov 0.7 \
--min-count 4 \
--min-freq 0.01 \
--miss-frac 0.5 \
--names sample1,sample2 \
> testMpileup2Sync/joined.vcf
```
