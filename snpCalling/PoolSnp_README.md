
### SNP calling with PoolSNP based on SYNC file format

Call SNPs on joined SYNC file based on (1) minimum allele count (across all samples) or (2) minimum allele frequency (across all samples). Note, that only positions with a proportion of missing data <= than the miss-frac threshold are retained. The max-cov argument is only used for the VCF header.

```bash
python ~/DEST/PoolSNP4Sync/PoolSnp.py \
--sync test_masked.sync.gz \
--min-cov 4 \
--max-cov 0.7 \
--min-count 4 \
--min-freq 0.01 \
--miss-frac 0.5 \
--names sample1,sample2 \
> test.vcf
```
