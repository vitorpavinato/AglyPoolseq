## SNP calling with SNAPE software

As an alternative to Pool-SNP, SNPs can be called using [SNAPE-pooled](https://github.com/EmanueleRaineri/snape-pooled/blob/master/README.md) (Raineri et al. 2012).
SNAPE-pooled computes the probability distribution ofor the frequency of the minor allele in a certain population, at a certain position in the genome.

This script will split the mpileup file into chromosomes, so SNAPE can be run for each chromosome with the specific chromosome number (i.e. when samples are males).
Then, SNAPE outputs are concatenated and parsed into VCF format. SNAPE-pooled parameters such as _Nucleotide Diversity_ can be changed.

The final VCF format includes all information per SNP provided by SNAPE:

```bash
RD:AD:RQ:AQ:PROB:P:MEAN
RD: Reference nucleotides number
AD: Alternative nucleotides number
RQ: Average quality of the reference nucleotides
AQ: Average quality of the alternative nucleotides
PROB: $1 - p (0)$ where $p (f)$ is the probability distribution function for the minor allele freqeuncy
P: $p (1)$
MEAN: $E (f)$ mean value of $f$
```
### REFERENCES
_SNP calling by sequencing pooled samples_ (Raineri E, Ferretti L, Esteve-Codina A, Nevado B, Heath S, Pérez-Enciso M. , BMC Bioinformatics. 2012 Sep 20;13:239. doi: 10.1186/1471-2105-13-239)
