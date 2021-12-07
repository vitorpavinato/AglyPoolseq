## README

Repeatmasker needs to be installed, e.g. using homebrew

```bash
brew install repeatmasker
```

First create and change into directory

```bash
mkdir ~/AglyPoolseq/RepeatMasker/ref
cd ~/AglyPoolseq/RepeatMasker/ref
```

Download reference genome
```bash
curl -O https://zenodo.org/record/3453468/files/Aphis_glycines_4.v2.1.scaffolds.fa.gz
gunzip -c Aphis_glycines_4.v2.1.scaffolds.fa.gz > A_glycines_b4_r2.1.fasta
```
repeat mask reference genome
```bash
RepeatMasker \
-pa 20 \
--lib ~/AglyPoolseq/RepeatMasker/ref/repeatmasker.insecta.fa \
--gff \
--q \
~/AglyPoolseq/RepeatMasker/ref/A_glycines_b4_r2.1.fasta
```
