#!/bin/sh

#  FilterMAF.sh
#  
#
#  Created by Martin Kapun on 26.10.20.
#  

python3 /Users/mkapun/Documents/GitHub/DEST/SNAPE/FilterSNAPEByMAF.py \
--VCF /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz \
--MAF 0.05 \
| gzip > /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.MAF005.ann.vcf.gz


gunzip -c /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | head -100000 | gzip > ~/Desktop/test.vcf.gz

python3 /Users/mkapun/Documents/GitHub/DEST/utils/identifySNPsPerPop.py \
--VCF ~/Desktop/test.vcf.gz \
--MAFs 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2 \
--out ~/Desktop/test.dist


python3 /Users/mkapun/Documents/GitHub/DEST/utils/identifySNPsPerPop.py \
--VCF /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz \
--MAFs 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2 \
--out /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.SNP.dist

python3 /Users/mkapun/Documents/GitHub/DEST/utils/identifySNPsPerPop.py \
--VCF /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.PoolSNP.001.50.ann.vcf.gz \
--MAFs 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2 \
--out /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.PoolSNP.001.50.SNP.dist

python3 /Users/mkapun/Documents/GitHub/DEST/SNAPE/PrivateSNPMutations.py \
--VCF ~/Desktop/test.vcf.gz \
--private /Volumes/MartinResearch2/DEST/analyses/PrivateSNPs/test \
--MAFs 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2 \
--out ~/Desktop/test.dist

python3 /Users/mkapun/Documents/GitHub/DEST/SNAPE/PrivateSNPMutations.py \
--VCF /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz \
--private /Volumes/MartinResearch2/DEST/analyses/PrivateSNPs/scratch\ 2/aob2x/dest/privateSNPs \
--MAFs 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2 \
--out /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.mut

