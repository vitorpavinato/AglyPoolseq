#!/bin/sh

#  Script.sh
#  
#
#  Created by Martin Kapun on 09.05.20.
#

cd /Users/mkapun/Documents/GitHub/DEST/mappingPipeline/scripts

for i in 1 2

do

## convert mpileup to SYNC file and preserve Indel positions and the coverage distribution as python objects

python Mpileup2Sync_onSteroids.py \
--mpileup testMpileup2Sync/test${i}.txt \
--ref testMpileup2Sync/ref.fa.txt \
--output testMpileup2Sync/out${i}

## create masked sync file and BED file with masked positions

python MaskSYNC.py \
--sync testMpileup2Sync/out${i}.sync.gz \
--output testMpileup2Sync/out${i} \
--indel testMpileup2Sync/out${i}.indel \
--coverage testMpileup2Sync/out${i}.cov \
--mincov 4 \
--maxcov 0.7

done

## combine the two masked SYNC files

command=paste
for i in testMpileup2Sync/out*_masked.sync.gz; do
command="$command <(gzip -cd $i )"
done
eval $command | cut -f 1-4,8 | gzip > testMpileup2Sync/joined_masked.sync.gz

## call SNPs on joined sync file

python /Users/mkapun/Documents/GitHub/DEST/PoolSNP4Sync/PoolSnp.py \
--sync testMpileup2Sync/joined_masked.sync.gz \
--min-cov 4 \
--max-cov 0.7 \
--min-count 4 \
--min-freq 0.01 \
--miss-frac 0.5 \
--names sample1,sample2 \
> testMpileup2Sync/joined.vcf
