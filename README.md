# DEST
Scripts for mapping and quality control of DEST dataset

## getSimContaminationLevel.sh
Runs samtools idxstats on mapped bam files to extract out the number of reads mapping to sim and mel.

## populationInfo directory
Has supplemental data from the DrosEU, DrosRTEC, and DPGP files to make a unified meta-datafile; attaches GHCND station based on lat. and long.

The basic script to generate the meta-data is `makeJointSampleInfo.R`
The otuput file is `samps.csv` file
