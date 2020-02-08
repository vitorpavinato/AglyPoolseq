# DEST
  > Scripts for mapping and quality control of DEST dataset

## Make list of polymorphic sites
  > `add_DGN_data\`: The folder name is sort of a misnomer because the scripts in this folder identify the DGN polymorphic sites,
  > and merge with the DrosRTEC. RENAME FOLDER AND DIRECTORY PATHS WITHIN. <br/>

## Make DGN SYNC file
  > `DGN_sync\`

## getSimContaminationLevel.sh
* Runs samtools idxstats on mapped bam files to extract out the number of reads mapping to sim and mel.

## populationInfo directory
* Has supplemental data from the DrosEU, DrosRTEC, and DPGP files to make a unified meta-datafile; attaches GHCND station based on lat. and long.

* The basic script to generate the meta-data is `makeJointSampleInfo.R`
The otuput file is `samps.csv`

## get_dpgp.sh
* Download dpgp data and perform in-silico pooling to tack haplotypes into pooled sequencing dataset
