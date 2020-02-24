# Whole Genome Sync methods

## DGN:
### getting DGN data to this stage requires 1-3, 15 of `DEST/add_DGN_data/README.md`
## Parse DGN data ###
  ### 0. Download all DGN data
  > Needs a tab delimited file with jobID, prefix, path to DGN bz2 file: DEST/add_DGN_data/dgn.list <br/>
  > Note that job 4 will fail. Why? Because 4 is the fourth line on DGN website for the DSPR. I don't think that we need to include that one.<br/>
  > RUN: `sbatch --array=1-8 /scratch/aob2x/dest/DEST/add_DGN_data/downloadDGN.sh`<br/>
  > OUT: /scratch/aob2x/dest/dgn/rawData<br/>

  ### 1. Unpack
  > Each tarball is a bit different so the unpack scripts are different for each 1-8 (minus 4), from above. <br/>
  > `sbatch --array=1-8 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/>
  ~~> These are only for debugging: <br/> ~~
  ~~> DPGP2 RUN & DONE: `sbatch --array=1 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> DPGP3 RUN & DONE: `sbatch --array=2 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> DGRP RUN & DONE: `sbatch --array=3 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> CLARK RUN & DONE: `sbatch --array=5 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> NUZHDIN RUN & DONE: `sbatch --array=6 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> POOL RUN: `sbatch --array=7 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> SIMULANS RUN & DONE: `sbatch --array=8 /scratch/aob2x/dest/DEST/add_DGN_data/unpack.sh` <br/> ~~
  ~~> OUT: /scratch/aob2x/dest/dgn/wideData/<br/> ~~

  ### 2. Wide to long
  > should be 4725 jobs <br/>
  > RUN: `cd /scratch/aob2x/dest/dgn/wideData/; ls * | tr '\t' '\n' | awk '{print NR"\t"$0}' > /scratch/aob2x/dest/dgn/dgn_wideFiles.delim` <br/>
  > RUN: `sbatch --array=1-$( tail -n1 /scratch/aob2x/dest/dgn/dgn_wideFiles.delim | cut -f1 ) /scratch/aob2x/dest/DEST/add_DGN_data/wide2long.sh` <br/>
  > A quick check to make sure things look good:
  > `w2l_check.R`

  ### 3. Format reference genome
  > RUN: `sbatch /scratch/aob2x/dest/DEST/add_DGN_data/getRefGenome.sh` <br/>
  > Check that length of reference geomes are the same as the DGP long files: <br/>
  > `wc -l /scratch/aob2x/dest/referenceGenome/r5/2L.long` <br/>
  > `wc -l /scratch/aob2x/dest/dgn/longData/B_B04_Chr2L.seq.long` <br/>
  > `wc -l /scratch/aob2x/dest/referenceGenome/r5/3L.long` <br/>
  > `wc -l /scratch/aob2x/dest/dgn/longData/B_B04_Chr3L.seq.long` <br/>
  > the reference geome is one line longer than the DGN data but that is due to a trailing empty line. No prob. <br/>

  ### 15. BGzip compress whole genome, per population SYNC datasets
  ### borrows code from step #3
  > RUN: `nJobs=$( tail -n1 /scratch/aob2x/dest/dgn/pops.delim | cut -f1 )` <br/>
  > RUN: `sbatch --array=1-${nJobs} /scratch/aob2x/dest/DEST/add_DGN_data/makePopGenomeSync.sh` <br/>


## Pooled data:
### TBD





## Tools
### 1. Extract polymorphic sites for chr%in%pop
> RUN `sbatch /scratch/aob2x/dest/DEST/wholeGenomeSync/getPolymorphicSites.sh`
