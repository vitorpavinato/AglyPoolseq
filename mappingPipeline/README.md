# Modified scripts to download, map, call polymorphism in pooled sequencing data-sets for Aphis glycines

## Description
> This set of scripts provides a pipeline to build wholeGenomeSync files for each population sample from raw FASTQ data and defines a Dockerfile to build a docker image which can act as a standalone tool to run the pipeline.

### 0. Define working directory
```bash
wd=/fs/scratch/PAS1715/aphidpool
```

### 1. Download data from SRA (specify 72 hour time limit)
```bash
sbatch --array=1-$( wc -l < ${wd}/DEST-AglyPoolseq/populationInfo/labPools2.csv ) \
${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/downloadSRA.sh
```

### 2. Check that data are in the correct FASTQ format
Double check that all downloaded data are in Fastq33. Uses script from [here](https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py). </br>

```bash
sbatch ${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/check_fastq_encoding.sh
grep -v "1.8" ${wd}/fastq/qualEncodings.delim
```

### 3. Build singularity container from docker image in Ohio Super Computer (OSC)
```bash
cd ${wd}
sinteractive -n 1 -A <<project_name>>
export SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR 
module load singularity
singularity pull docker://vpavinato/aglypoolseq:latest
exit
```

### 4. Run the singularity container across list of populations
```bash
sbatch --array=2-$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | wc -l ) \
${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/runDocker.sh
```
