# Modified scripts for downloading, mapping, and variant calling of Aphis glycines pool-seq data

## Description
> This set of scripts provides a pipeline to build wholeGenomeSync files for each population sample from raw FASTQ data and defines a Dockerfile to build a docker image which can act as a standalone tool to run the pipeline (forked from https://github.com/alanbergland/DEST).

### 0. Define working directory
```bash
wd=/fs/scratch/PAS1715/aphidpool
```

### 1. Download data from SRA
```bash
sbatch --array=2-$( wc -l < ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv ) \
${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/download_SRA.sh
```

### 2. Check that data are in the correct FASTQ format
Double check that all downloaded data are in Fastq33. Uses script from [here](https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py). </br>

```bash
sbatch ${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/check_fastq_encoding.sh
grep -v "1.8" ${wd}/fastq/qualEncodings.delim
```

### 3. Build singularity container from docker image in SLURM cluster interactively
```bash
cd ${wd}
sinteractive -n 1 -A <<project_name>>
export SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR 
module load singularity
singularity pull docker://vpavinato/aglypoolseq:latest
exit
```

### 4. Or build it with a sbatch script
```bash
cd ${wd}
sbatch ${wd}/DEST-AglyPoolseq/mappingPipeline/docker_pull_osc.sh
```

### 5. Run the singularity container across list of populations (technical replicates files)
```bash
sbatch --array=2-$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | wc -l ) \
${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/runDocker.sh
```

### 6. Aggregate BAM files from technical replicates 
This script combines BAM files from replicated runs and re-process the combined BAM by sorting reads by coordinate, removing duplicated reads and re-aligning reads around indels.
 It assumes that all technical replicated BAMs are sorted and contains reads aligned to the hologenome.
```bash
sbatch --array=2-$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | wc -l ) \ 
${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/runAggregateBams.sh
```

### 7. Run the singularity container across list of populations (aggregated BAM files)
```bash
sbatch --array=2-$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | wc -l ) \
${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/runDocker.sh
```
