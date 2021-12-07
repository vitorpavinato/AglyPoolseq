# Scripts to calculate summary Popoolation summary statistics

## Description
>  This bash script produce a mpileup file for each sample and use this file to calculate summary statistics implemented in Popoolation

### 0. Define working directory
```bash
wd=/fs/scratch/PAS1715/aphidpool
```

### 1. Run the script
```bash
sbatch --array=2-$( cat ${wd}/AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | wc -l ) \ 
${wd}/AglyPoolseq/analyses/aggregated_data/run_popoolation.sh
```
