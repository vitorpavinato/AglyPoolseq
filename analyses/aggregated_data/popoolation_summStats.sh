#!bin/bash

check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2"
    exit $2
  fi
  echo "Checked step $1"
}

### USAGE: sh popoolation_summStats.sh <<inputDir>> <<outputDir>> <<poolname>> <poolsize>>

#######
#### CONFIGURATION
#######

### Global Configuration (only if you need to specify another NSL_TABLE)
POPOOLATIONPATH=/users/PAS1554/vitorpavinato/local/popoolation/1.2.2
CODON_TABLE=$POPOOLATIONPATH/syn-nonsyn/codon-table.txt
NSL_TABLE=$POPOOLATIONPATH/syn-nonsyn/nsl_p1.txt

### Specific Configuration
GTF_CDSs=/fs/scratch/PAS1715/aphidpool/ref/raw/A_glycines_b4_r2.1_CDSs.gtf
GTF_EXONS=/fs/scratch/PAS1715/aphidpool/ref/raw/A_glycines_b4_r2.1_exons.gtf
REF=/fs/scratch/PAS1715/aphidpool/ref/raw/A_glycines_b4_r2.1.fasta

### BED FILE OF RETAINED SNPS
BED_SNPs=/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/popoolation_sumstats/snps_ordered4mpileup_filtering.001.2.13Ago2021.bed

### INPUT/OUTPUT
INPUTDIR=${1}
OUTPUTDIR=${2}
POOLNAME=${3}
POOLSIZE=${4}

INPUT_FILE=${INPUTDIR}/${POOLNAME}/${POOLNAME}.agly.bam
OUTPUTPATH=${OUTPUTDIR}/${POOLNAME}

### CHECK IF OUTPUT DIR EXISTS 
[ ! -d ${OUTPUTPATH} ] && mkdir -p  ${OUTPUTPATH}

### CREATE FILTERED MPILEUP FILE
echo "running samtools mpileup"
samtools mpileup --positions ${BED_SNPs} -B -Q 25 -f ${REF} ${INPUT_FILE}  > ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup
wait

check_exit_status "mpileup" $?

### RUN POPOOLATION SUMMARY STATISTICS

## RUN POPOOLATION SYN-NONSYN AT POSITION
echo "running syn-nonsyn at position"
perl $POPOOLATIONPATH/syn-nonsyn/Syn-nonsyn-at-position.pl --pileup ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --fastq-type 'sanger' --measure pi --gtf $GTF_CDSs --codon-table $CODON_TABLE --nonsyn-length-table $NSL_TABLE --output ${OUTPUTPATH}/${POOLNAME}.SynNonsyn.out --snp-output ${OUTPUTPATH}/${POOLNAME}.SynNonsyn.snp.out --region-output ${OUTPUTPATH}/${POOLNAME}.SynNonsyn.region.out --pool-size ${POOLSIZE} --min-count 2 --min-coverage 4 --max-coverage 50 --min-qual 25 
wait

check_exit_status "syn_nonsyn" $?

### RUN VARIANCE SLIDING FOR PI W_theta AND TAJIMA'S D
echo "running sliding window pi"
perl $POPOOLATIONPATH/Variance-sliding.pl --input ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --output ${OUTPUTPATH}/${POOLNAME}.window.pi --measure pi --pool-size ${POOLSIZE} --fastq-type 'sanger' --min-count 2 --min-coverage 4 --max-coverage 50 --min-covered-fraction 0 --min-qual 25 --window-size 5000 --step-size 5000
wait

echo "running sliding window theta"
perl $POPOOLATIONPATH/Variance-sliding.pl --fastq-type 'sanger' --input ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --output ${OUTPUTPATH}/${POOLNAME}.window.theta --measure theta --window-size 10000 --step-size 10000 --min-count 1 --min-coverage 4 --pool-size ${POOLSIZE} --min-covered-fraction 0 --dissable-corrections
wait

echo "running sliding window TajimaD"
perl $POPOOLATIONPATH/Variance-sliding.pl --fastq-type 'sanger' --input ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --output ${OUTPUTPATH}/${POOLNAME}.window.TjD --measure D --window-size 10000 --step-size 10000 --min-count 1 --min-coverage 4 --pool-size ${POOLSIZE} --min-covered-fraction 0 --dissable-corrections
wait

check_exit_status "windowSummerStats" $?

### RUN VARIANCE AT POSITION FOR PI W_theta AND TAJIMA'S D
echo "running exons pi"
perl $POPOOLATIONPATH/Variance-at-position.pl --fastq-type 'sanger' --pileup ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --pool-size ${POOLSIZE} --min-coverage 4 --min-count 1  --gtf $GTF_EXONS --output ${OUTPUTPATH}/${POOLNAME}.exons.pi --measure pi --snp-output ${OUTPUTPATH}/${POOLNAME}.snp.output.exons.pi --min-covered-fraction 0 --dissable-corrections
wait

echo "running exons theta"
perl $POPOOLATIONPATH/Variance-at-position.pl --fastq-type 'sanger' --pileup ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --pool-size ${POOLSIZE} --min-coverage 4 --min-count 1  --gtf $GTF_EXONS --output ${OUTPUTPATH}/${POOLNAME}.exons.theta --measure theta --snp-output ${OUTPUTPATH}/${POOLNAME}.snp.output.exons.theta --min-covered-fraction 0 --dissable-corrections
wait

echo "running exons Tajima's D"
perl $POPOOLATIONPATH/Variance-at-position.pl --fastq-type 'sanger' --pileup ${OUTPUTPATH}/${POOLNAME}.filtered.mpileup --pool-size ${POOLSIZE} --min-coverage 4 --min-count 1  --gtf $GTF_EXONS --output ${OUTPUTPATH}/${POOLNAME}.exons.TjD --measure D --snp-output ${OUTPUTPATH}/${POOLNAME}.snp.output.exons.TjD --min-covered-fraction 0 --dissable-corrections
wait

check_exit_status "geneSummerStats" $?

echo "done"
