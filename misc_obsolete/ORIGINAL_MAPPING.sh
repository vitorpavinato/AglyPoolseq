#!/bin/sh

#  pipeline.sh
#
#
#  Created by Martin Kapun on 10/09/15.
#

read1=$1 #The full path to Read1
read2=$2 #The full path to Read2
name=$3 #The name of the output-file, e.g. Barcelona-spring
out=$4 #The full path to the outputfolder
mapper=$5 # Either choose "bwa" or bbmap

### set working directory
BASEDIR=$(dirname $0)
cd $BASEDIR

### run FASTQC

#mkdir $out

mkdir $out/fastqc

perl scripts/FastQC.app/Contents/Resources/Java/fastqc $read1 $read2 -o $out/fastqc

# ---------- >  1) trimming and adapter removal

### trim NEBNext Ultra DNA Library Prep Kit adapter sequences

### a) install cutadapt

#cd scripts/cutadapt-1.8.3

#python setup.py build_ext -i

export PATH=$PATH:scripts/cutadapt-1.8.3/bin

### b) run cutadapt

mkdir -p $out/trimmed

cutadapt \
-q 18 \
--minimum-length 75 \
-o $out/trimmed/$name-1.fq.gz \
-p $out/trimmed/$name-2.fq.gz \
-b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
-B CAAGCAGAAGACGGCATACGAGAT \
-O 15 \
-n 3 \
$read1 \
$read2

#-b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#-b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#-b GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#-b CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
#-O 15
#-n 3

## OPTIONS:
#-b: 5' or 3' adapter
#-O 15: minimum overlap length between sequence and adapter
#-n 3: repeat three times
#-q 18: minimum base-quality 18
#--minimun-length 75: minimum read-length 75 (if one read smaller -> both fragments discarded)
#-o $out/trimmed/1_S1_L001_R1_001.fastq.gz
#-p $out/trimmed/1_S1_L001_R2_001.fastq.gz: indicate that paired data
#data/1_S1_L001_R1_001.fastq.gz: input
#data/1_S1_L001_R2_001.fastq.gz: input 2

### c) test again with fastqc

mkdir -p $out/fastqc/trimmed
#
perl scripts/FastQC.app/Contents/Resources/Java/fastqc $out/trimmed/$name-1.fq.gz $out/trimmed/$name-2.fq.gz -o $out/fastqc/trimmed

## $out/mapping reads
#
# ---------- >  2) $out/mapping reads

##  install samtools

#cd scripts/samtools-1.3
#
#make

export PATH=$PATH:scripts/samtools-1.3

### The hologenome consists of the following genomes:

#D_melanogaster_r6.12
#S_cerevisiae
#W_pipientis
#P_entomophila
#C_intestini
#A_pomorum
#G_morbifer
#P_burhodogranariea
#P_alcalifaciens
#P_rettgeri
#E_faecalis
#L_brevis
#L_plantarum

### a) make index and map reads

mkdir $out/mapping

### either map with BWA

if [[ $mapper = bwa ]]

then

### install bwa

#cd scripts/bwa-0.7.7

#make

export PATH=$PATH:scripts/bwa-0.7.15
echo 'BWA!!!!!!!!!!!'
#bwa index reference/holo_dmel_6.12.fa
bwa mem -M -t 20 reference/holo_dmel_6.12.fa $out/trimmed/$name-1.fq.gz $out/trimmed/$name-2.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > $out/mapping/$name.bam

## OPTIONS
#-M: flag all primary local alignemnts of SMALLER substrings of a read as secondary
#-T 24: use 24 cores

#-S: input is SAM
#-b: output is BAM
#-h: keep header
#-q 20: only retain reads with MQ>20
#-F 0x100: exclude secondary alignments

### or map with BBMAP

else

scripts/bbmap-35.50/bbmap.sh in1=$out/trimmed/$name-1.fq.gz in2=$out/trimmed/$name-2.fq.gz out=$out/mapping/$name.bam ref=reference/holo_dmel_6.12.fa

fi

# ---------- >  2) depuclicate reads using Picard

mkdir $out/mapping/dedup-report

## a) sort

java \
-Xmx20g \
-Dsnappy.disable=true \
-Djava.io.tmpdir=/mnt/ssd/tmp \
-jar scripts/picard-tools-1.109/SortSam.jar \
I=$out/mapping/${name}.bam \
O=$out/mapping/${name}-sort.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

## b) de-duplicate

java \
-Xmx20g \
-Dsnappy.disable=true \
-Djava.io.tmpdir=/mnt/ssd/tmp \
-jar scripts/picard-tools-1.109/MarkDuplicates.jar \
REMOVE_DUPLICATES=true \
I=$out/mapping/$name-sort.bam \
O=$out/mapping/$name-dedup.bam \
M=$out/mapping/dedup-report/$name.txt \
VALIDATION_STRINGENCY=SILENT

## c) remove unsorted BAM file

#rm $out/mapping/$name.bam
#rm $out/mapping/$name-sort.bam

## ---------- >  3) realign around indels with GATK

## a) create DICT file

java -jar scripts/picard-tools-1.109/CreateSequenceDictionary.jar \
REFERENCE=reference/holo_dmel_6.12.fa \
OUTPUT=reference/holo_dmel_6.12.dict

## b) create FASTA index file

samtools faidx reference/holo_dmel_6.12.fa

## c) add read group to BAM files

java -jar -Xmx10g -Djava.io.tmpdir=/mnt/ssd/tmp \
scripts/picard-tools-1.109/AddOrReplaceReadGroups.jar \
INPUT=$out/mapping/$name-dedup.bam \
OUTPUT=$out/mapping/$name-dedup_rg.bam \
SORT_ORDER=coordinate \
RGID=$name \
RGLB=$name \
RGPL=illumina \
RGSM=sample \
RGPU=name \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT

## d) remove BAM without ReadGroup

rm $out/mapping/$name-dedup.bam

## e) generate target List of InDel positions

mkdir $out/mapping/realign_list

java -Djava.io.tmpdir=/mnt/ssd/tmp \
-jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R reference/holo_dmel_6.12.fa \
-I $out/mapping/$name-dedup_rg.bam \
-o $out/mapping/realign_list/$name.list

## f) re-align around InDels

java -Xmx20g -Djava.io.tmpdir=/mnt/ssd/tmp \
-jar scripts/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R reference/holo_dmel_6.12.fa \
-I $out/mapping/$name-dedup_rg.bam \
-targetIntervals $out/mapping/realign_list/$name.list \
-o $out/mapping/$name-dedup_rg_InDel.bam

rm $out/mapping/$name-dedup_rg.bam

# ---------- >  4) test for simulans contamination

mkdir $out/sim_cont

samtools mpileup \
-B \
-f reference/holo_dmel_6.12.fa \
$out/mapping/$name-dedup_rg_InDel.bam | \
python scripts/python/contamination_in_pileup.py - \
reference/sim_mel-new.div \
> $out/sim_cont/$name-hist.txt

## plot sim-specific allele frequencies as histograms
echo """
AFL=read.table('$out/sim_cont/$name-hist.txt',header=F)
pdf('$out/sim_cont/$name-hist.pdf')
SimAlleles=AFL[,3]
hist(SimAlleles,breaks=100)
abline(v=mean(SimAlleles),col='red')
legend('topright',legend=paste('mean Freq:',mean(SimAlleles),sep=' '))
dev.off()
""" > $out/sim_cont/$name-hist.r

Rscript $out/sim_cont/$name-hist.r
