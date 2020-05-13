#!/bin/sh


## USAGE: sbatch SNAPE_calling.sh mel_mpileup.txt prueba 80
## Final outputs: prueba-SNAPE_formated_header.vcf


##### To be changed by the user
#module load snape-pooled/1
#module load Python/2.7.15-foss-2018b


pileup=$1
name=$2
chromosomes=$3

awk '$1=="2L" {print $0}' $pileup > 2L_$pileup
awk '$1=="2R" {print $0}' $pileup > 2R_$pileup
awk '$1=="3L" {print $0}' $pileup > 3L_$pileup
awk '$1=="3R" {print $0}' $pileup > 3R_$pileup
awk '$1=="4" {print $0}' $pileup > 4_$pileup
awk '$1=="X" {print $0}' $pileup > X_$pileup
awk '$1=="Y" {print $0}' $pileup > Y_$pileup
awk '$1=="mitochondrion_genome" {print $0}' $pileup > mitochondrion_genome_$pileup

snape-pooled -nchr $chromosomes -theta 0.005 -D 0.01 -priortype informative -fold unfolded < 2L_$pileup > 2L-$name-SNAPE.txt
snape-pooled -nchr $chromosomes -theta 0.005 -D 0.01 -priortype informative -fold unfolded < 2R_$pileup > 2R-$name-SNAPE.txt
snape-pooled -nchr $chromosomes -theta 0.005 -D 0.01 -priortype informative -fold unfolded < 3L_$pileup > 3L-$name-SNAPE.txt
snape-pooled -nchr $chromosomes -theta 0.005 -D 0.01 -priortype informative -fold unfolded < 3R_$pileup > 3R-$name-SNAPE.txt
snape-pooled -nchr $chromosomes -theta 0.005 -D 0.01 -priortype informative -fold unfolded < 4_$pileup > 4-$name-SNAPE.txt
snape-pooled -nchr $(($chromosomes/2)) -theta 0.005 -D 0.01 -priortype informative -fold unfolded < X_$pileup > X-$name-SNAPE.txt
snape-pooled -nchr $(($chromosomes/2)) -theta 0.005 -D 0.01 -priortype informative -fold unfolded < Y_$pileup > Y-$name-SNAPE.txt
snape-pooled -nchr $(($chromosomes/2)) -theta 0.005 -D 0.01 -priortype informative -fold unfolded < mitochondrion_genome_$pileup > mitochondrion_genome-$name-SNAPE.txt

cat 2L-$name-SNAPE.txt 2R-$name-SNAPE.txt 3L-$name-SNAPE.txt 3R-$name-SNAPE.txt 4-$name-SNAPE.txt X-$name-SNAPE.txt Y-$name-SNAPE.txt mitochondrion_genome-$name-SNAPE.txt > $name_SNAPE.txt

python SNAPE_to_VCF.py $name-SNAPE.txt

## STEP 4: Filter masked sites

## STEP 5: Merged individual VCFs
