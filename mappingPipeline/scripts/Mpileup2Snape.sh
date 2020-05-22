#!/bin/bash

mpileup=$1
sample=$2
theta=$3
D=$4
priortype=$5
fold=$6

awk '{if (last != $1) close(last); print >> $1; last = $1}' $mpileup

for chr in {2L,2R,3L,3R,4}; do
  snape-pooled -nchr 2 -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
done

for chr in {X,Y}; do
  snape-pooled -nchr 1 -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
done

for chr in {2L,2R,3L,3R,4,X,Y}; do
  rm $chr
done

cat *-$sample-SNAPE.txt  > ${sample}_SNAPE.txt

rm *-$sample-SNAPE.txt
