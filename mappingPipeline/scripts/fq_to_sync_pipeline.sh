#!/bin/bash

read1=$1
read2=$2
sample=$3
output=$4
threads=$5

if [ -f "$read1" ]; then
  if [ -f "$read2" ]; then

    mkdir -p $output/$sample/

    mkdir $output/$sample/fastqc
    mkdir $output/$sample/fastqc/trimmed
    fastqc $read1 $read2 -o $output/$sample/fastqc

    cutadapt \
    -q 18 \
    --minimum-length 75 \
    -o $output/$sample/trimmed1.fq.gz \
    -p $output/$sample/trimmed2.fq.gz \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    $read1 $read2

    fastqc $output/$sample/trimmed1.fq.gz $output/$sample/trimmed2.fq.gz -o $output/$sample/fastqc/trimmed

    #Automatically uses all available cores
    bbmerge.sh in1=$output/$sample/trimmed1.fq.gz in2=$output/$sample/trimmed2.fq.gz out=$output/$sample/merged.fq.gz outu1=$output/$sample/1_un.fq.gz outu2=$output/$sample/2_un.fq.gz

    rm $output/$sample/trimmed*

    bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_dmel_6.12.fa $output/$sample/1_un.fq.gz $output/$sample/2_un.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/merged_un.bam

    rm $output/$sample/1_un.fq.gz
    rm $output/$sample/2_un.fq.gz

    bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_dmel_6.12.fa $output/$sample/merged.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/merged.bam

    rm $output/$sample/merged.fq.gz

    java -jar $PICARD MergeSamFiles I=$output/$sample/merged.bam I=$output/$sample/merged_un.bam SO=coordinate USE_THREADING=true O=$output/$sample/sorted_merged.bam

    rm $output/$sample/merged.bam
    rm $output/$sample/merged_un.bam

    java -jar $PICARD MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=$output/$sample/sorted_merged.bam \
    O=$output/$sample/dedup.bam \
    M=$output/$sample/mark_duplicates_report.txt \
    VALIDATION_STRINGENCY=SILENT

    rm $output/$sample/sorted_merged.bam

    samtools index $output/$sample/dedup.bam

    java -jar $GATK -T RealignerTargetCreator \
    -nt $threads \
    -R /opt/hologenome/holo_dmel_6.12.fa \
    -I $output/$sample/dedup.bam \
    -o $output/$sample/hologenome.intervals

    java -jar $GATK \
    -T IndelRealigner \
    -R /opt/hologenome/holo_dmel_6.12.fa \
    -I $output/$sample/dedup.bam \
    -targetIntervals $output/$sample/hologenome.intervals \
    -o $output/$sample/contaminated_realigned.bam

    rm $output/$sample/dedup.bam*

    # samtools index $output/$sample/contaminated_realigned.bam

    #Number of reads mapping to simulans and mel
    # grep -v "sim_" $output/$sample/${sample}_original_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > $output/$sample/num_mel.txt
    # grep "sim_" $output/$sample/${sample}_original_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > $output/$sample/num_sim.txt

    #Filter out the simulans contaminants
    mel_chromosomes="2L 2R 3L 3R 4 X Y mitochondrion_genome"
    sim_chromosomes="sim_2L sim_2R sim_3L sim_3R sim_4 sim_X sim_mtDNA"

    samtools view -@ $threads $output/$sample/contaminated_realigned.bam $mel_chromosomes -b > $output/$sample/mel.bam
    samtools view -@ $threads $output/$sample/contaminated_realigned.bam $sim_chromosomes -b > $output/$sample/sim.bam

    # python /opt/DEST/mappingPipeline/scripts/fix_bam.py --contaminated $output/$sample/contaminated_realigned.bam --detect $output/$sample/mel_and_sim.sam --prefix sim_ --output $output/$sample/filtered

    mv $output/$sample/contaminated_realigned.bam  $output/$sample/original.bam

    # samtools index $output/$sample/mel.bam
    # samtools idxstats $output/$sample/mel.bam > $output/$sample/mel_idxstats.txt
    samtools mpileup $output/$sample/mel.bam -f /opt/hologenome/raw/D_melanogaster_r6.12.fasta > $output/$sample/mel_mpileup.txt

    python3 /opt/DEST/mappingPipeline/scripts/Mpileup2Sync.py --mpileup $output/$sample/mel_mpileup.txt --ref /opt/hologenome/raw/D_melanogaster_r6.12.fasta > $output/$sample/genomewide.sync

    for contig in {2L,2R,3L,3R,4,X,Y,mitochondrion_genome}; do
        python3 /opt/DEST/PoolSNP4Sync/max-cov.py --sync $output/$sample/genomewide.sync --contig $contig --cutoff 0.95 --out $output/$sample/data_${contig}.cov
    done

  fi
fi
