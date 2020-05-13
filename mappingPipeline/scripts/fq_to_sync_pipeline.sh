#!/bin/bash

read1=$1
read2=$2
sample=$3
output=$4
threads=$5

if [ -f "$read1" ]; then
  if [ -f "$read2" ]; then

    mkdir -p $output/$sample/

    mkdir $output/$sample/${sample}_fastqc
    mkdir $output/$sample/${sample}_fastqc/trimmed
    fastqc $read1 $read2 -o $output/$sample/${sample}_fastqc

    cutadapt \
    -q 18 \
    --minimum-length 75 \
    -o $output/$sample/${sample}_trimmed1.fq.gz \
    -p $output/$sample/${sample}_trimmed2.fq.gz \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    $read1 $read2

    fastqc $output/$sample/${sample}_trimmed1.fq.gz $output/$sample/${sample}_trimmed2.fq.gz -o $output/$sample/${sample}_fastqc/trimmed

    #Automatically uses all available cores
    bbmerge.sh in1=$output/$sample/${sample}_trimmed1.fq.gz in2=$output/$sample/${sample}_trimmed2.fq.gz out=$output/$sample/${sample}_merged.fq.gz outu1=$output/$sample/${sample}_1_un.fq.gz outu2=$output/$sample/${sample}_2_un.fq.gz

    rm $output/$sample/${sample}_trimmed*

    bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_dmel_6.12.fa $output/$sample/${sample}_1_un.fq.gz $output/$sample/${sample}_2_un.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}_merged_un.bam

    rm $output/$sample/${sample}_1_un.fq.gz
    rm $output/$sample/${sample}_2_un.fq.gz

    bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_dmel_6.12.fa $output/$sample/${sample}_merged.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}_merged.bam

    rm $output/$sample/${sample}_merged.fq.gz

    java -jar $PICARD MergeSamFiles I=$output/$sample/${sample}_merged.bam I=$output/$sample/${sample}_merged_un.bam SO=coordinate USE_THREADING=true O=$output/$sample/${sample}_sorted_merged.bam

    rm $output/$sample/${sample}_merged.bam
    rm $output/$sample/${sample}_merged_un.bam

    java -jar $PICARD MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=$output/$sample/${sample}_sorted_merged.bam \
    O=$output/$sample/${sample}_dedup.bam \
    M=$output/$sample/${sample}_mark_duplicates_report.txt \
    VALIDATION_STRINGENCY=SILENT

    rm $output/$sample/${sample}_sorted_merged.bam

    samtools index $output/$sample/${sample}_dedup.bam

    java -jar $GATK -T RealignerTargetCreator \
    -nt $threads \
    -R /opt/hologenome/holo_dmel_6.12.fa \
    -I $output/$sample/${sample}_dedup.bam \
    -o $output/$sample/${sample}_hologenome.intervals

    java -jar $GATK \
    -T IndelRealigner \
    -R /opt/hologenome/holo_dmel_6.12.fa \
    -I $output/$sample/${sample}_dedup.bam \
    -targetIntervals $output/$sample/${sample}_hologenome.intervals \
    -o $output/$sample/${sample}_contaminated_realigned.bam

    rm $output/$sample/${sample}_dedup.bam*

    # samtools index $output/$sample/${sample}_contaminated_realigned.bam

    #Number of reads mapping to simulans and mel
    # grep -v "sim_" $output/$sample/${sample}_${sample}_original_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > $output/$sample/${sample}_num_mel.txt
    # grep "sim_" $output/$sample/${sample}_${sample}_original_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > $output/$sample/${sample}_num_sim.txt

    #Filter out the simulans contaminants
    mel_chromosomes="2L 2R 3L 3R 4 X Y mitochondrion_genome"
    sim_chromosomes="sim_2L sim_2R sim_3L sim_3R sim_4 sim_X sim_mtDNA"

    samtools view -@ $threads $output/$sample/${sample}_contaminated_realigned.bam $mel_chromosomes -b > $output/$sample/${sample}_mel.bam
    samtools view -@ $threads $output/$sample/${sample}_contaminated_realigned.bam $sim_chromosomes -b > $output/$sample/${sample}_sim.bam

    # python /opt/DEST/mappingPipeline/scripts/fix_bam.py --contaminated $output/$sample/${sample}_contaminated_realigned.bam --detect $output/$sample/${sample}_mel_and_sim.sam --prefix sim_ --output $output/$sample/${sample}_filtered

    mv $output/$sample/${sample}_contaminated_realigned.bam  $output/$sample/${sample}_original.bam
    rm $output/$sample/${sample}_contaminated_realigned.bai

    # samtools index $output/$sample/${sample}_mel.bam
    # samtools idxstats $output/$sample/${sample}_mel.bam > $output/$sample/${sample}_mel_idxstats.txt
    samtools mpileup $output/$sample/${sample}_mel.bam -f /opt/hologenome/raw/D_melanogaster_r6.12.fasta > $output/$sample/${sample}_mel_mpileup.txt

    gzip $output/$sample/${sample}_mel_mpileup.txt

    python3 /opt/DEST/mappingPipeline/scripts/Mpileup2Sync.py --mpileup $output/$sample/${sample}_mel_mpileup.txt --ref /opt/hologenome/raw/D_melanogaster_r6.12.fasta > $output/$sample/${sample}_genomewide.sync

    samtools bgzip $output/$sample/${sample}_genomewide.sync

    for contig in {2L,2R,3L,3R,4,X,Y,mitochondrion_genome}; do
        python3 /opt/DEST/PoolSNP4Sync/max-cov.py --sync $output/$sample/${sample}_genomewide.sync --contig $contig --cutoff 0.95 --out $output/$sample/${sample}_data_${contig}.cov
    done

  fi
fi
