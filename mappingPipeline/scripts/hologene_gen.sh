#!/bin/bash

################################################
# Script to construct reference "hologenome"   #
# for D. melanogaster and associated microbes  #
# modified from Casey M. Bergman (Manchester)  #
# 	Martin Kapun 17/10/2016		               #
################################################

################################################
#    Script for Aphid glycines "hologenome"    #
#    Modified by Vitor Pavinato 09/03/2021     #
################################################

# download genomes for microbes known to be associated with A. glycines from NCBI A_gly_BT4_v2.1
# unzip file (the fasta-file contains all chromosomes)
# clean-up header a bit to improve readability

OUTDIR=/opt/hologenome/raw
mkdir -p ${OUTDIR}
cd ${OUTDIR}

## A. glycines B4 genome v2.1
#curl -O https://ftp.ncbi.nih.gov/genomes/genbank/invertebrate/Aphis_glycines/latest_assembly_versions/GCA_009928515.1_A_gly_BT4_v2.1/GCA_009928515.1_A_gly_BT4_v2.1_genomic.fna.gz
#gunzip -c GCA_009928515.1_A_gly_BT4_v2.1_genomic.fna.gz | sed 's/ A.*//g' > A_glycines_b4_r2.1.fasta
# 941 contigs

# A. glycines B4 genome v2.1 - from ZENODO 10.5281/zenodo.3453468
curl -O https://zenodo.org/record/3453468/files/Aphis_glycines_4.v2.1.scaffolds.fa.gz
gunzip -c Aphis_glycines_4.v2.1.scaffolds.fa.gz > A_glycines_b4_r2.1.fasta
# 941 contigs

# download genomes for microbes known to be associated with A. glycines from NCBI
# also download genomes from Mathers (2019) publication, which contain assembled genomes for
# Buchenera aphidicola, Wolbachia, and for the mitochondrial genome of A. glycines
# All were recovered from the assembly (Mathers 2019).
# if necessary, combine multi-contig draft assemblies into single fasta file
# rename fasta headers to include species name

# Buchnera aphidicola - Aphis glycines - from NCBI
curl -O https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Buchnera_aphidicola/all_assembly_versions/GCF_001280225.1_ASM128022v1/GCF_001280225.1_ASM128022v1_genomic.fna.gz
gunzip -c GCF_001280225.1_ASM128022v1_genomic.fna.gz | sed 's/>* .*/_B_aphidicola/g' > B_aphidicola_ragly.1.fasta
# 3 contigs

# Arsenophonus endosymbiont - Aphis craccivora - from NCBI
curl -O https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Arsenophonus_endosymbiont_of_Aphis_craccivora/all_assembly_versions/GCF_013460135.1_ASM1346013v1/GCF_013460135.1_ASM1346013v1_genomic.fna.gz
gunzip -c GCF_013460135.1_ASM1346013v1_genomic.fna.gz | sed 's/>* .*/_Arsenophonus_endosymbiont/g' > Arsenophonus_endosymbiont_racra.1.fasta
# 2 contigs

# Candidatus Regiella insecticola 5.15 - from NCBI
curl -O https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Candidatus_Regiella_insecticola/all_assembly_versions/GCF_000284655.1_ASM28465v1/GCF_000284655.1_ASM28465v1_genomic.fna.gz
gunzip -c GCF_000284655.1_ASM28465v1_genomic.fna.gz | sed 's/>* .*/_C_R_insecticola/g' > C_R_insecticola_r5.15.fasta
# 562 contigs

# Wolbachia pipientis - from NCBI
#curl -O https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Wolbachia_pipientis/all_assembly_versions/GCF_000242415.1_ASM24241v2/GO_TO_CURRENT_VERSION/GCF_000242415.2_ASM24241v3_genomic.fna.gz
#gunzip -c GCF_000242415.2_ASM24241v3_genomic.fna.gz | sed '/^>/s//>W_pipientis_/g' | sed 's/ W.*//g' > W_pipientis_r3.fasta
# 156 contigs

# Wolbachia - Aphis glycines - from ZENODO 10.5281/zenodo.3453468
curl -O https://zenodo.org/record/3453468/files/Aphis_glycines_4_Wolbachia_v1.fa
cat Aphis_glycines_4_Wolbachia_v1.fa | sed '/_pilon/s/$/_Wolbachia_agly/' > Wolbachia_endosymbiont_ragly.1.fasta
# 9 contig

# Aphis glycines mitochondrial DNA - from ZENODO 10.5281/zenodo.3453468
curl -O https://zenodo.org/record/3453468/files/Aphis_glycines_4_Mitochondrial_genome.fa
cat Aphis_glycines_4_Mitochondrial_genome.fa | sed '/_pilon/s/$/_Mitochondrial_genome_agly/' > Mitochondrial_genome_ragly.1.fasta
# 1 contig


## remove raw gzipped files

rm *gz
rm *.fa

# gzip *

# create A. glycines biotype 4 "Hologenome" from individual species fasta files
cat ${OUTDIR}/A_glycines_b4_r2.1.fasta ${OUTDIR}/Mitochondrial_genome_ragly.1.fasta ${OUTDIR}/B_aphidicola_ragly.1.fasta ${OUTDIR}/Arsenophonus_endosymbiont_racra.1.fasta ${OUTDIR}/C_R_insecticola_r5.15.fasta ${OUTDIR}/Wolbachia_endosymbiont_ragly.1.fasta > ../holo_agly_b4_r2.1.fa
