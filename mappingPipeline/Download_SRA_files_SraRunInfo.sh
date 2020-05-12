
IFS=','
i=1

### read columns split by comma
while read Run ReleaseDate LoadDate spots bases spots_with_mates avgLength size_MB AssemblyName download_path Experiment LibraryName LibraryStrategy LibrarySelection LibrarySource LibraryLayout InsertSize InsertDev Platform Model SRAStudy BioProject Study_Pubmed_id ProjectID Sample BioSample SampleType TaxID ScientificName SampleName g1k_pop_code source g1k_analysis_group Subject_ID Sex Disease Tumor Affection_Status Analyte_Type Histological_Type Body_Site CenterName Submission dbgap_study_accession Consent RunHash ReadHash

do

## ignore header row
test $i -eq 1 && ((i=i+1)) && continue

## make data folder
mkdir -p ~/Desktop/mapping/${LibraryName}
cd ~/Desktop/mapping/${LibraryName}


/Volumes/MartinResearch1/NewInvMarkers/scripts/sratoolkit.2.10.5-mac64/bin/prefetch -O ~/Desktop/mapping/${LibraryName} ${Run}

/Volumes/MartinResearch1/NewInvMarkers/scripts/sratoolkit.2.10.5-mac64/bin/fasterq-dump \
--split-files \
--outfile ~/Desktop/mapping/${LibraryName}/${LibraryName} \
-e 10 \
~/Desktop/mapping/${LibraryName}/${Run}/${Run}.sra


gzip ~/Desktop/mapping/${LibraryName}/${LibraryName}\_1.fastq
gzip ~/Desktop/mapping/${LibraryName}/${LibraryName}\_2.fastq

done < /Users/mkapun/Documents/GitHub/DEST/populationInfo/drosEU_SraRunInfo.csv

