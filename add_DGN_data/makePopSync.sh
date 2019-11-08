#!/usr/bin/env bash
#
#SBATCH -J makePopSync # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/makePopSync.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/makePopSync.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab



############################
### paste per chromosome ###
############################

paste_per_chr () {
	paste -d','  ${2}/${3}*${1}*long > \
	${2}/${3}_${1}.csv
}
export -f paste_per_chr

#nohup parallel --gnu -j1 paste_per_chr ::: 2L 2R 3L 3R X ::: /mnt/spicy_2/dest/african/dpgp3_sequences_temp ::: dpgp3 &
nohup parallel --gnu -j1 paste_per_chr ::: 2L 2R 3L 3R X ::: /mnt/spicy_2/dest/african/dpgp2_sequences_temp ::: CO GA GU NG &




################################
### tack in reference genome ###
################################

w2l_ref () {
	#echo ${1}

	grep -v ">" ${1} | sed 's/\(.\)/\1\n/g' | grep -v '^$' > ${1}.long
}
export -f w2l_ref


#parallel --gnu -j1 w2l_ref ::: /mnt/spicy_2/dest/reference/chr2L.fa /mnt/spicy_2/dest/reference/chr2R.fa /mnt/spicy_2/dest/reference/chr3L.fa /mnt/spicy_2/dest/reference/chr3R.fa /mnt/spicy_2/dest/reference/chrX.fa


##########################
### csv to sync format ###
##########################

csv2sync () {


fn=${1}
#chr=$( echo ${fn} | rev | cut -f1 -d'/' | rev | sed 's/dpgp3_//g' | sed 's/.csv//g' )
chr=$( echo ${fn} | grep -oE '_2L|_2R|_3L|_3R|_X' | sed 's/_//g' )

paste -d' ' /mnt/spicy_2/dest/reference/chr${chr}.fa.long ${fn} | awk -F' ' -v chr=${chr} '
{
nN=gsub(/N/,"",$2)"\t"
nA=gsub(/A/,"",$2)"\t"
nC=gsub(/C/,"",$2)"\t"
nT=gsub(/T/,"",$2)"\t"
nG=gsub(/G/,"",$2)"\t"

nObs=nA+nC+nT+nG

if(nObs>0) {
if((nA/nObs > 0 && nA/nObs < 1) || (nC/nObs > 0 && nC/nObs < 1) || (nT/nObs > 0 && nT/nObs < 1) || (nG/nObs > 0 && nG/nObs<1)) {
print chr"\t"NR"\t"toupper($1)"\t"nN"\t"nA"\t"nC"\t"nT"\t"nG
}
}
}' > ${1}.sync
}
export -f csv2sync

#nohup parallel --gnu -j1 csv2sync ::: $( ls /mnt/spicy_2/dest/african/dpgp3_sequences_temp/dpgp3_*.csv ) &
nohup parallel --gnu -j1 csv2sync ::: $( ls /mnt/spicy_2/dest/african/dpgp2_sequences_temp/*.csv ) &
