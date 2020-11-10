#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 0.5 hours
##SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083


# sbatch /scratch/aob2x/dest/DEST/misc_obsolete/move.fixedSamps.pipelineOutput.sh
# sacct -j 18771065


#ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_fall
#ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_spring
#ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_fall
#ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_spring

rm -fr /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_fall
rm -fr /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_spring
rm -fr /project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_fall
rm -fr /project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_spring

rsync -a /scratch/aob2x/tempOut/AT_gr_12_fall /project/berglandlab/DEST/dest_mapped/pipeline_output/.
rsync -a /scratch/aob2x/tempOut/AT_gr_12_spring /project/berglandlab/DEST/dest_mapped/pipeline_output/.
rsync -a /scratch/aob2x/tempOut/ES_ba_12_fall /project/berglandlab/DEST/dest_mapped/pipeline_output/.
rsync -a /scratch/aob2x/tempOut/ES_ba_12_spring /project/berglandlab/DEST/dest_mapped/pipeline_output/.

ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_fall
ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_spring
ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_fall
ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_spring
