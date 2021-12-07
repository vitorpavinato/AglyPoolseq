module load parallel


renameFun () {
  #param="UA_Pir_14_26;UA_Pyr_14_26"
  param=${1}
  old=$( echo ${param} | cut -f1 -d'.' )
  new=$( echo ${param} | cut -f2 -d'.' )

  echo ${old}
  echo ${new}

  for name in /project/berglandlab/DEST/dest_mapped/pipeline_output/${old}/${old}*; do
      #nam=/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pir_14_26/UA_Pir_14_26.original.bam

      newname=$( echo ${name} | sed "s/\/${old}\./\/$new\./g" )
      echo ${name}
      echo ${newname}

      mv $name $newname
  done

  mv \
  /project/berglandlab/DEST/dest_mapped/pipeline_output/${old} \
  /project/berglandlab/DEST/dest_mapped/pipeline_output/${new}

}
export -f renameFun

parallel renameFun ::: UA_Pir_14_26.UA_Pyr_14_26 UA_Pir_15_21.UA_Pyr_15_21 UA_Pyr_16_48.UA_Pir_16_48


mv \
/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pyr_14_26/UA_Pir_14_26_fastqc \
/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pyr_14_26/UA_Pyr_14_26_fastqc


mv \
/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pyr_15_21/UA_Pir_15_21_fastqc \
/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pyr_15_21/UA_Pyr_15_21_fastqc


mv \
/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pir_16_48/UA_Pyr_16_48_fastqc \
/project/berglandlab/DEST/dest_mapped/pipeline_output/UA_Pir_16_48/UA_Pir_16_48_fastqc
