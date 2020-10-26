




















## some scripts to help download population specific files
```bash
wget -i <fileName>
```

### all SYNC files
```bash
ls /project/berglandlab/DEST/dest_mapped/*/*/*masked.sync.gz  | grep -v "SNAPE" |
    sed 's/\/project\/berglandlab\/DEST\//http:\/\/berglandlab.uvadcos.io\//g' > /scratch/aob2x/dest/DEST/utils/syncURLs.txt
```


### all SNAPE files
```bash
ls /project/berglandlab/DEST/dest_mapped/*/*/*SNAPE.complete.masked.sync.gz |
sed 's/\/project\/berglandlab\/DEST\//http:\/\/berglandlab.uvadcos.io\//g' > /scratch/aob2x/dest/DEST/utils/snapeURLs.txt
```

### 9 remade SYNC & SNAPE files
```bash
ls /project/berglandlab/DEST/dest_mapped/extra_pipeline_output/*/*masked.sync.gz  | grep -v "SNAPE" |
    sed 's/\/project\/berglandlab\/DEST\//http:\/\/berglandlab.uvadcos.io\//g' > /scratch/aob2x/dest/DEST/utils/nineExtra_syncURLs.txt

ls /project/berglandlab/DEST/dest_mapped/extra_pipeline_output/*/*SNAPE.complete.masked.sync.gz  |
    sed 's/\/project\/berglandlab\/DEST\//http:\/\/berglandlab.uvadcos.io\//g' > /scratch/aob2x/dest/DEST/utils/nineExtra_snapeURLs.txt
```
