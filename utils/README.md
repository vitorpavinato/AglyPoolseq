## some scripts to help download population specific files
```bash
wget -i <fileName>
```

### all SYNC files
```bash
ls /project/berglandlab/DEST/dest_mapped/*/*masked.sync.gz /project/berglandlab/DEST/dest_mapped/*/*/*masked.sync.gz |
    sed 's/\/project\/berglandlab\/DEST\//http:\/\/berglandlab.uvadcos.io\//g' > /scratch/aob2x/dest/DEST/utils/syncURLs.txt
```


### all SNAPE files
```bash
ls /project/berglandlab/DEST/dest_mapped/*/*SNAPE.complete.masked.sync.gz /project/berglandlab/DEST/dest_mapped/*/*/*SNAPE.complete.masked.sync.gz |
    sed 's/\/project\/berglandlab\/DEST\//http:\/\/berglandlab.uvadcos.io\//g' > /scratch/aob2x/dest/DEST/utils/snapeURLs.txt
```
