# Population information & basic figures

## Description
>  This directory contains scripts to generate meta-data file for the DEST dataset. It first pulls together the meta-data files for the drosRTEC, drosEU, and dgn datasets from their respective supplemental data files. Also pulls in ghcnd data.

## File structure set up

## Make meta-data file ###
  > RUN: `makeJointSampleInfo.R`
  > Outputs `DEST/populationInfo/samps.csv`

## Make figures
  ## Map figure
  > input.txt is tab-delimited with four columns (Lat,Lon,Season[spring,fall,frost,NA] and Year) \br
