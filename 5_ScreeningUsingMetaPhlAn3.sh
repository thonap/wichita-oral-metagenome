#!/bin/bash

## Running MetaPhlAn3 on ancient shotgun metagenomic data
## Written by Tanvi Honap

############################################
## Modify before running ## 

## Set the path for the folder containing the analysis-ready reads for each sample
## Each file is named as $NAME.AR.fastq.gz where $NAME is the unique sample ID
## This folder also contains a file containing the sample IDs; this file can be named SampleList.txt 

path=/path/to/folder/containing/sample/files

############################################

## The following commands will be run on all sample files using a while loop

cat $path/SampleList.txt | while read NAME;

 do

## 1. CutAdapt3 is used to trim two nucleotides from each end of the reads which are most likely to contain damage in pUDG-treated libraries
## This step is optional and can be skipped

  cutadapt3 -u +2 -u -2 -o $path/$NAME.noDamage.AR.fastq.gz $path/$NAME.AR.fastq.gz


## 2. MetaPhlAn3 is used for metagenomic screening
## Use minimum read length of 30 for ancient DNA

  metaphlan  --nproc 10 --input_type fastq --read_min_len 30 --bowtie2out $path/$NAME.metaphlan3.bowtie2.bz2 -o $path/$NAME.metaphlan3.profile.txt $path/$NAME.noDamage.AR.fastq.gz

 done

## End ##