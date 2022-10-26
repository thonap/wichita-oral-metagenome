#!/bin/bash

## Script for using AdapterRemoval2 to process shotgun metagenomic reads
## Written by Tanvi Honap

############################################
## Modify before running ## 

## Set the path for the folder containing the raw FASTQ reads (R1 and R2) for each sample
## This folder also contains a file containing the sample IDs; this file can be named SampleList.txt 

path=/path/to/folder/containing/sample/files

############################################

## The following commands will be run on all sample files using a while loop

cat $path/SampleList.txt | while read NAME;
 
 do
  
  AdapterRemoval --threads 2 --file1 $path/$NAME.R1.fastq.gz --file2 $path/$NAME.R2.fastq.gz --gzip --minalignmentlength 10 --trimns --trimqualities --minquality 30 --collapse --minlength 30 --output1 $path/$NAME.trimmed.R1.fastq.gz --output2 $path/$NAME.trimmed.R2.fastq.gz --outputcollapsed $path/$NAME.AR.fastq.gz --settings $path/$NAME.AdapterRemoval.txt --adapter1 AGATCGGAAGAGCACACGTGTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT 

## --trimns: removes stretches of Ns from 5-prime and 3-prime ends
## --trimqualities and --minquality: removes stretches of low qual bases as set by minquality option
## --collapse: merges overlapping reads into one and recalculates the quality
## --minlength: merged reads shorter than specified length are removed
## --output1 and output2: trimmed R1 and R2 files
## --outputcollapsed: final output file containing trimmed and merged (analysis-ready) reads
## --minalignmentlength: minimum length of overlap required for merging
## --settings: txt file containing parameters used and overall statistics
## --adapter1 and --adapter2: specifies Adapter1 and Adapter2 sequences

## Remove unnecessary files
  rm $path/your_output*

 done

## At the end of this script, you will have files for trimmed R1 and R2 reads and also merged reads.