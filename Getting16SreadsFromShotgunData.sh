#!/bin/bash

## Script for acquiring bacterial 16S rRNA gene sequences from shotgun metagenomic data 
## Written by Tanvi Honap

############################################
## Modify before running ## 

## Define the path to the folder containing the GreenGenes Bowtie2 indexes (gg_ref files) as well as the 97_otus.fasta file
## This folder should contain the 97_otus.fasta file available from QIIME: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
## This folder should also contain Bowtie2 indices for the 97_otus.fasta file.
## To generate these indices for the first time, run $ bowtie2-build 97_otus.fasta gg_ref
## gg_ref is the index prefix used here; can be changed as desired. 

ref=/path/to/folder/containing/GreenGenes/database/files

## Set the path for the folder containing the analysis-ready reads for each sample
## Each file is named as $NAME.AR.fastq.gz where $NAME is the unique sample ID
## This folder also contains a file containing the sample IDs; this file can be named SampleList.txt 

path=/path/to/folder/containing/sample/files

############################################

## The following commands will be run on all sample files using a while loop

cat $path/SampleList.txt | while read NAME;
 
 do
    
## Map analysis-ready reads to the Greengenes database using Bowtie2
## -p specifies number of threads
## --no-unal specifies that unmapped reads should not be output

  bowtie2 -p 5 -x $ref/gg_ref -U $path/$NAME.AR.fastq.gz --no-unal -S $path/$NAME.AR.sam

## Convert to BAM file and sort it using SAMtools
## Display SAM file as a BAM file (b); input is SAM (S), include the header (h), use 5 threads (@)

  samtools view -bSh -@ 5 $path/$NAME.AR.sam > $path/$NAME.AR.bam

  samtools sort -o $path/$NAME.AR.sorted.bam $path/$NAME.AR.bam 

## Remove duplicates using DeDup

  $dedup -i $path/$NAME.AR.sorted.bam -m -o .

## Dedup automatically adds "_rmdup" to the end of the output BAM file, for e.g. $NAME.AR.sorted_rmdup.bam 

## Generate a FASTA file 
## In the BAM file, column 1 is the name of the read and column 10 is the actual sequence of the read
## Use a grep command to grep all lines which do not begin with @
## Using awk, write out a ">" symbol, followed by name of read - this will be the FASTA header
## Then on a newline, print the sequence

  samtools view $path/$NAME.AR.sorted_rmdup.bam | grep -v "^@"|awk '{print ">"$1"\n"$10}' > $path/$NAME.fa

## Use an awk command to rename the FASTA header lines
## We want to have the name of the sample in the header line followed by a number for the read
## The number will be assigned starting from 1 for the first read and so on... 

  awk '/^>/{print ">'"${NAME}"'_"++i; next} {print}' < $path/$NAME.fa > $path/$NAME_16S.fa

## Add the 16S sequences from the sample to a MasterFile containing 16S sequences from all samples together
  cat $path/$NAME_16S.fa >> $path/All_samples_16S.fa
 
## Remove unnecessary files
  rm $path/$NAME.AR.sam
  rm $path/$NAME.AR.bam
  rm $path/$NAME.AR.sorted.bam
  rm $path/$NAME.fa

 done

## At the end of this script, you will have a FASTA file called All_samples_16S.fa which contains those reads mapping to the GreenGenes database of 16S rRNA gene sequences across all samples. 
## This file will be used for OTU-picking.  