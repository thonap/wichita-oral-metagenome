#!/bin/bash

## Script for mapping analysis-ready reads to the human genome
## Written by Tanvi Honap

############################################
## Modify before running ## 

## Define the path to the folder containing the human genome reference (FASTA file and bowtie2 indices)
## To generate these indices for the first time, run $ bowtie2-build hg19.fasta hg19
## hg19 is the index prefix used here; can be changed as desired. 
## Assumes that Bowtie2 indices and FASTA file have the same prefix

ref=/path/to/folder/

## Set an output suffix which will be added to every output file 
out=hg19

## Set the path for the folder containing the analysis-ready reads for each sample
## Each file is named as $NAME.AR.fastq.gz where $NAME is the unique sample ID
## This folder also contains a file containing the sample IDs; this file can be named SampleList.txt 

path=/path/to/folder/containing/sample/files

############################################

## The following commands will be run on all sample files using a while loop

cat $path/SampleList.txt | while read NAME;

 do
  
## Map analysis-ready (AR) reads to the huma genome using Bowtie2
## -p specifies number of threads
## --no-unal specifies that unmapped reads should not be output

  bowtie2 -p 5 -x $ref -U $path/$NAME.AR.fastq.gz -S $path/$NAME.$out.sam --no-unal

## Convert SAM to BAM file using SAMTools
## input is SAM (-S), display the previous output as a BAM file (-b), include the header (-h), use five threads (-@ 5)

  samtools view -bSh -@ 5 $path/$NAME.$out.sam > $path/$NAME.$out.bam

## Filter out low quality reads
## display the previous output as a BAM file (-b), include the header (-h), filter out reads with mapping quality less than 37 (-q 37)

  samtools view -@ 5 -bh -q 37 $path/$NAME.$out.bam > $path/$NAME.$out.q37.bam

## Sort alignments by leftmost coordinates

  samtools sort -o $path/$NAME.$out.q37.sorted.bam $path/$NAME.$out.q37.bam 

## Remove duplicates using DeDup
## Dedup automatically adds "_rmdup" to the end of the output BAM file, for e.g. $NAME.$out.q37.sorted_rmdup.bam

  $dedup -i $path/$NAME.$out.q37.sorted.bam -m -o .

## Index the file
  samtools index $path/$NAME.$out.q37.sorted_rmdup.bam

## Run MapDamage2 
## The --rescale option can be used to generate a rescaled BAM file, which rescales the quality scores for likely damaged positions in the reads
## -q denotes quiet

  mapDamage -i $path/$NAME.$out.q37.sorted_rmdup.bam -r $ref.fasta --rescale -q -d $path/$NAME.$out.MapDamage

## Run Qualimap to get coverage estimates
  
  qualimap bamqc -bam $path/$NAME.$out.MapDamage/$NAME.$out.q37.sorted_rmdup.rescaled.bam -outdir $path/$NAME.$out.Qualimap --java-mem-size=4G -outformat pdf > /dev/null 
  
 done

## End ##
