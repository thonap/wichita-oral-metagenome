#!/bin/bash

## Script for mapping analysis-ready reads to a microbial genome
## Written by Tanvi Honap

############################################
## Modify before running ## 

## Define the path to the folder containing the microbial genome reference (FASTA file)

ref=/path/to/folder/containing/referencefasta

## Index the reference using both bwa and samtools
## If this has been done previously, comment out the following two lines 

bwa index $ref
samtools faidx $ref

## Set an output suffix which will be added to every output file 
out=microbe1

## Set the path for the folder containing the analysis-ready reads for each sample
## Each file is named as $NAME.AR.fastq.gz where $NAME is the unique sample ID
## This folder also contains a file containing the sample IDs; this file can be named SampleList.txt 

path=/path/to/folder/containing/sample/files

############################################

## The following commands will be run on all sample files using a while loop

cat $path/SampleList.txt | while read NAME;

 do

## 1. Mapping analysis-ready reads to a reference genome using BWA
##  -l 1024 disables seed for usage with aDNA reads
##  -n denotes edit distance and allows more or less mismatches; increase to 0.1 to make mapping more strict
## -q denotes mapping quality 
## -t specifies number of threads

  bwa aln -l 1024 -n 0.1 -q 37 -t 8 $ref $path/$NAME.AR.fastq.gz > $path/$NAME.$out.sai

  bwa samse $ref $path/$NAME.$out.sai $path/$NAME.AR.fastq.gz > $path/$NAME.$out.sam

## 2. Filtering, sorting, removing duplicates using SAMtools
## samtools view -bSh displays SAM file as a BAM file (b); input is SAM (S),and includes the header (h)

  samtools view -bSh $path/$NAME.$out.sam > $path/$NAME.$out.bam

## Filter out unmappped and low quality reads
## samtools view -bh -F4 displays the previous output as a BAM file (b), and includes the header (h), and alignments containing the 4 flag (0x4 segment unmapped)

  samtools view -bh -F4 $path/$NAME.$out.bam > $path/$NAME.$out.mapped.bam

## Sort alignments by leftmost coordinates

  samtools sort -o $path/$NAME.$out.mapped.sorted.bam $path/$NAME.$out.mapped.bam

## Remove PCR duplicates using DeDup

  java -jar $dedup -i $path/$NAME.$out.mapped.sorted.bam -m -o .

## Index alignment

  samtools index $path/$NAME.$out.mapped.sorted_rmdup.bam 

## 5. Generate damage patterns with mapDamage and estimate coverage statistics using Qualimap2
## The --rescale option can be used to generate a rescaled BAM file, which rescales the qual scores for likely damaged positions in the reads

  mapDamage -i $path/$NAME.$out.mapped.sorted_rmdup.bam -r $ref --rescale -q -d $path/$NAME\_MapDamage_results

   mv $path/$NAME\_MapDamage_results/Fragmisincorporation_plot.pdf $path/$NAME\_MapDamage_results/$NAME.$out.Fragmisincorporation_plot.pdf

  mv $path/$NAME\_MapDamage_results/$NAME.$out.mapped.sorted_rmdup.rescaled.bam $path/$NAME.$out.mapped.sorted_rmdup.rescaled.bam

  samtools index $path/$NAME.$out.mapped.sorted_rmdup.rescaled.bam

  $qualimap bamqc -bam $path/$NAME.$out.mapped.sorted_rmdup.rescaled.bam --java-mem-size=4G -outdir $path/$NAME\_RescaledBAM_QualimapStats -outformat pdf > /dev/null

## 5. Use bamUtils to clip first and second bases from each read to avoid damage affecting variant calling

  $bam trimBam $path/$NAME.$out.mapped.sorted_rmdup.rescaled.bam $path/$NAME.$out.mapped.sorted_rmdup.rescaled.clipped.bam --ignoreStrand -L 2 -R 2

## 6. Use samtools mpileup and VarScan2 to generate VCF file
## samtools mpileup -a option is used to generate calls at ALL sites, including zero read-depth ones

## Generate the VCF file using VarScan mpileup2cns
## --min-coverage = minimum coverage at a site to make a call 
## --min-reads2 = minimum number of reads covering the variant allele to call it as a SNP 
## --min-avg-qual = minimum quality required at a base to make a call 
## --min-var-freq = minimum percent of reads covering variant allele to call it a SNP
## --min-freq-for-hom = minimum percent of reads covering variant allele to call it a homozygous SNP

  samtools mpileup -aa -f $ref $path/$NAME.$out.mapped.sorted_rmdup.rescaled.clipped.bam | java -jar $VarScan mpileup2cns --min-coverage 5 --min-reads2 3 --min-avg-qual 37 --min-var-freq 0.2 --min-freq-for-hom 0.9 --p-value 1 --strand-filter 0 --output-vcf > $path/$NAME.$out.raw.vcf

 done

## End ##