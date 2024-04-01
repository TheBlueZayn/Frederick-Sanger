#!/bin/bash

#create a new working directory for Project 3 and change into it using a single line of code
mkdir Project_3 && cd Project_3

#create a new folder (data) and a sub-folder (fastq) for fastq files
mkdir data
mkdir data/fastq

#download the forward and reverse strand into the fastq folder
wget -P data/fastq https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1 && wget -P data/fastq https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1

#change the name of the files
cd data/fastq && mv 'ERR8774458_1.fastq.gz?download=1' ERR8774458_1.fastq.gz && mv 'ERR8774458_2.fastq.gz?download=1' ERR8774458_2.fastq.gz && cd ../../

#create a new folder for the reference genome and download the reference genome into it
mkdir data/ref && wget -P data/ref https://zenodo.org/records/10886725/files/Reference.fasta?download=1

#rename the refernce genome
cd data/ref && mv Reference.fasta?download=1 Reference.fasta && cd ../../

#make directory for quality control output
mkdir QC_Reports

#perform quality control using fastqc and send result to the output folder(QC_Reports) 
fastqc data/fastq/*.fastq.gz -o QC_Reports

#use multiqc to aggregate fastqc results
multiqc QC_Reports/*_fastqc.zip -o QC_Reports

#Make new directory for trimmed files and trim using fastp
mkdir Trimmed_reports && fastp -i data/fastq/ERR8774458_1.fastq.gz -I data/fastq/ERR8774458_2.fastq.gz -o Trimmed_reports/ERR8774458_1_trimmed.fastq.gz -O Trimmed_reports/ERR8774458_2_trimmed.fastq.gz 

#Index the refrerence genome using bwa
bwa index data/ref/Reference.fasta

#Make new directory for Aligned sequences 
mkdir Aligned

#map the reads using bwa
bwa mem data/ref/Reference.fasta Trimmed_reports/*.fastq.gz -o Aligned/Aligned.sam

#convert sam file to bam file
samtools view -h -b Aligned/Aligned.sam -o Aligned/Aligned.bam

#Check mapping statistics
samtools flagstat Aligned/Aligned.bam

#Filter the mapped reads using flags i.e removing reads that do not match ref seq and reads which have their mates unmapped to the ref seq
samtools view -b -F 0xc Aligned/Aligned.bam -o Aligned/Aligned.filtered.bam

#sort filtered file by name
samtools sort -@ 8 -n Aligned/Aligned.filtered.bam -o Aligned/Aligned.sorted.n.bam

# Apply fixmate
samtools fixmate -m Aligned/Aligned.sorted.n.bam Aligned/Aligned.fixmate.bam

#sort by coordinate
samtools sort -@ 8 Aligned/Aligned.fixmate.bam -o Aligned/Aligned.sorted.c.bam

#Mark then remove duplicate
samtools markdup -r -@ 8 Aligned/Aligned.sorted.c.bam Aligned/Aligned.final.bam

#Index the bam ready file (Aligned.final.bam)
samtools index Aligned/Aligned.final.bam

#Make a new directory for variant calling files
mkdir VCF

#call variants using freebayes
freebayes -f data/ref/Reference.fasta -b Aligned/Aligned.final.bam --vcf VCF/Result.vcf

#compress file using bgzip 
bgzip VCF/Result.vcf

#index the VCf file using bcftools
bcftools index VCF/Result.vcf.gz

#Perform some basic statistics
#count the number of variants identified
zgrep -v -c '^#' VCF/Result.vcf.gz

#count the number of indels and snps
bcftools view -v snps VCF/Result.vcf.gz | grep -v -c '^#'
bcftools view -v indels VCF/Result.vcf.gz | grep -v -c '^#'

#create snp and indel files in the VCF directory
bcftools view -v snps VCF/Result.vcf.gz -Oz -o VCF/snps.vcf.gz
bcftools view -v indels VCF/Result.vcf.gz -Oz -o VCF/indels.vcf.gz

#count the number of variants in snps and indel files
zgrep -v -c '^#' VCF/snps.vcf.gz
zgrep -v -c '^#' VCF/indels.vcf.gz

#filter and select variants with quality score greater or equal to 30
bcftools filter -i "QUAL>=30" VCF/Result.vcf.gz -Oz -o VCF/quality_score.vcf.gz && grep -v -c '^#' VCF/quality_score.vcf.gz



