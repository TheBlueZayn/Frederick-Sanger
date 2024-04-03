#!/bin/bash

1. # make a directory for your pipeline
 mkdir variant-calling && cd variant-calling

2. # make subdirectory for your datasets
 mkdir data
mkdir data/fastq
mkdir data/ref
cd data/fastq

3. # once in your data/fastq directory download your data
cd data/fastq
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz

4. # go to  data/ref directory and download reference data
cd ..
cd ref
wget https://zenodo.org/records/10886725/files/Reference.fasta

5. # download and run fastqc on the forward and reverse strands

sudo apt-get install fastqc
fastqc ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz

6. # download and run multiqc to aggregate data

sudo apt-get install multiqc
multiqc ERR8774458_1_fastqc.zip ERR8774458_2_fastqc.zip

7. # download fastp and trim your reads

sudo apt-get install fastp

fastp -i ERR8774458_1.fastq.gz -I ERR8774458_2.fastq.gz -o output_R1_trimmed.fastq.gz -O output_R2_trimmed.fastq.gz --detect_adapter_for_pe  --trim_poly_g --html report.html --json report.json

8. #install and use bwa for indexing data and mapping

sudo apt-get install bwa

bwa index Reference.fasta


bwa mem data/ref/Reference.fasta data/fastq/output_R1_trimmed.fastq.gz data/fastq/output_R2_trimmed.fastq.gz > aligned_reads.sam

8. # convert sam to bam file

9. # install samtools and use it to convert sam files to bam files and sort reads

sudo apt-get install samtools

samtools flagstat aligned_reads.sam

samtools view -bS aligned_reads.sam | samtools sort -o sorted_reads.bam

10. # remove duplicates and index using samtools


samtools rmdup sorted_reads.bam deduplicated_reads.bam

samtools index deduplicated_reads.bam

11 # install freebayes/bcftools and call variants

sudo apt-get install freebayes

freebayes -f Reference.fasta ~/variant-calling/deduplicated_reads.bam > ~/variant-calling/variants.vcf

12. # install tabix to zip vcf file and bcftools for analysis of the called variants

sudo apt-get install bcftools
sudo apt-get install tabix

bgzip variants.vcf

bcftools index variants.vcf.gz


12. # get the total no of variants, total no of snps and indels respectively 

zgrep -v -c '^#' variants.vcf.gz 

bcftools view -v snps variants.vcf.gz | grep -v -c '^#

bcftools view -v indels VCF/Result.vcf.gz | grep -v -c '^#'















 

 
