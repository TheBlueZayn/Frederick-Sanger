#!/bin/bash

#make the working directory 
mkdir variant-calling
ls
mkdir data
mkdir data/fastq
mkdir data/ref
cd data/fastq/

#download the working datasets (forward and reverse strand)
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
ls -lh

#download reference strand
cd ..
cd ref
wget https://zenodo.org/records/10886725/files/Reference.fasta
ls -lh
cd ..
cd ..

#visualize the downloaded data
less data/ref/Reference.fasta
zcat data/fastq/ERR8774458_1.fastq.gz | less -s
zcat data/fastq/ERR8774458_2.fastq.gz | less -s 

#quality control to fastq files
fastqc data/fastq/*.fastq.gz
#download html files 
ls -lh data/fastq

#multiqc ro aggregate the fastqc reports
mkdir multiqc
multiqc data/fastq/*_fastqc.zip -o multiqc/multiqc_data

#trimming using fastp
mkdir filtered_reads
#sudo apt-get install fastp
fastp -i data/fastq/ERR8774458_1.fastq.gz -I data/fastq/ERR8774458_2.fastq.gz -o filtered_reads/ERR8774458_1.trim.fastq -O filtered_reads/ERR8774458_2.trim.fastq -q 30
mv fastp.html filtered_reads/fastp.html
mv fastp.json filtered_reads/fastp.json

#fastq to trimmed reads
fastqc filtered_reads/*.trim.fastq
ls filtered_reads/

#genome mapping 
#1. index to ref genome
ls data/ref
bwa index data/ref/Reference.fasta
ls data/ref
#2. aligment of reads with reference genome using bwa mem
mkdir SAMS
#sudo apt-get install bwa
bwa mem data/ref/Reference.fasta filtered_reads/ERR8774458_1.trim.fastq filtered_reads/ERR8774458_2.trim.fastq > SAMS/ERR8774458.sam
ls -lh SAMS

#convert sam to bam file
mkdir BAMS
samtools view -S -b SAMS/ERR8774458.sam > BAMS/ERR8774458.aligned.bam

#sort and index bam file 
samtools sort BAMS/ERR8774458.aligned.bam -o BAMS/ERR8774458.aligned.sorted.bam
bwa index BAMS/ERR8774458.aligned.sorted.bam
#now the .aligned.sorted.bam is indexed also

#statistical analysis 
samtools flagstat BAMS/ERR8774458.aligned.sorted.bam
#NO DUPLICATES so removing duplicate reads steps will be neglected

#variant calling using bcftools
mkdir BCFresults
bcftools mpileup -O b -o BCFresults/ERR8774458_raw.bcf -f data/ref/Reference.fasta BAMS/ERR8774458.aligned.sorted.bam
ls -lh BCFresults/

#detect SNVs in VCF format
mkdir VCFresults
bcftools call -m -v -o VCFresults/ERR8774458_variants.vcf BCFresults/ERR8774458_raw.bcf
ls -lh VCFresults/
less VCFresults/ERR8774458_variants.vcf 

#filter SNVs output in vcf format 
vcfutils.pl varFilter VCFresults/ERR8774458_variants.vcf  > VCFresults/ERR8774458_final_variants.vcf
ls -lh VCFresults
less VCFresults/ERR8774458_variants.vcf
#how many variants in the vcf file 
echo SNVs Count are $(grep -v "#" VCFresults/ERR8774458_final_variants.vcf | wc -l)

#how many snps and indels in the vcf file 
echo SNPs count are $(bcftools view -v snps VCFresults/ERR8774458_final_variants.vcf|grep -v -c '^#')
echo Indels count are $(bcftools view -v indels VCFresults/ERR8774458_final_variants.vcf|grep -v -c '^#') 
