#!/bin/bash

#make a new drirectory
mkdir variantcall_loop 
cd variantcall_loop

#download datasets
mkdir datasets
cd datasets
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
#download reference 
cd ..
mkdir ref
cd ref 
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta
bwa index reference.fasta
#download reference seq for ERR8774458
wget https://zenodo.org/records/10886725/files/Reference.fasta
bwa index Reference.fasta
#quality control to fastq files
cd ..
mkdir datasets/fastqc_reports 
fastqc datasets/*.fastq.gz -o datasets/fastqc_reports
 
#multiqc to aggregate
multiqc datasets/fastqc_reports/*_fastqc.zip -o datasets/multiqc_data

#make results directories
mkdir filtered_reads
mkdir filtered_reads/fastqc_trim_reports
mkdir SAMS
mkdir BAMS
mkdir BCFresults
mkdir VCFresults

#for loop, define sample
samples=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

for sample in "${samples[@]}"; do
  fastp \
    -i "datasets/${sample}_R1.fastq.gz" \
    -I "datasets/${sample}_R2.fastq.gz" \
    -o "filtered_reads/${sample}_R1.trim.fastq" \
    -O "filtered_reads/${sample}_R2.trim.fastq" \
    # Add other fastp options as needed
    # Add further processing commands here (e.g., move files)
  fastqc "filtered_reads/${sample}_R1.trim.fastq" -o "filtered_reads/fastqc_trim_reports"
  fastqc "filtered_reads/${sample}_R2.trim.fastq" -o "filtered_reads/fastqc_trim_reports"
  bwa mem ref/reference.fasta "filtered_reads/${sample}_R1.trim.fastq" "filtered_reads/${sample}_R2.trim.fastq" > SAMS/${sample}.sam
  samtools view -S -b SAMS/${sample}.sam > BAMS/${sample}.aligned.bam
  samtools sort BAMS/${sample}.aligned.bam -o BAMS/${sample}.aligned.sorted.bam
  bwa index BAMS/${sample}.aligned.sorted.bam
  bcftools mpileup -O b -o BCFresults/${sample}_raw.bcf -f ref/reference.fasta BAMS/${sample}.aligned.sorted.bam
  bcftools call -m -v -o VCFresults/${sample}_variants.vcf BCFresults/${sample}_raw.bcf
done

##analysis of ERR8774458 ##
fastp -i datasets/ERR8774458_1.fastq.gz -I datasets/ERR8774458_2.fastq.gz -o filtered_reads/ERR8774458_1.trim.fastq -O filtered_reads/ERR8774458_2.trim.fastq
fastqc filtered_reads/*.trim.fastq
bwa mem ref/Reference.fasta filtered_reads/ERR8774458_1.trim.fastq filtered_reads/ERR8774458_2.trim.fastq > SAMS/ERR8774458.sam
samtools view -S -b SAMS/ERR8774458.sam > BAMS/ERR8774458.aligned.bam
samtools sort BAMS/ERR8774458.aligned.bam -o BAMS/ERR8774458.aligned.sorted.bam
bwa index BAMS/ERR8774458.aligned.sorted.bam
bcftools mpileup -O b -o BCFresults/ERR8774458_raw.bcf -f ref/Reference.fasta BAMS/ERR8774458.aligned.sorted.bam
bcftools call -m -v -o VCFresults/ERR8774458_variants.vcf BCFresults/ERR8774458_raw.bcf
