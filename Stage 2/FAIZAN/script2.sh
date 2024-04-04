#!/bin/bash

# Go to home directory and create necessary directories
 cd ~
 mkdir -p rawdata/sequence rawdata/reference QCreports Trimreports Alignmentreports sortedsequence Variant

 # Download the reference genome
 echo "Downloading Reference Genome..."
 wget -O rawdata/reference/reference.fasta https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

 # Index the reference genome
 echo "Indexing Reference Genome..."
 bwa index rawdata/reference/reference.fasta

 # Define datasets
 declare -a samples=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

 # Loop through each sample processing
 for sample in "${samples[@]}"; do
# Download URL
wget -O "rawdata/sequence/${sample}_R1.fastq.gz" "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/${sample}_R1.fastq.gz"     
wget -O "rawdata/sequence/${sample}_R2.fastq.gz" "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/${sample}_R2.fastq.gz"

# Run FastQC
 fastqc "rawdata/sequence/${sample}_R1.fastq.gz" "rawdata/sequence/${sample}_R2.fastq.gz" -o QCreports/

# Run fastp for trimming
fastp -i "rawdata/sequence/${sample}_R1.fastq.gz" -I "rawdata/sequence/${sample}_R2.fastq.gz" -o "Trimreports/${sample}_trim_R1.fastq.gz" -O "Trimreports/${sample}_trim_R2.fastq.gz"

# Align with bwa mem, convert to BAM, sort, and index
bwa mem rawdata/reference/reference.fasta "Trimreports/${sample}_trim_R1.fastq.gz" "Trimreports/${sample}_trim_R2.fastq.gz" | samtools view -Sb - | samtools sort -o "sortedsequence/${sample}_sorted.bam"
samtools index "sortedsequence/${sample}_sorted.bam"

# Variant calling
bcftools mpileup -Ou -f rawdata/reference/reference.fasta "sortedsequence/${sample}_sorted.bam" | bcftools call -mv -Ov -o "Variant/${sample}_variants.vcf"

# Convert VCF to CSV for variants
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' "Variant/${sample}_variants.vcf" > "Variant/${sample}_variant.csv"

 done

