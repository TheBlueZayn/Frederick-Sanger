#!/bin/bash

echo "VARIANT CALLING PIPELINE"
echo "Input details as required"

#INPUTS
read -p "input name of sample gene: " Name
read -p "input link to forward read: " Read1
read -p "input link to reverse read: " Read2
read -p "input link to reference fasta file: " Reference

#Make required directories
mkdir -p "$Name" "$Name/data" "$Name/data/readseq" "$Name/data/refseq" "$Name/QC_Reports" "$Name/Trimmed_reports" "$Name/Aligned_reports" "$Name/VCF"

#Download each sample gene read sequences and reference sequences
echo "Downloading forward and reverse read sequences..."
wget -O "$Name/data/readseq/${Name}_R1.fastq.gz" "$Read1"
wget -O "$Name/data/readseq/${Name}_R2.fastq.gz" "$Read2"

echo "Downloading reference fasta file..."
wget -O "$Name/data/refseq/reference.fasta" "$Reference"

echo "Data download complete."

#Perform quality control using fastqc and output to QC_Reports
echo "Checking read quality"
fastqc "$Name/data/readseq/"*.fastq.gz -o "$Name/QC_Reports"

#Aggregate fastqc reports
echo "Aggregating fastqc reports"
multiqc "$Name/QC_Reports" -o "$Name/QC_Reports"

echo "Quality control results aggregated"

#Using FastP, trim the reads
echo "Trimming faulty reads"
fastp -i "$Name/data/readseq/${Name}_R1.fastq.gz" -I "$Name/data/readseq/${Name}_R2.fastq.gz" -o "$Name/Trimmed_reports/${Name}_R1.trim.fastq.gz" -O "$Name/Trimmed_reports/${Name}_R2.trim.fastq.gz"

# Indexing the reference fasta file
echo "Indexing the reference fasta file"
bwa index ACBarrie/data/refseq/reference.fasta
echo "Reference indexing complete"

# Mapping read sequences to reference sequence
echo "Mapping read sequences to reference sequence"
bwa mem ACBarrie/data/refseq/reference.fasta ACBarrie/data/readseq/ACBarrie_R1.fastq.gz ACBarrie/data/readseq/ACBarrie_R2.fastq.gz -o ACBarrie/Aligned_reports/ACBarrie.aligned.sam
echo "Sequences successfully mapped"

# Converting sam files to bam files
echo "Converting sam files to bam files"
samtools view -h -b ACBarrie/Aligned_reports/ACBarrie.aligned.sam -o ACBarrie/Aligned_reports/ACBarrie.aligned.bam
echo "Sam files successfully converted to bam files"

# Checking mapping statistics
echo "Checking mapping statistics"
samtools flagstat ACBarrie/Aligned_reports/ACBarrie.aligned.bam

# Sorting bam file by coordinates
echo "Sorting bam file by coordinates"
samtools sort -@ 8 ACBarrie/Aligned_reports/ACBarrie.aligned.bam -o ACBarrie/Aligned_reports/ACBarrie.aligned.sorted.bam

# Indexing the sorted bam file
echo "Indexing the sorted bam file"
samtools index ACBarrie/Aligned_reports/ACBarrie.aligned.sorted.bam

# Calling variants
echo "Calling variants"
freebayes -f ACBarrie/data/refseq/reference.fasta ACBarrie/Aligned_reports/ACBarrie.aligned.sorted.bam > ACBarrie/VCF/ACBarrie_Result.vcf
echo "Variant calling successfully carried out"

# Compressing result file
echo "Compressing result file"
bgzip ACBarrie/VCF/ACBarrie_Result.vcf
echo "Result compression successful"

# Indexing the vcf file
echo "Indexing the vcf file"
bcftools index ACBarrie/VCF/ACBarrie_Result.vcf.gz

# Count number of variants identified
echo "Number of variants identified are:"
bcftools view -H ACBarrie/VCF/ACBarrie_Result.vcf.gz | wc -l

# Count number of indels
echo "Number of indels"
bcftools view -v indels ACBarrie/VCF/ACBarrie_Result.vcf.gz | wc -l

# Count number of snps
echo "Number of snps"
bcftools view -v snps ACBarrie/VCF/ACBarrie_Result.vcf.gz | wc -l

