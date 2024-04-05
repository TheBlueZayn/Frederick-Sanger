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
echo "aggregating fastqc reports"
multiqc "$Name/QC_Reports/"*_fastqc.zip -o "$Name/QC_Reports"

echo "quality control results aggregated"

#Using FastP, trim the reads
echo "Trimming faulty reads"
fastp -i "$Name/data/readseq/"*1.fastq.gz -I "$Name/data/readseq/"*2.fastq.gz -o "$Name/Trimmed_reports/$Name_R1.trim.fastq.gz" -O "$Name/Trimmed_reports/$Name_R2.trim.fastq.gz"

#Index reference data file for easch sample genome using BWA
echo "Indexing the refernce fasta file"
bwa index "$Name/data/refseq/Reference.fasta"
echo "refernce indexing complete"

#Map read sequence using BWA
echo "mappig read sequences to reference sequence"
bwa mem "$Name/data/refseq/Reference.fasta" "$Name/Trimmed_reports/"*.trim.fastq.gz -o "$Name/Aligned_reports/$Name.aligned.sam"
echo "Sequence successfully mapped" -o "$Name/Aligned_reports/$Name.aligned.sam"

#Convert sam file to bam file
echo "converting sam files to bam files"
samtools view -h -b "$Name/Aligned_reports/$Name.aligned.sam" -o "$Name/Aligned_reports/$Name.aligned.bam"
echo " Sam fils successfully converted to bam files"

#Check mapping statistics
echo "checking mapping statistics"
samtools flagstat "$Name/Aligned_reports/$Name.aligned.bam"

#Filter maped reads using flags
samtools view -b -F 0xc "$Name/Aligned_reports/$Name.aligned.bam" -o "$Name/Aligned_reports/$Name.aligned.filt.bam"

#Sort filtered files by name
samtools sort -@ 8 -n "$Name/Aligned_reports/$Name.aligned.filt.bam" -o "$Name/Aligned_reports/$Name.aligned.sorted.n.bam"

#Apply fixmate
samtools fixmate -m "$Name/Aligned_reports/$Name.aligned.sorted.n.bam" -o "$Name/Aligned_reports/$Name.aligned.fixmate.bam"

#sort by corrdinate
samtools sort -@ 8 "$Name/Aligned_reports/$Name.aligned.fixmate.bam" -o "$Name/Aligned_reports/$Name.aligned.sorted.c.bam"

#Mark then remove duplicates
samtools markdup -r -@ 8 "$Name/Aligned_reports/$Name.aligned.sorted.c.bam" "$Name/Aligned_reports/$Name.aligned.final.bam"

echo "filerting and sorting completed"

#Index the bam ready file
echo "indexing bam ready file"
samtools index "$Name/Aligned_reports/$Name.aligned.final.bam"
echo "indexing complete"

#Call variants using frrebayes
echo "Variants calling in process"
freebayes -f "$Name/data/refseq/Reference.fasta" -b "$Name/Aligned_reports/$Name.aligned.final.bam" --vcf "$Name/VCF/$Name_Result.vcf"
echo "Variant calling sucessfully carried out"

#Compress result files using bgzip
echo "compressing result file"
bgzip "$Name/VCF/$Name_Result.vcf"
echo "result compression successful"

#Index the vcf files using bcftools
bcftools index "$Name/VCF/$Name_Result.vcf.gz"

#count number of variants identified
echo "Number of variants indentified are:"
zgrep -v -c '^#' "$Name/VCF/$Name_Result.vcf.gz"

#Count number of indels and snps
echo "Number of indels"
bcftools view -v indels "$Name/VCF/$Name_Result.vcf.gz" -Oz -o "$Name/VCF/$Name_indels.vcf.gz"
echo "Number of snps"
bcftools view -v snps "$Name/VCF/$Name_Result.vcf.gz" -Oz -o "$Name/VCF/$Name_snps.vcf.gz"








