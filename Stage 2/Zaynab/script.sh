#!/usr/bin/env bash
# Create directories for variant calling
mkdir qc_report trimmed results
samples=("ACBarrie", "Alsen", "Baxter", "Chara", "Drysdale")

for s in "${samples[@]}"; do 
# Perform quality control check on both reads, output to qc_folder folder
    fastqc "$PWD/${s}_R1.fastq.gz" "$PWD/${s}_R2.fastq.gz" -o qc_report 
# Aggregate fastqc reports     
    multiqc qc_report/*_fastqc.zip -o qc_report
# Trim faulty reads with fastp, output to trimmed folder
    fastp \
    -i "$PWD/${s}_R1.fastq.gz"\
    -I "$PWD/${s}_R2.fastq.gz"\
    -o "trimmed/${s}_R1.trim.fastq.gz"\
    -O "trimmed/${s}_R1.trim.fastq.gz"\
    --html trimmed/"${s}_fastp.html"
# Map reads to reference genome with bwa
    bwa index "$PWD/*.fasta"
# Align reads to reference genome
    bwa mem "$PWD/*.fasta" "trimmed/${s}_R1.trim.fastq.gz" "trimmed/${s}_R2.trim.fastq.gz" > "results/${s}.aligned.sam"
# Convert sam file to bam file 
    samtools view -S -b "results/${s}.aligned.sam" > "results/${s}.aligned.bam"

# Sort sequemce
    samtools sort -o "results/${s}.aligned.sorted.bam" "results/${s}.aligned.bam"
# Index sorted bam file with samtools
    samtools index "results/${s}.aligned.sorted.bam"
# Generate pileup format for bam file
    bcftools mpileup -O b -o "results/${s}.bcf" -f "$PWD/*.fasta" "results/${s}.aligned.sorted.bam"
# Identify variants using bcftools call and generately variant (vcf) file
    bcftools call -m -v -o "results/${s}.variants.vcf" "results/${s}.bcf"
done

