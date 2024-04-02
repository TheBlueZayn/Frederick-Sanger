#!/usr/bin/env bash

# Create directories for variant calling
mkdir qc_report
mkdir trimmed 
mkdir results
SAMPLES=("ACBarrie", "Alsen", "Baxter", "Chara", "Drysdale")

for s in "${samples[@]}";

do 
    fastqc *_fastq.gz -o qc_report 
    multiqc qc_report/*_fastqc.zip -o qc_report
    fastp \
    -i "$PWD/${s}_R1.fastq.gz"\
    -I "${s}_R2.fastq.gz"\
    -o "trimmed/${s}_R1.trim.fastq.gz"\
    -O "trimmed/${s}_R1.trim.fastq.gz"\
    --html trimmed/"${s}_fastp.html"

    bwa index "$PWD/*.fasta"

    bwa mem "$PWD/*.fasta" "trimmed/${s}_R1.trim.fastq.gz" "trimmed/${s}_R2.trim.fastq.gz" > "results/${s}.aligned.sam"

    samtools view -S -b "results/${s}.aligned.sam" > "results/${s}.aligned.bam"


    samtools sort -o "results/${s}.aligned.sorted.bam" "results/${s}.aligned.bam"

    samtools index "results/${s}.aligned.sorted.bam"


    bcftools mpileup -O b -o "results/${s}.bcf" -f "$PWD/*.fasta" "results/${s}.aligned.sorted.bam"

    bcftools call -m -v -o "results/${s}.variants.vcf" "results/${s}.bcf"

    done



    bgzip "results/${s}.variants.vcf"

    bcftools "results/${s}.variants.vcf.gz"
    done


# Download reads and reference genome
wget -P data/fastq https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget -P data/ref https://zenodo.org/records/10886725/files/Reference.fasta

# Perform quality control check on both reads, output to qc_folder folder
mkdir qc_report
fastqc data/fastq/*.fastq.gz -o qc_report

# Aggregate fastqc reports 
multiqc qc_report/*_fastqc.zip -o qc_report

# Trim faulty reads with fastp, output to trimmed folder
mkdir trimmed 
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --html trimmed/ERR8774458.fastp.html --json --thread 2 -i data/fastq/ERR8774458_1.fastq.gz -I data/fastq/ERR8774458_2.fastq.gz -o trimmed/ERR8774458_1.fastq.gz -O trimmed/ERR8774458_2.fastq.gz

# Map reads to reference genome with bwa
mkdir -p results/sam results/bam results/bcf results/vcf
# Index reference with bwa
bwa index data/ref/Reference.fasta
# Align reads to reference genome
bwa mem data/ref/Reference.fasta trimmed/ERR8774458_1.fastq.gz trimmed/ERR8774458_2.fastq.gz > results/sam/ERR8774458.aligned.sam

# Convert sam file to bam file 
samtools view -S -b results/sam/ERR8774458.aligned.sam > results/bam/ERR8774458.aligned.bam

# Sort sequemce
samtools sort -o results/bam/ERR8774458.aligned.sorted.bam results/bam/ERR8774458.aligned.bam

# Index sorted bam file with samtools
samtools index results/bam/ERR8774458.aligned.sorted.bam

#Generate statistics about the sorted bam file
samtools flagstat results/bam/ERR8774458.aligned.sorted.bam

# Generate pileup format for bam file
bcftools mpileup -O b -o results/bcf/ERR8774458_raw.bcf -f data/ref/Reference.fasta results/bam/ERR8774458.aligned.sorted.bam

# Identify SNVs using bcftools call and generately variant (vcf) file
bcftools call --ploidy 1 -m -v -o results/vcf/ERR8774458.variants.vcf results/bcf/ERR8774458_raw.bcf

# Compress variant file
bgzip results/vcf/ERR8774458.variants.vcf
# index variant file
bcftools index results/vcf/ERR8774458.variants.vcf.gz

# Counting all variants and storing value in a text file
echo "Number of variants: $(zgrep -v -c "^#" results/vcf/ERR8774458.variants.vcf.gz)" >> variants.txt 

# Counting single nucleotides polymorphisms (SNPs) and Insertions/deletions (indels)
echo "Number of snps: $(bcftools view -v snps results/vcf/ERR8774458.variants.vcf.gz|grep -v -c "^#")" >> variants.txt
echo "Number of indels: $(bcftools view -v indels results/vcf/ERR8774458.variants.vcf.gz|grep -v -c "^#")" >> variants.txt

# Filter snps and indels into separate files
bcftools view -v snps results/vcf/ERR8774458.variants.vcf.gz -Oz -o results/vcf/snps.vcf.gz
bcftools view -v indels results/vcf/ERR8774458.variants.vcf.gz -Oz -o results/vcf/indels.vcf.gz


# Convert the vcf file to  csv file for further analysis and visualisation
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' results/vcf/ERR8774458.variants.vcf.gz > variant.csv
# Convert snps and indel vcf files to csv files 
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' results/vcf/snps.vcf.gz > snps.csv && bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' results/vcf/indels.vcf.gz > indels.csv


