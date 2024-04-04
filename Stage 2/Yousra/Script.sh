#!/bin/bash

#  demand user input
read -p "Input name of the gene:" name
read -p "Input url to forward read:" R1
read -p "Input url to backward read:" R2
read -p "Input url to reference fasta data:" ref


# create the required directories to store your data

mkdir -p $name/data $name/Qc_reports $name/trimmed_results $name/Results

# download data sets
wget -P $name/data $R1
wget -P $name/data $R2
wget -P $name/data $ref

# quality control

fastqc $name/data/*.fastq.gz -o $name/Qc_reports

# Multiqc to aggregate reports

multiqc $name/Qc_reports/*.fastqc.zip -o $name/Qc_reports

# perform trimming using fastp and save output in the trimmed_results directory
fastp -i $name/data/"$name"*1.fastq.gz -I $name/data/"$name"*2.fastq.gz -o $name/trimmed_results/"$name"_R1.trimmed.fastq.gz -O $name/trimmed_results/"$name"_R2.trimmed.fastq.gz --html $name/trimmed_results/"$name"_fastp.html --json $name/trimmed_results/"$name"_fastp.json

# perform indexing using bwa
bwa index $name/data/reference.fasta

#perform read alignments
bwa mem $name/data/reference.fasta $name/trimmed_results/"$name"_R1.trimmed.fastq.gz $name/trimmed_results/"$name"_R2.trimmed.fastq.gz > $name/Results/"$name".sam

# use samtools to convert sam files to bam
samtools view -S -b $name/Results/"$name".sam > $name/Results/"$name".bam

# perform sorting of bam files
samtools sort -o $name/Results/"$name".sorted.bam $name/Results/"$name".bam

# perform indexing of bam files

samtools index $name/Results/"$name".sorted.bam

# use bcftools to pileup format 
bcftools mpileup -O b -o $name/Results/"$name".bcf -f $name/data/reference.fasta $name/Results/"$name".sorted.bam

# Call variants using bcftool and create vcf file
bcftools call -m -v -o $name/Results/"$name".variants.vcf $name/Results/"$name".bcf

echo "successfully done!"
















 

 
