#!/bin/bash

# Go to home directory and create necessary directories
  cd ~
  mkdir -p dataset/fastq dataset/ref1 dataset/ref2 QCreports trimreports alignmentreports sorted variant

# Change directory from home to dataset/fastq directory and download the forward and reverse sequence 
  cd dataset/fastq && wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz

# Change to home directory and download the reference genomes
  cd ~ && wget -O dataset/ref1/ref_1.fasta https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta && wget -O dataset/ref2/ref_2.fasta https://zenodo.org/records/10886725/files/Reference.fasta

# Index the reference genomes
  bwa index dataset/ref1/ref_1.fasta && bwa index dataset/ref2/ref_2.fasta

# Create a loop script
  nano loop.sh
  
# Define datasets
        samples=("ACBarrie" 
                 "Alsen" 
                 "Baxter" 
                 "Chara" 
                 "Drysdale"
                 )

	# Loop through each sample processing
	    for sample in "${samples[@]}"; do

	# Run FastQC
		  fastqc "dataset/fastq/${sample}_R1.fastq.gz" "dataset/fastq/${sample}_R2.fastq.gz" -o ~/QCreports/

	# Run fastp for trimming
		  fastp -i "dataset/fastq/${sample}_R1.fastq.gz" -I "dataset/fastq/${sample}_R2.fastq.gz" -o "trimreports/${sample}_trim_R1.fastq.gz" -O "trimreports/${sample}_trim_R2.fastq.gz"

	# Align the reference genome
		  bwa mem dataset/ref2/ref_2.fasta "trimreports/${sample}_trim_R1.fastq.gz" "trimreports/${sample}_trim_R2.fastq.gz" > alignmentreports/${sample}_align.sam

	# Convert the sam file to a bam
		  samtools view -S -b alignmentreports/${sample}_align.sam > alignmentreports/${sample}_align.bam

	# Sort bam files
		  samtools sort alignmentreports/${sample}_align.bam -o sorted/${sample}_Sortedalign.bam

	# Index the sorted bam files
		  samtools index sorted/${sample}_Sortedalign.bam

		  # Variant calling
		  bcftools mpileup -Ou -f dataset/ref2/ref_2.fasta "sorted/${sample}_Sortedalign.bam" |bcftools call -mv -Ov -o "variant/${sample}_variants.vcf"

	  done

# Define datasets
        samples=("ERR8774458")

	# Loop through each sample processing
	 for sample in "${samples[@]}"; do

		 # Run FastQC
		   fastqc "dataset/fastq/${Sample}_1.fastq.gz" "dataset/fastq/${Sample}_2.fastq.gz" -o ~/QCreports/

		 # Run fastp for trimming
		     fastp -i "dataset/fastq/${Sample}_1.fastq.gz" -I "dataset/fastq/${Sample}_2.fastq.gz" -o "trimreports/${Sample}_trim_1.fastq.gz" -O "trimreports/${Sample}_trim_2.fastq.gz"
		 
   # Align to reference genome
		   bwa mem dataset/ref1/ref_1.fasta "trimreports/${Sample}_trim_1.fastq.gz" "trimreports/${Sample}_trim_2.fastq.gz" > alignmentreports/${Sample}_align.sam

		 # Convert the sam file to a bam
		   samtools view -S -b alignmentreports/${Sample}_align.sam > alignmentreports/${Sample}_align.bam

		 # Sort bam files
		   samtools sort alignmentreports/${Sample}_align.bam -o sorted/${Sample}_sortedalign.bam

		 # Index the sorted bam files
		   samtools index sorted/${Sample}_sortedalign.bam

		 # Variant calling
		  bcftools mpileup -Ou -f dataset/ref1/ref_1.fasta "sorted/${Sample}_sortedalign.bam" |bcftools call -mv -Ov -o "variant/${Sample}_variants.vcf"

	  done

# Give permission to execute
  chmod u+x loop.sh

# Execute loop.sh
  ./loop.sh
