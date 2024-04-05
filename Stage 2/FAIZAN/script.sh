#!/usr/bin/env bash

# Go to home directory and create necessary directories
  cd ~
   mkdir -p rawdata/sequence rawdata/reference1 rawdata/reference2 QCreports Trimreports Alignmentreports sortedsequence Variant

   # Change directory from home  to sequence directory Download the forward and reverse sequence in sequence directory
    cd rawdata/sequence  && wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz

    # Download the reference genomes
    cd ~ && wget -O rawdata/reference1/reference1.fasta https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta && wget -O rawdata/reference2/reference2.fasta https://zenodo.org/records/10886725/files/Reference.fasta

     # Index the reference genomes
      bwa index rawdata/reference1/reference1.fasta &&  bwa index rawdata/reference2/reference2.fasta

       # Define datasets
        samples=("ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale")

	# Loop through each sample processing
	 for sample in "${samples[@]}"; do

		 # Run FastQC
		  fastqc "rawdata/sequence/${sample}_R1.fastq.gz" "rawdata/sequence/${sample}_R2.fastq.gz" -o ~/QCreports/

		  # Run fastp for trimming
		  fastp -i "rawdata/sequence/${sample}_R1.fastq.gz" -I "rawdata/sequence/${sample}_R2.fastq.gz" -o "Trimreports/${sample}_trim_R1.fastq.gz" -O "Trimreports/${sample}_trim_R2.fastq.gz"

		  # Align t reference genome
		  bwa mem rawdata/reference1/reference1.fasta "Trimreports/${sample}_trim_R1.fastq.gz" "Trimreports/${sample}_trim_R2.fastq.gz" > Alignmentreports/${sample}_Alignsequence.sam

		  #Convert the sam file to a bam
		  samtools view -S -b Alignmentreports/${sample}_Alignsequence.sam > Alignmentreports/${sample}_Alignsequence.bam

		  #Sort bam files
		  samtools sort Alignmentreports/${sample}_Alignsequence.bam -o sortedsequence/${sample}_SortedAlignsequence.bam

		  # Index the sorted bam files
		  samtools index sortedsequence/${sample}_SortedAlignsequence.bam

		  # Variant calling
		  bcftools mpileup -Ou -f rawdata/reference1/reference1.fasta "sortedsequence/${sample}_SortedAlignsequence.bam" | bcftools call -mv -Ov -o "Variant/${sample}_variants.vcf"

	  done

# Define datasets
        samples=("ERR8774458")

	# Loop through each sample processing
	 for sample in "${samples[@]}"; do

		 # Run FastQC
		  fastqc "rawdata/sequence/${sample}_1.fastq.gz" "rawdata/sequence/${sample}_2.fastq.gz" -o ~/QCreports/

		  # Run fastp for trimming
		  fastp -i "rawdata/sequence/${sample}_1.fastq.gz" -I "rawdata/sequence/${sample}_2.fastq.gz" -o "Trimreports/${sample}_trim_1.fastq.gz" -O "Trimreports/${sample}_trim_2.fastq.gz"

		  # Align to reference genome
		  bwa mem rawdata/reference2/reference2.fasta "Trimreports/${sample}_trim_1.fastq.gz" "Trimreports/${sample}_trim_2.fastq.gz" > Alignmentreports/${sample}_Alignsequence.sam

		  #Convert the sam file to a bam
		  samtools view -S -b Alignmentreports/${sample}_Alignsequence.sam > Alignmentreports/${sample}_Alignsequence.bam

		  #Sort bam files
		  samtools sort Alignmentreports/${sample}_Alignsequence.bam -o sortedsequence/${sample}_SortedAlignsequence.bam

		  # Index the sorted bam files
		  samtools index sortedsequence/${sample}_SortedAlignsequence.bam

		  # Variant calling
		  bcftools mpileup -Ou -f rawdata/reference2/reference2.fasta "sortedsequence/${sample}_SortedAlignsequence.bam" | bcftools call -mv -Ov -o "Variant/${sample}_variants.vcf"

	  done

