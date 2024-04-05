#!/bin/bash

# Define reference genome URL
reference_url="https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"

# Define a list of datasets URLs and corresponding reverse strands (paired-end sequencing)
# Each dataset should have a corresponding reverse strand file
datasets=("https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz"
          "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz"
          "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz"
          "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz"
          "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz")

reverse_datasets=("https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz"
                   "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz"
                   "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz"
                   "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz"
                   "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz")

# Directory to store downloaded datasets, FastQC results, MultiQC results, and output VCF files
dataset_dir="datasets"
fastqc_dir="fastqc_results"
multiqc_dir="multiqc_results"
vcf_dir="vcf_files"

# Create the directories if they don't exist
mkdir -p "$dataset_dir" "$fastqc_dir" "$multiqc_dir" "$vcf_dir"

# Download reference genome
wget -P "$dataset_dir" "$reference_url"

# Extract reference genome filename
reference_file="${reference_url##*/}"

# Build BWA index
bwa index "$dataset_dir/$reference_file"

# Perform quality control using FastQC for each dataset
for fastq_file in "${datasets[@]}"; do
    fastq_basename=$(basename "$fastq_file")
    fastqc "$fastq_file" -o "$fastqc_dir"
done

# Generate a MultiQC report for all FastQC results
multiqc "$fastqc_dir" -o "$multiqc_dir"

# Iterate over each dataset URL
for i in "${!datasets[@]}"; do
    # Download forward and reverse strand datasets
    forward_file="${datasets[i]}"
    reverse_file="${reverse_datasets[i]}"
    wget -P "$dataset_dir" "$forward_file" "$reverse_file"

    # Step 1: Perform read alignment using BWA for forward reads
    forward_alignment_output="${forward_file##*/}.sam"
  bwa mem "$dataset_dir/$reference_file" "$dataset_dir/${forward_file##*/}" > "$dataset_dir/$forward_alignment_output"

    # Step 2: Perform read alignment using BWA for reverse reads
    reverse_alignment_output="${reverse_file##*/}.sam"
    bwa mem "$dataset_dir/$reference_file" "$dataset_dir/${reverse_file##*/}" > "$dataset_dir/$reverse_alignment_output"

    # Step 3: Merge alignments of forward and reverse reads
    merged_alignment_output="${forward_file##*/}_merged.bam"
    samtools merge -f "$dataset_dir/$merged_alignment_output" "$dataset_dir/$forward_alignment_output" "$dataset_dir/$reverse_alignment_output"

    # Step 4: Sort and index the merged alignment file
    sorted_alignment_output="${merged_alignment_output%.bam}_sorted.bam"
    samtools sort -o "$dataset_dir/$sorted_alignment_output" "$dataset_dir/$merged_alignment_output"
    samtools index "$dataset_dir/$sorted_alignment_output"
 
    # Step 5: Variant calling using samtools mpileup and bcftools call
    variant_output="${sorted_alignment_output%.bam}.vcf"
    samtools mpileup -uf "$dataset_dir/$reference_file" "$dataset_dir/$sorted_alignment_output" | bcftools call -mv > "$vcf_dir/$variant_output"

    echo "Variants for ${forward_file##*/} and ${reverse_file##*/} are in $vcf_dir/$variant_output"

echo "variants called successfully"
done




















 

 
