#!/bin/bash

# List of URLs for the datasets
urls=(
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz"
    "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz"
    "https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta"
)

# Directory to save the downloaded files
save_dir="datasets"

# Create the directory if it doesn't exist
mkdir -p "$save_dir"

# Downloading the datasets 
for url in "${urls[@]}"; do
    # Extract filename from the URL
    filename=$(basename "$url")
    save_path="$save_dir/$filename"

    # Download the file
    wget -O "$save_path" "$url"

    echo "Downloaded $filename"
done

echo "All datasets downloaded successfully."

# Directory to save fastqc and multiqc results
qc_dir="qc_results"

# Create the directory if it doesn't exist
mkdir -p "$qc_dir"

# Run FastQC for each fastq.gz file
for file in "$save_dir"/*.fastq.gz; do
    fastqc -o "$qc_dir" "$file"
done

echo "FastQC analysis completed successfully."

# Run MultiQC to aggregate FastQC results
multiqc "$qc_dir" -o "$qc_dir"

echo "MultiQC analysis completed successfully."

# Indexing reference file
reference_file="$save_dir/reference.fasta"

# Indexing reference file
bwa index "$reference_file"

# Run Trimmomatic for each pair of raw reads
trim_dir="trimmed_reads"
mkdir -p "$trim_dir"

for file in "$save_dir"/*_R1.fastq.gz; do
    input_r1="$file"
    input_r2="${file/_R1/_R2}"
    filename_r1=$(basename "$input_r1")
    filename_r2=$(basename "$input_r2")
    output_r1="$trim_dir/trimmed_${filename_r1}"
    output_r2="$trim_dir/trimmed_${filename_r2}"

    trimmomatic PE -phred33 "$input_r1" "$input_r2" "$output_r1" /dev/null "$output_r2" /dev/null \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "Trimmed reads for $filename_r1 and $filename_r2 saved as trimmed_${filename_r1} and trimmed_${filename_r2}"
done

# Directory to save BWA alignment results
alignment_dir="alignment_results"
mkdir -p "$alignment_dir"

# Run BWA alignment for each pair of trimmed reads
for file in "$trim_dir"/*_R1.fastq.gz; do
    input_r1="$file"
    input_r2="${file/_R1/_R2}"
    filename_r1=$(basename "$input_r1")
    filename_r2=$(basename "$input_r2")
    basename_r1="${filename_r1%_R1.fastq.gz}"
    output_sam="$alignment_dir/${basename_r1}.sam"

    bwa mem -t 4 "$reference_file" "$input_r1" "$input_r2" > "$output_sam"

    echo "BWA alignment for $filename_r1 and $filename_r2 completed."
done


# Directory to save BAM files
bam_dir="bam_files"
mkdir -p "$bam_dir"

# Convert SAM to BAM for each SAM file
for sam_file in "$alignment_dir"/*.sam; do
    filename=$(basename "$sam_file")
    output_bam="$bam_dir/${filename%.sam}.bam"

    # Sort and convert SAM to BAM
    samtools view -@ 4 -bS "$sam_file" | samtools sort -@ 4 -o "$output_bam" -

    echo "Converted $filename to BAM format."
done

# Directory to save marked BAM files
marked_bam_dir="marked_bam_files"
mkdir -p "$marked_bam_dir"

# Run Picard MarkDuplicates for each BAM file
for bam_file in "$bam_dir"/*.bam; do
    output_marked_bam="$marked_bam_dir/$(basename "$bam_file" .bam)_marked.bam"
    metrics_file="$marked_bam_dir/$(basename "$bam_file" .bam)_mark_duplicates_metrics.txt"

    picard MarkDuplicates \
        INPUT="$bam_file" \
        OUTPUT="$output_marked_bam" \
        METRICS_FILE="$metrics_file" \
        REMOVE_DUPLICATES=false \
        ASSUME_SORTED=true

    echo "Marked duplicates in $(basename "$bam_file")."
done
# Index the marked BAM files
for marked_bam_file in "$marked_bam_dir"/*.bam; do
    output_index="$marked_bam_dir/$(basename "$marked_bam_file").bai"    

    samtools index "$marked_bam_file" "$output_index"

    echo "Indexed $(basename "$marked_bam_file")."
done
# Directory to save variant calling results
variants_dir="variant_calling_results"
mkdir -p "$variants_dir"

# Run bcftools mpileup and call variants for each indexed BAM file
for marked_bam_file in "$marked_bam_dir"/*.bam; do
    output_vcf="$variants_dir/$(basename "$marked_bam_file" .bam).vcf"

    bcftools mpileup -Ov -f "$reference_file" "$marked_bam_file" | bcftools call -Ov -mv > "$output_vcf"

    echo "Variant calling for $(basename "$indexed_bam_file") completed."
done


