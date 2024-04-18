#!/bin/bash

# Exit immediately if any command returns a non-zero status
set -e
# Exit if any command in a pipeline fails
set -o pipefail

# Make necessary directories
mkdir -p data/{raw_reads,trimmed_reads,reference,raw_reads/fastqc_raw}
mkdir -p results/{amr_detection,read_assembly,seq_sero,abricate}

#Define the URL for the reference genome and the desired filename
reference_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NZ_CP044177.1&rettype=fasta&retmode=text"
reference_name="data/reference/reference.fasta"

# Download the reference genome
if wget -O "$reference_name" "$reference_url"; then
echo "Downloaded reference genome."
else
    echo "Failed to download reference genome."
    exit 1
fi

# Define an array of sample URLs
sample_urls=(
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/006/ERR9516306/ERR9516306_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/006/ERR9516306/ERR9516306_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/008/ERR9516178/ERR9516178_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/008/ERR9516178/ERR9516178_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/009/ERR9516309/ERR9516309_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/009/ERR9516309/ERR9516309_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/001/ERR9516211/ERR9516211_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/001/ERR9516211/ERR9516211_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/079/SRR11927779/SRR11927779_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/079/SRR11927779/SRR11927779_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/005/ERR9516265/ERR9516265_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/005/ERR9516265/ERR9516265_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/003/SRR11903203/SRR11903203_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/003/SRR11903203/SRR11903203_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/004/ERR9516204/ERR9516204_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/004/ERR9516204/ERR9516204_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/009/ERR9516269/ERR9516269_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/009/ERR9516269/ERR9516269_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/007/ERR9516307/ERR9516307_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/007/ERR9516307/ERR9516307_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/009/ERR9516279/ERR9516279_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/009/ERR9516279/ERR9516279_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/002/ERR9516242/ERR9516242_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/002/ERR9516242/ERR9516242_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR985/003/ERR9855643/ERR9855643_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR985/003/ERR9855643/ERR9855643_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/004/ERR9516234/ERR9516234_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/004/ERR9516234/ERR9516234_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/005/ERR9516335/ERR9516335_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR951/005/ERR9516335/ERR9516335_2.fastq.gz"

)

# Download samples
for url in "${sample_urls[@]}"; do
   if curl -L "$url" -o "data/raw_reads/$(basename $url)"; then
  echo "Downloaded raw reads."
    else
        echo "Failed to download raw reads."
        exit 1
    fi
done


# Main analysis

# Quality control using Fastqc
for file in data/raw_reads/*.fastq.gz; do
  if fastqc "$file" -o data/raw_reads/fastqc_raw/; then
  echo "Finished FastQC for raw reads."
  else
      echo "FastQC failed for raw reads."
      exit 1
  fi
done

# Trimming raw reads using Fastp
for file1 in data/raw_reads/*_1.fastq.gz; do
    file2="${file1/_1.fastq.gz/_2.fastq.gz}"
    
    # Check if the corresponding file for paired-end reads exists
    if [ ! -f "$file2" ]; then
        echo "Missing second file for paired-end reads: $(basename $file1)"
        exit 1
    fi

    sample=$(basename "$file1" _1.fastq.gz)
    
    if fastp -i "$file1" \
             -I "$file2" \
             -o "data/trimmed_reads/${sample}_trim_1.fastq.gz" \
             -O "data/trimmed_reads/${sample}_trim_2.fastq.gz"; then
        echo "Finished trimming: ${sample}"
    else
        echo "Trimming failed: ${sample}"
        exit 1
    fi
done

# Assembly using SPAdes
for file1 in data/trimmed_reads/*_trim_1.fastq.gz; do
    file2="${file1/_trim_1.fastq.gz/_trim_2.fastq.gz}"

    # Check if the corresponding file for paired-end reads exists
    if [ ! -f "$file2" ]; then
        echo "Missing second file for paired-end reads: $(basename $file1)"
        exit 1
    fi

    sample=$(basename "$file1" _trim_1.fastq.gz)

    if spades.py --isolate -1 "$file1" \
                 -2 "$file2" \
                 -o "results/read_assembly/${sample}_assembly"; then
        echo "Finished assembly using SPAdes: ${sample}"
    else
        echo "Assembly using SPAdes failed: ${sample}"
        exit 1
    fi
done


# AMR detection using Abricate
# Loop through all folders in the results/read_assembly directory
for assembly in results/read_assembly/*; do
    # Extract sample name (assuming the folder name is the sample name)
    sample_name=$(basename "$assembly" _assembly)

    # Check if contigs.fasta exists in the folder
    if [ -f "$assembly/scaffolds.fasta" ]; then
        abricate --db plasmidfinder "$assembly/scaffolds.fasta" --csv > "results/abricate/${sample_name}.csv"
        echo "Processed ${sample_name}"
    else
        echo "scaffolds.fasta not found in ${sample_name}"
    fi
done


# AMR prediction using AMRFinder
# Loop through all folders in the results/read_assembly directory
for assembly in results/read_assembly/*; do
    # Extract sample name (assuming the folder name is the sample name)
    sample=$(basename "$assembly" _assembly)

    # Check if contigs.fasta exists in the folder
    if [ -f "$assembly/contigs.fasta" ]; then
        # Run amrfinder
        if amrfinder -n "$assembly/contigs.fasta" -O Salmonella -o "results/amr_detection/${sample}_amrfinder.csv"; then
            echo "Finished AMR prediction using AMRFinder for sample: ${sample}"
        else
            echo "AMR prediction using AMRFinder failed for sample: ${sample}"
            exit 1
        fi
    else
        echo "contigs.fasta not found in ${sample}"
    fi
done

#!/bin/bash

# Define an array of sample paths
samples=(results/read_assembly/*_assembly/contigs.fasta)

# Define the output directory for SeqSero2 results
output_dir="results/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over each sample
for sample in "${samples[@]}"; do
    # Get the sample name without the path and extension
    sample_name=$(basename "$(dirname "$sample")")

    # Define the output file path
    output_file="${output_dir}${sample_name}_contig.fasta"

    # Copy the scaffold.fasta file to the new directory and rename it
    cp "$sample" "$output_file"
done

# Put them all in a new directory under results
mkdir results/contig_files && mv results/*.fasta results/contig_files


# Define an array of sample paths
samples=(results/read_assembly/*_assembly/contigs.fasta)

# Define the output directory for copied contig files
output_dir="results/contig_files/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over each sample to copy contig files
for sample in "${samples[@]}"; do
    # Get the sample name without the path and extension
    sample_name=$(basename "$(dirname "$sample")")

    # Define the output file path
    output_file="${output_dir}${sample_name}_contig.fasta"

    # Copy the contig.fasta file to the new directory and rename it
    cp "$sample" "$output_file"
done

# Define an array of contig samples
contig_samples=(results/contig_files/*_contig.fasta)

# Loop over each contig sample to run SeqSero2
for contig_sample in "${contig_samples[@]}"; do
    # Get the sample name without the extension
    sample_name=$(basename "$contig_sample" _contig.fasta)

    # Run SeqSero2 using the contig file
    SeqSero2_package.py -m k -t 4 -i "$contig_sample" -d "seq_sero/${sample_name}"
done

