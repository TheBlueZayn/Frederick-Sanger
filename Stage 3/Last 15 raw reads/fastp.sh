#!/bin/bash

# Loop through all FASTQ files in the current directory
for file in *_1.fastq.gz; do
  # Extract sample name from filename (modify as needed)
  sample_name=${file%%_1.fastq.gz}
  
  # fastp command with sample-specific filenames, removing "processed"
  fastp -i "$file" -I "${sample_name}_2.fastq.gz" \
       -o "~/trim.reads/${sample_name}_1.trim.fastq.gz" \
       -O "~/trim.reads/${sample_name}_2.trim.fastq.gz"
done
