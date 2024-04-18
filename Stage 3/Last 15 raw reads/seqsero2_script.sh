#!/bin/bash

for folder in */; do
  # Remove trailing slash from folder name
  sample_name="${folder%/}"

  # Print the extracted sample name
  SeqSero2_package.py -m k -t 4 -i ${sample_name}/contigs.fasta -n ${sample_name}
done
