#!/bin/bash

# Loop through all folders in the current directory
for folder in */; do
    # Extract sample name (adjust based on your naming convention)
    sample_name="${folder%/}"

  amrfinder -n ${sample_name}/contigs.fasta -O Salmonella --name ${sample_name} -o AMRoutput/${sample_name}.csv
done

rm AMRoutput/plasmidfinder_output.csv
rm AMRoutput/SeqSero2results.csv
rm AMRoutput/AMRoutput.csv
