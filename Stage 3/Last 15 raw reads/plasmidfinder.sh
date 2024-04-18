#!/bin/bash

# Loop through all folders in the current directory
for folder in */; do
    # Extract sample name (adjust based on your naming convention)
    sample_name="${folder%/}"

  abricate --db plasmidfinder ${sample_name}/contigs.fasta --csv  > plasmidfinder_output/${sample_name}.csv
  
done

rm plasmidfinder_output/plasmidfinder_output.csv
rm plasmidfinder_output/SeqSero2results.csv
