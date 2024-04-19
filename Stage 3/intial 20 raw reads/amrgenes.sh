#!/usr/bin/env bash

# Define an array of samples
samples=("ERR9516208" "ERR9516260" "ERR9516213" "ERR9516215" "ERR9516190" "ERR9516264" "ERR9516186" "ERR9516341" "ERR9516233" "ERR9516202" "ERR9516219" "ERR9516192" "ERR9516194" "SRR13221299" "ERR9516319" "ERR9516292" "ERR9516266" "SRR11927768" "ERR9516298" "SRR11927771")

# Loop over each sample
for sample in "${samples[@]}"; do
    # Run amrfinderplus 

amrfinder -n "spadesdata/${sample}_spades/contigs.fasta" --organism Salmonella -o "amrgenes/${sample}.csv"

    # Check if the output CSV exists and is not empty
    if [ -s "amrgene/${sample}.csv" ]; then
        # Add a new column with the sample name
        awk -v sample_name="$sample" 'BEGIN {FS=OFS=","} 
        NR==1 {print $0, "Sample Name"} 
        NR>1 {print $0, sample_name}' "amrgenes/${sample}.csv" > "amrgenes/${sample}_mod.csv"

        # Replace the original CSV with the modified one
        mv "amrgenes/${sample}_mod.csv" "amrgenes/${sample}.csv"
    else
        echo "Warning: No output for ${sample}, skipping."
    fi
done
