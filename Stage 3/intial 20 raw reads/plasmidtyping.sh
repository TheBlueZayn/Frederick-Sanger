#!/usr/bin/env bash

# Define an array of samples
samples=("ERR9516208" "ERR9516260" "ERR9516213" "ERR9516215" "ERR9516190" "ERR9516264" "ERR9516186" "ERR9516341" "ERR9516233" "ERR9516202" "ERR9516219" "ERR9516192" "ERR9516194" "SRR13221299" "ERR9516319" "ERR9516292" "ERR9516266" "SRR11927768" "ERR9516298" "SRR11927771")

# Loop over each sample
for sample in "${samples[@]}"; do
    # Run abricate using plasmidfinder database
    abricate --db plasmidfinder "spadesdata/${sample}_spades/contigs.fasta" --csv > plasmidtyping/${sample}.csv
done
