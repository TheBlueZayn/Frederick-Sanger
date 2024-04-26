#!/bin/bash

#Creat reference directory and move to that directory
mkdir reference
cd reference

#Download reference genome
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=GCF_008727535.1&rettype=fasta&retmode=text" -O reference/refgenome.fasta

#Move to home directory Creat rawreads directory and move to that directory
cd ~
mkdir rawreads
cd rawreads

# Define accessions number
accessions=(ERR9516208" "ERR9516260" "ERR9516213" "ERR9516215" "ERR9516190" "ERR9516264" "ERR9516186" "ERR9516341" "ERR9516233" "ERR9516202" "ERR9516219" "ERR9516192" "ERR9516194" "SRR13221299" "ERR9516319" "ERR9516292" "ERR9516266" "SRR11927768" "ERR9516298" "SRR11927771" "ERR9516306" "ERR9516178" "ERR9516309" "ERR9516211" "SRR11927779" "ERR9516265" "SRR11903203" "ERR9516204" "ERR9516269" "ERR9516307" "ERR9516279" "ERR9516242" "ERR9855643" "ERR9516234" "ERR9516335" "ERR9516228" "ERR9516338" "ERR9516226" "ERR9516287" "SRR8294664" "ERR9516349" "SRR8327695" "ERR9516354" "SRR8159839" "ERR9516356" "SRR8159842" "ERR9516189" "SRR7514404" "ERR9516191" "ERR3635344")

for acc in "${accessions[@]}"; do
    fastq-dump --split-files --gzip "$acc"
done

#Move to home directory and creat other directories
cd ~ 
mkdir qcdata trimreads spadesdata serotyping plasmidtyping amrgenes treedata

#Define samples
samples=(ERR9516208" "ERR9516260" "ERR9516213" "ERR9516215" "ERR9516190" "ERR9516264" "ERR9516186" "ERR9516341" "ERR9516233" "ERR9516202" "ERR9516219" "ERR9516192" "ERR9516194" "SRR13221299" "ERR9516319" "ERR9516292" "ERR9516266" "SRR11927768" "ERR9516298" "SRR11927771" "ERR9516306" "ERR9516178" "ERR9516309" "ERR9516211" "SRR11927779" "ERR9516265" "SRR11903203" "ERR9516204" "ERR9516269" "ERR9516307" "ERR9516279" "ERR9516242" "ERR9855643" "ERR9516234" "ERR9516335" "ERR9516228" "ERR9516338" "ERR9516226" "ERR9516287" "SRR8294664" "ERR9516349" "SRR8327695" "ERR9516354" "SRR8159839" "ERR9516356" "SRR8159842" "ERR9516189" "SRR7514404" "ERR9516191" "ERR3635344")
for sample in "${samples[@]}"; do

#Quality control check of rawreads
    fastqc rawreads/${sample}_1.fastq.gz rawreads/${sample}_2.fastq.gz -o ~/qcdata/

#Trim the rawreads
    fastp -i "rawreads/${sample}_1.fastq.gz" -I "rawreads/${sample}_2.fastq.gz" -o "trimreads/${sample}_trim_1.fastq.gz" -O "trimreads/${sample}_trim_2.fastq.gz"

#Creat contigs files
    spades.py -1 "trimreads/${sample}_trim_1.fastq.gz" -2 "trimreads/${sample}_trim_2.fastq.gz" -o "spadesdata/${sample}_spades"

#Run SeqSero2 to serotype
   SeqSero2_package.py -m k -t 4 -i "spadesdata/${sample}_spades/contigs.fasta" -d "serotyping/${sample}.tsv"

#Run abricate with plasmidfinder database
   abricate --db plasmidfinder "spadesdata/${sample}_spades/contigs.fasta" --csv > plasmidtyping/${sample}.csv

#Run ncbi-amrfinderplus with amrfinder database
     amrfinder -n "spadesdata/${sample}_spades/contigs.fasta" --organism Salmonella -o "amrgenes/${sample}.csv"

# Check if the output CSV exists and is not empty
    if [ -s "amrgenes/${sample}.csv" ]; then

# Add a new column with the sample name
    awk -v sample_name="$sample" 'BEGIN {FS=OFS=","} 
    NR==1 {print $0, "Sample Name"} 
    NR>1 {print $0, sample_name}' "amrgenes/${sample}.csv" > "amrgenes/${sample}_mod.csv"

# Replace the original CSV with the modified one
    mv "amrgenes/${sample}_mod.csv" "amrgenes/${sample}.csv"
    else
    echo "Warning: No output for ${sample}, skipping."
    fi

#Run snippy for variant calling and core genome alignment
   snippy --cpus 4 \
           --outdir "snippydata/${sample}" \
           --ref reference/refgenome.fna \
           --R1 "trimreads/${sample}_trim_1.fastq.gz" \
           --R2 "trimreads/${sample}_trim_1.fastq.gz"

done

#Move to snippydata directory Run snippy core to creat single core genome alignment file
cd snippydata
snippy-core --ref ~/reference/refgenome.fasta "ERR9516208" "ERR9516260" "ERR9516213" "ERR9516215" "ERR9516190" "ERR9516264" "ERR9516186" "ERR9516341" "ERR9516233" "ERR9516202" "ERR9516219" "ERR9516192" "ERR9516194" "SRR13221299" "ERR9516319" "ERR9516292" "ERR9516266" "SRR11927768" "ERR9516298" "SRR11927771" "ERR9516306" "ERR9516178" "ERR9516309" "ERR9516211" "SRR11927779" "ERR9516265" "SRR11903203" "ERR9516204" "ERR9516269" "ERR9516307" "ERR9516279" "ERR9516242" "ERR9855643" "ERR9516234" "ERR9516335" "ERR9516228" "ERR9516338" "ERR9516226" "ERR9516287" "SRR8294664" "ERR9516349" "SRR8327695" "ERR9516354" "SRR8159839" "ERR9516356" "SRR8159842" "ERR9516189" "SRR7514404" "ERR9516191" "ERR3635344"

cd ~

#Move to tree directory and Run raxml-ng tool
raxml-ng --msa ~/snippydata/core.full.aln --model GTR+G --tree pars{10} --prefix RAxML_result --threads 4 --search

echo "Task completed successfully."
