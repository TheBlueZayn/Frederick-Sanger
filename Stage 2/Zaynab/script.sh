#!bin/bash

mkdir variant_calling && cd variant_calling

mkdir data
mkdir data/fastq
mkdir data/ref

wget -P data/fastq {$} {$}

wget -P data/ref {$}

mkdir qc_report
fastqc data/fastq/*.fastq.gz -o qc_report

multiqc qc_report/*_fastqc.zip -o qc_report

fastp 



