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

fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --html trimmed/ERR8774458.fastp.html --json --thread 2 -i data/fastq/ERR8774458_1.fastq.gz -I data/fastq/ERR8774458_2.fastq.gz -o trimmed/ERR8774458_1.fastq.gz -O trimmed/ERR8774458_2.fastq.gz



