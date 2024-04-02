# Create a new directory and change to it
mkdir variant_calling && cd variant_calling 

# Create a directory for fastq data
mkdir data/fastq

# Download fastq data
wget -P data/fastq https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz

# Create a directory for reference data
mkdir data/ref

# Download reference data
wget -P data/ref https://zenodo.org/records/10886725/files/Reference.fasta

# Conduct quality control on the fastq data
fastqc data/fastq/*.fastq.gz

# Aggregate the fastq data reports
multiqc data/fastq/*_fastqc.zip

# Create a 'trimmed' directory
mkdir trimmed

# Trim fastq data and output in the trimmed folder
fastp -i data/fastq/ERR8774458_1.fastq.gz -I data/fastq/ERR8774458_2.fastq.gz -o trimmed/trim_1.fastq.gz -O trimmed/trim_2.fastq.gz

# 
