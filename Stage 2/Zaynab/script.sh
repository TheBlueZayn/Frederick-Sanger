#!bin/bash

# Create directory for variant calling
mkdir variant_calling && cd variant_calling

# Create directory for data
mkdir -p data/fastq data/ref

# Download reads and reference genome
wget -P data/fastq https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget -P data/ref https://zenodo.org/records/10886725/files/Reference.fasta

# Perform quality check on reads
mkdir qc_report
fastqc data/fastq/*.fastq.gz -o qc_report

# Aggregate fastqc reports
multiqc qc_report/*_fastqc.zip -o qc_report

# Trim faulty reads with fastp
mkdir trimmed 
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --html trimmed/ERR8774458.fastp.html --json --thread 2 -i data/fastq/ERR8774458_1.fastq.gz -I data/fastq/ERR8774458_2.fastq.gz -o trimmed/ERR8774458_1.fastq.gz -O trimmed/ERR8774458_2.fastq.gz

# Map reads to reference genome
mkdir -p results/sam results/bam results/bcf results/vcf
bwa index data/ref/Reference.fasta
bwa mem data/ref/Reference.fasta trimmed/ERR8774458_1.fastq.gz trimmed/ERR8774458_2.fastq.gz > results/sam/ERR8774458.aligned.sam


samtools view -S -b results/sam/ERR8774458.aligned.sam > results/bam/ERR8774458.aligned.bam


samtools sort -o results/bam/ERR8774458.aligned.sorted.bam results/bam/ERR8774458.aligned.bam

samtools index results/bam/ERR8774458.aligned.sorted.bam

samtools flagstat results/bam/ERR8774458.aligned.sorted.bam


bcftools mpileup -O b -o results/bcf/ERR8774458_raw.bcf -f data/ref/Reference.fasta results/bam/ERR8774458.aligned.sorted.bam

bcftools call --ploidy 1 -m -v -o results/vcf/ERR8774458.variants.vcf results/bcf/ERR8774458_raw.bcf

# Compress variant file
bgzip results/vcf/ERR8774458.variants.vcf
# index variant file
bcftools index results/vcf/ERR8774458.variants.vcf.gz



# Counting all variants
zgrep -v -c "^#" results/vcf/ERR8774458.variants.vcf.gz > variants.txt
echo "Number of variants: $(zgrep -v -c "^#" results/vcf/ERR8774458.variants.vcf.gz)" > variants.txt 

# View single nucleotides polymorphisms (SNPs) and Insertions/deletions
bcftools view -v snps results/vcf/ERR8774458.variants.vcf.gz|grep -v -c "^#"
bcftools view -v indels results/vcf/ERR8774458.variants.vcf.gz|grep -v -c "^#"

# Filter snps and indels into separate files
bcftools view -v snps results/vcf/ERR8774458.variants.vcf.gz -Oz -o results/vcf/snps.vcf.gz
bcftools view -v indels results/vcf/ERR8774458.variants.vcf.gz -Oz -o results/vcf/indels.vcf.gz




samtools view -b -F 0xc results/bam/ERR8774458.aligned.bam -o results/bam/ERR8774458.aligned.filtered.bam
samtools sort -@ 8 -n results/bam/ERR8774458.aligned.filtered.bam -o results/bam/ERR8774458.aligned.sorted.n.bam
samtools fixmate -m results/bam/ERR8774458.aligned.sorted.n.bam results/bam/ERR8774458.aligned.fixmate.bam

samtools sort -@ 8 results/bam/ERR8774458.aligned.fixmate.bam -o results/bam/ERR8774458.aligned.sorted.p.bam
samtools markdup -r -@ 8 results/bam/ERR8774458.aligned.sorted.p.bam results/bam/ERR8774458.aligned.dedup.bam

freebayes -f data/ref/Reference.fasta -b results/bam/ERR8774458.aligned.dedup.bam --vcf results/vcf/ERR8774458.vcf
bgzip results/vcf/ERR8774458.vcf

bcftools index results/vcf/ERR8774458.vcf.gz