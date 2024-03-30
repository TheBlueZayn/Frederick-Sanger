#!/usr/bin/env bash

#Any directory you are come to home directory directly
cd ~

#Creat one new directory name rawdata and two subdirectories names sequence and reference
mkdir rawdata && mkdir rawdata/sequence rawdata/reference

# Change directory from home  to sequence directory Download the forward and reverse sequence in sequence directory
cd rawdata/sequence  &&  wget wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz rawdata/sequence

# Change directory from sequence  to refrence directory Download the reference sequence in reference directory
cd ~/rawdata/reference && wget https://zenodo.org/records/10886725/files/Reference.fasta

# Change to home directory and creat QCreports directory for the results of quality control
cd ~ && mkdir QCreports

#Run fastqc command for quality control and output results in QCreports directory
fastqc rawdata/sequence/ERR8774458_1.fastq.gz rawdata/sequence/ERR8774458_2.fastq.gz -o QCreports/

#Creat multiQCreport directory and use  MultiQC to summarize the QCreports and output results in multiQCreport directory
mkdir multiQCreport && multiqc  ~/QCreports  -o multiQCreport/

# Creat Trimreports directory trim the forward and reverse sequences by fastp tool and output results with name trimforward and trimreverse in Trimreports directory
mkdir Trimreports && cd Trimreports && fastp -i ~/rawdata/sequence/ERR8774458_1.fastq.gz -I ~/rawdata/sequence/ERR8774458_2.fastq.gz -o trimforward.fastq.gz -O trimreverse.fastq.gz

# Change to home directory index reference sequence by bwa index tool and output results into reference directory
cd ~  && bwa index rawdata/reference/Reference.fasta  rawdata/reference/

# Creat Alignmentreports directory align the trimsequence to reference sequence by bwa mem tool and output results in Aligmentreports directory
mkdir Alignmentreports && bwa mem rawdata/reference/Reference.fasta Trimreports/trimforward.fastq.gz Trimreports/trimreverse.fastq.gz > Alignmentreports/Alignsequence.sam

#Convert the sam file to a bam file by samtools and output the results into Alignmentreports directory
samtools view -S -b Alignmentreports/Alignsequence.sam > Alignmentreports/Alignsequence.bam

#Creat sortedsequence directory  and change directory to sortedsequence Sort bam file using samtoolssort and output results in name sortedAlignsequence
mkdir sortedsequence && cd sortedsequence && samtools sort ~/Alignmentreports/Alignsequence.bam -o SortedAlignsequence.bam

# Index the sorted bam file using samtools index tools
samtools index SortedAlignsequence.bam

# Variant calling by bcftools
cd ~ && mkdir Variant && cd ~/Variant && bcftools mpileup -Ou -f ~/rawdata/reference/Reference.fasta ~/sortedsequence/SortedAlignsequence.bam | bcftools call -mv -Ov -o ~/Variant/variants.vcf

# Convert the vcf file to  csv file for further analysis and visualisation
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' variants.vcf > variant.csv

# Creat snps and indels files from variant vcf file using bcftools
bcftools view -v snps variants.vcf > snps.vcf && bcftools view -v indels variants.vcf  >  indels.vcf

# Convert snps and indel vcf files to csv files using bcftools
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' ~/Variant/snps.vcf  > snps.csv && bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' ~/Variant/indels.vcf  >  indels.csv

# Count the number of snps and indels from csv files
wc -l snps.csv && wc -l indels.csv

