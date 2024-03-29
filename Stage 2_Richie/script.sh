# Stage_2_script #

  mkdir  stage_2
  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Downloaded the forward, reverse and reference sequences
  wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
  wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
  wget https://zenodo.org/records/10886725/files/Reference.fasta
  
  # FASTQC checking and validate QC report
  fastqc *.fastq.gz
  
  # Trimming done using trimmomatic tools to remoe bad quality reads
  conda install trimmomatic
  /usr/bin/TrimmomaticPE
  trimmomatic PE ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz ERR8774458_1_paired.fq.gz ERR8774458_1_unpaired.fq.gz ERR8774458_2_paired.fq.gz ERR8774458_2_unpaired.fq.gz LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:100
  
  # Reference sequence is indexed
  conda install bwa
  bwa index Reference.fasta 
 
  # Mapping the reads to the reference genome
  bwa mem Reference.fasta ERR8774458_1_paired.fq.gz ERR8774458_2_paired.fq.gz >ERR8774458_aligned.sam
  
  # Convert the alignment file from SAM to BAM format
  samtools view -O BAM -o ERR8774458_aligned.bam ERR8774458_aligned.sam
  
  # Sorting the alignment file
  samtools sort -T temp -O bam -o ERR8774458_aligned_sorted.bam ERR8774458_aligned.bam
  
  # Indexing the sorted bam file
  samtools index ERR8774458_aligned_sorted.bam
  
  # Mark the duplicates
  samtools index Reference.fasta
  conda install picard
  picard MarkDuplicates I=ERR8774458_aligned_sorted.bam O=ERR8774458_aligned_markdup.bam M=metrics.txt

  # In order to check how many reads were marked as duplicates we can use the following
   command:
   [grep -A 2 “^##METRICS” metrics.txt]
  
  # Index the above markduplicate file #
  samtools index ERR8774458_aligned_markdup.bam
  
  #Variant calling using bcftools (earlier this was possible using samtools but now it has been removed)
  bcftools mpileup -f Reference.fasta ERR8774458_aligned_markdup.bam >ERR8774458_variants.bcf
  
  # Converting from .bcf file to .vcf file
  bcftools call -mv -O v -o ERR8774458_variants.vcf ERR8774458_variants.bcf
  
  # To see the first 100 lines of the "ERR8774458_variants.vcf" file
  head -n 100 ERR8774458_variants.vcf
  
