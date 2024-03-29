# Stage_2_script #

  608  mkdir  stage_2
  # Downloaded the forward, reverse and reference sequences
  762  wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
  763  ls
  764  wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
  765  ls
  766  wget https://zenodo.org/records/10886725/files/Reference.fasta
  # QC quality checked
  770  fastqc *.fastq.gz
  # Trimming done using trimmomatic tools to remoe bad quality reads
  805  conda install trimmomatic
  806  /usr/bin/TrimmomaticPE
  807  trimmomatic PE ERR8774458_1.fastq.gz ERR8774458_2.fastq.gz ERR8774458_1_paired.fq.gz ERR8774458_1_unpaired.fq.gz ERR8774458_2_paired.fq.gz ERR8774458_2_unpaired.fq.gz LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:100
  811  conda install bwa
  # Reference sequence is indexed
  812  bwa index Reference.fasta 
  813  ls
  # Mapping the reads to the reference genome
  814  bwa mem Reference.fasta ERR8774458_1_paired.fq.gz ERR8774458_2_paired.fq.gz >ERR8774458_aligned.sam
  # Convert the alignment file from SAM to BAM format
  818  samtools view -O BAM -o ERR8774458_aligned.bam ERR8774458_aligned.sam
  # Sorting the alignment file
  820  samtools sort -T temp -O bam -o ERR8774458_aligned_sorted.bam ERR8774458_aligned.bam
  # Indexing the sorted bam file
  822  samtools index ERR8774458_aligned_sorted.bam
  824  samtools index Reference.fasta
  # Mark the duplicates
  # In order to check how many reads were marked as duplicates we can use the following
   command:
   [grep -A 2 “^##METRICS” metrics.txt]
  843  conda install picard
  844  picard MarkDuplicates I=ERR8774458_aligned_sorted.bam O=ERR8774458_aligned_markdup.bam M=metrics.txt
  # Index the above markduplicate file #
  846  samtools index ERR8774458_aligned_markdup.bam
  #Variant calling using bcftools (earlier this was possible using samtools but now it has been removed)
  851  bcftools mpileup -f Reference.fasta ERR8774458_aligned_markdup.bam >ERR8774458_variants.bcf
  # Converting from .bcf file to .vcf file
  853  bcftools call -mv -O v -o ERR8774458_variants.vcf ERR8774458_variants.bcf
  # To see the first 100 lines of the "ERR8774458_variants.vcf" file
  855  head -n 100 ERR8774458_variants.vcf
  
