# Receive inputs
read -p "Input name of the gene:" name
read -p "Input link to forward read:" R1
read -p "Input link to reverse read:" R2
read -p "Input link to reference fasta data:" reference

# Create directory of name of gene and sub directories
mkdir -p $name/data/ref $name/qc_report $name/trimmed $name/results 

# Download datasets and store in data folder
wget -O "$name"_R1.fastq.gz -P $name/data $R1
wget -O "$name"_R2.fastq.gz -P $name/data $R2
wget -O reference.fasta -P $name/data/ref $reference


# # Perform quality control check on both reads, output to qc_folder folder
# fastqc $name/data/*.fastq.gz -o $name/qc_report

# # Aggregate fastqc reports 
# multiqc $name/qc_report/*_fastqc.zip -o $name/qc_report

# # Trim faulty reads with fastp, output to trimmed folder
# fastp -i $name/data/"$name"*1.fastq.gz -I $name/data/"$name"*2.fastq.gz -o $name/trimmed/"$name"_R1.trimm.fastq.gz -O $name/trimmed/"$name"_R2.trimm.fastq.gz --html $name/trimmed/"$name"_fastp.html --json $name/trimmed/"$name"_fastp.json

# # Map reads to reference genome with bwa
# bwa index $name/data/ref/reference.fasta

# # Align reads to reference genome
# bwa mem $name/data/ref/reference.fasta $name/trimmed/"$name"_R1.trimm.fastq.gz $name/trimmed/"$name"_R2.trimm.fastq.gz > $name/results/"$name".aligned.sam

# # Convert sam file to bam file
# samtools view -S -b $name/results/"$name".aligned.sam > $name/results/"$name".aligned.bam

# # Sort sequemce
# samtools sort -o $name/results/"$name".aligned.sorted.bam $name/results/"$name".aligned.bam

# # Index sorted bam file with samtools
# samtools index $name/results/"$name".aligned.sorted.bam

# # Generate pileup format for bam file
# bcftools mpileup -O b -o $name/results/"$name".bcf -f $name/data/ref/reference.fasta $name/results/"$name".aligned.sorted.bam

# # Identify variants using bcftools call and generately variant (vcf) file
# bcftools call -m -v -o $name/results/"$name".variants.vcf $name/results/"$name".bcf

# # Compress variant file
# # bgzip $name/results/"$name".variants.vcf


