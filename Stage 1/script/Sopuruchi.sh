# Creating a folder titled my name
$ mkdir Sopuruchi_Mba
# Creating a new folder titled biocomputing and changing to it
$ mkdir biocomputing; cd biocomputing
# Downloading 3 files
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
# Moving .fna file to my folder
$ mv wildtype.fna /home/sopuruchilizbeth/Sopuruchi_Mba
# Deleting the duplicate file
$ rm wildtype.gbk1
# Confirming if the file is mutant or wild type 
$ if grep -q "tatatata" wildtype.fna; then echo "mutant"; else echo "wildtype"; fi
# Printing the lines that show it is a mutant into a new file
$ grep "tatatata" wildtype.fna > mutant.fna
# Downloading the fasta format of my gene of choice (PDCD1)
$ wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NM_005018.3 &rettype=fasta&retmode=text" -O pdcd1gene.fasta 
#Counting the number of lines are in the FASTA file
$ grep -c -v "^>" pdcd1gene.fasta
# Counting the frequency of A
$ grep -v ">" pdcd1gene.fasta | tr -d '\n' | awk -F "A" '{count+=NF-1} END {print count}'
# Counting the frequency of G
$ grep -v ">" pdcd1gene.fasta | tr -d '\n' | awk -F "G" '{count+=NF-1} END {print count}'
# Counting the frequency of C
$ grep -v ">" pdcd1gene.fasta | tr -d '\n' | awk -F "C" '{count+=NF-1} END {print count}'
# Counting the frequency of T
$ grep -v ">" pdcd1gene.fasta | tr -d '\n' | awk -F "T" '{count+=NF-1} END {print count}'
# Calculating the %GC content of my gene
$ awk '!/^>/{n+=length($0); g+=gsub(/[GgCc]/,"") } END{print g/n * 100 "%"}' pdcd1gene.fasta
# Creating a nucleotide (.fasta) file titled my name
$ touch sopuruchi.fasta
# Adding the number of A, G, T, and C into the nucleotide file
$ echo 'A count is: 382' >> sopuruchi.fasta && echo 'G count is: 651' >> sopuruchi.fasta && echo 'T count is: 351' >> sopuruchi.fasta && echo 'C count is: 713' >> sopuruchi.fasta
$ clear
$ history 
