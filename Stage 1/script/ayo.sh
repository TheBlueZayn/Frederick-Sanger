#!/bin/bash

#Creating my working directory
mkdir ayo

#creating anotherdirectory called biocomputing, and changing to it with a single line of command
mkdir biocomputing cd biocomputing

#Downloading 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna && wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk && wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

#moving the .fna file to the folder titled your name directly
mv wildtype.fna ~/ayo/

#removing the duplicate .gbk file
rm wildtype.gpk.1

#confirming whether or not the file is a wildtype or mutant
grep 'tatatata' wildtype.fna

#Printing all the lines that show it is mutant into a new file
 grep ‘tatatata’ wildtype.fna > mutant.txt
#Selecting a favourite gene
echo EBNA2
#Downloading the FASTA format of the gene from NCBI
wget -O PP215920.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=PP215920.1&rettype=fasta&retmode=text"

# Counting the number of lines in the FASTA file (excluding the header)
grep -v '^>' PP215920.1.fasta | sed '/^\s*$/d' | wc -l

# Counting the occurrences of each nucleotide base in the FASTA file
grep -o 'A' PP215920.1.fasta | wc -l
grep -o 'G' PP215920.1.fasta | wc -l
grep -o 'C' PP215920.1.fasta | wc -l
grep -o 'T' PP215920.1.fasta | wc -l

# Calculating the %GC content of the gene
GC_content=$(awk '/^[^>]/ {total += length ($0); gc += gsub(/[GC]/,"",$0)} END {printf "%.2f\n", (gc / total) * 100}' PP215920.1.fasta) && echo "GC content: $GC_content%"

# Creating a nucleotide file titled with my name
vim ayo.fasta

# Echoing a sequence into the file
echo "Numbers of A: $(grep -o 'A' PP215920.1.fasta | wc -l)" >> ayo.fasta && echo "Numbers of G: $(grep -o 'G' PP215920.1.fasta | wc -l)" >> ayo.fasta && echo "Numbers of C: $(grep -o 'C' PP215920.1.fasta | wc -l)" >> ayo.fasta && echo "Numbers of T: $(grep -o 'T' PP215920.1.fasta | wc -l)" >> ayo.fasta

#Saving all codes used in this project to a file in my name
vim ayo.sh

#clear terminal
clear

#listing the files in two folders 
ls
