#! usr/bin/bash

# Making a folder titled my name
$ mkdir zaynab

# Making a folder named "biocomputing" and changing to it
$ mkdir biocomputing && cd biocomputing

# Downloading sample datasets
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# Moving file to my folder
$ mv wildtype.fna ../zaynab

# Deleting a duplicated file
$ rm wildtype.gbk.1

# Checking if file is a mutant or not
$ cd ~/Hackbio/zaynab && if grep -q "tatatata" wildtype.fna; then echo "mutant"; else echo "wildtype"; fi

# Moving mutant part into a new file
$ grep "tatatata" wildtype.fna > mutant.fna

# changing directory and downloading my gene of choice "CYP2D6" in FASTA format
cd biocomputing && wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=JX678711.1 &rettype=fasta&retmode=text" -O CYP2D6.fasta

# Counting the number of lines in the FASTA file
$ grep -c -v "^>" CYP2D6.fasta

# Counting frequency of "A"
$ grep -o "A" CYP2D6.fasta | wc -l

# Counting frequency of "G"
$ grep -o "G" CYP2D6.fasta | wc -l

# Counting frequency of "C"
$ grep -o "C" CYP2D6.fasta | wc -l

# counting frequency of "T"
$ grep -o "T" CYP2D6.fasta | wc -l

# Calculating the "GC" content
$ awk '!/^>/{n+=length($0); g+=gsub(/[GgCc]/,"") } END{print g/n * 100 "%"}' CYP2D6.fasta

# Creating a FASTA file titled my name
$ nano Zaynab.fasta

# Adding the number of "A", "G", "T", and "C" into the file
echo "Number of A: $(grep -o 'A' CYP2D6.fasta | wc -l)" >> Zaynab.fasta && echo "Number of G: $(grep -o 'G' CYP2D6.fasta | wc -l)" >> Zaynab.fasta && echo "Number of T: $(grep -o 'T' CYP2D6.fasta | wc -l)" >> Zaynab.fasta && echo "Number of C: $(grep -o 'C' CYP2D6.fasta | wc -l)" >> Zaynab.fasta
$ clear
$ history

