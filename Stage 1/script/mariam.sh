#!/bin/bash
mkdir mariam
mkdir biocomputing && cd "$_"
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
mv wildtype.fna /home/mariam/mariam
rm wildtype.gbk.1
grep 'tatata' wildtype.fna
grep 'tatata' wildtype.fna > mutant.fna
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${HOME}/edirect:${PATH}
esearch -db nucleotide -query "NG_008385.2" | efetch -format fasta > CYP2C9.fna
wc -l < CYP2C9.fna | awk '{print($1-1)}'
vim CYP2C9_noheader.fna 
sed 's/>NG_008385.2 Homo sapiens cytochrome P450 family 2 subfamily C member 9 (CYP2C9), RefSeqGene (LRG_1195) on chromosome 10/ /g' CYP2C9.fna > CYP2C9_noheader.fna
grep -o "A" CYP2C9_noheader.fna | wc -l
grep -o "G" CYP2C9_noheader.fna | wc -l
grep -o "C" CYP2C9_noheader.fna | wc -l
grep -o "T" CYP2C9_noheader.fna | wc -l
vim GCcontent.sh
chmod u+x GCcontent.sh
G=$(grep -o "G" CYP2C9_noheader.fna | wc -l)

echo "G count is: $G"

C=$(grep -o "C" CYP2C9_noheader.fna | wc -l)

echo "C count is: $C"

A=$(grep -o "A" CYP2C9_noheader.fna | wc -l)

echo "A count is: $A"

T=$(grep -o "T" CYP2C9_noheader.fna | wc -l)

echo "T count is: $T"

GC=$(($G + $C))

echo "summation of GC is: $GC"

total=$(( $G + $C + $A + $T ))

echo "total is: $total"
GCcontent= echo "scale=4; $GC / $total" | bc
vim mariam.fasta
A=$(grep -o "A" CYP2C9_noheader.fna | wc -l)

C=$(grep -o "C" CYP2C9_noheader.fna | wc -l)

T=$(grep -o "T" CYP2C9_noheader.fna | wc -l)

G=$(grep -o "G" CYP2C9_noheader.fna | wc -l)

echo "C count is:$C" >> mariam.fasta && echo "A count is:$A" >> mariam.fasta && echo "G count is:$G" >> mariam.fasta && echo "T count is:$T" >> mariam.fasta

vim mariam.fasta
history

