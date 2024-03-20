1. $ mkdir yousra_wali
2. $ mkdir  ir biocomputing && cd biocomputing
3. $ s://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk	https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
4. $ mv wildtype.fna /home/usra_akhss/yousra_wali
5. $ rm wildtype_copy.gbk
6. $ grep ‘tatatata’ wildtype.fna
7. $ grep ‘tatatata’ wildtype.fna > mutant.fna
8. $ wget" httpss://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NM_005228.5 &rettype=fasta&retmode=text" -O gene.fasta
9. $ wc-l < gene.fasta | awk '{print($1-1)}'
10. $ grep -v '^>' gene.fasta | tr -d '\n' | grep -o 'A' | wc –l 
11. $ grep -v '^>' gene.fasta | tr -d '\n' | grep -o 'G' | wc –l 
12. $ grep -v '^>' gene.fasta | tr -d '\n' | grep -o 'C' | wc –l 
13. $ grep -v '^>' gene.fasta | tr -d '\n' | grep -o 'T' | wc –l
14. $ awk '/^>/ {next} {total+=length($0); gc+=gsub(/[GCgc]/,"")} END {print (gc / total) * 100}' gene.fasta
15. $ vim yousra1.txt
16. $ echo G=2282, A=2736, C=2451, T=2436 >> yousra1.txt
17. $ history 
18. $ clear
19. $ ls
