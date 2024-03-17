$ mkdir FAIZAN
$ mkdir biocomputing && cd biocomputing
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
$ wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
$ mv wildtype.fna ../FAIZAN
$ rm wildtype.gbk.1
$ if grep -q "tatatata" wildtype.fna; then echo "mutant"; else echo "wildtype"; fi
grep "tatatata" wildtype.fna > mutant.fna
$ wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NM_001256799 &rettype=fasta&retmode=text" -O gadphgene.fasta
$ grep -c -v "^>" gadphgene.fasta
$ grep -o "A" gadphgene.fasta | wc -l
$ grep -o "G" gadphgene.fasta | wc -l
$ grep -o "C" gadphgene.fasta | wc -l
$ grep -o "T" gadphgene.fasta | wc -l
$ awk '!/^>/{n+=length($0); g+=gsub(/[GgCc]/,"") } END{print g/n * 100 "%"}' gadphgene.fasta
$ nano FAIZAN.fasta
$ C=$(grep -o "C" gadphgene.fasta | wc -l) && A=$(grep -o "A" gadphgene.fasta | wc -l) && G=$(grep -o "G" gadphgene.fasta | wc -l) && T=$(grep -o "T" gadphgene.fasta | wc -l) | echo "G count is:$G" >> FAIZAN.fasta && echo "A count is:$A" >> FAIZAN.fasta && echo "C count is:$C" >> FAIZAN.fasta && echo "T count is:$T" >> FAIZAN.fasta
$ history 
$ clear
$ ls
