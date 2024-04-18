#!/bin/bash

conda install fastqc
conda install fastp
conda install bioconda::sra-tools
conda install bioconda::samtools
conda install bioconda::bcftools
conda install bioconda::spades
conda install bioconda::seqsero2
conda install bioconda::plasmidfinder
download-db.sh
conda install bioconda::ncbi-amrfinderplus
amrfinder -u
conda install -c bioconda raxml-ng
