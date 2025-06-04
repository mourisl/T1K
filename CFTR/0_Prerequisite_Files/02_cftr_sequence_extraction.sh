#!/bin/bash

# Step 0: Prerequisite files
# Download GRCh38 DNA Chromosome 7 fasta file from https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/ (Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz)



# Step 1: Preparation of the Standardized Chromosome 7 FASTA File
sed 's/>7 .*/>7/' Homo_sapiens.GRCh38.dna.chromosome.7.fa > fixed_chr7.fa



# Step 2: Extraction of the CFTR-201 genomic region (CFTR.dna.fa)
#   Transcript_id: "ENST00000003084.11"
#   Upstream extension: 200 bps before the first UMI (position 117480025)
#      - Start 117480025 - 200 = 117479825
#   Downstream extension: 200 bps after the last UMI (position 117668665)
#      - End 117668665 + 200 = 117668865
samtools faidx fixed_chr7.fa 7:117479825-117668865 > CFTR-201.dna.fa && rm fixed_chr7.fa fixed_chr7.fa.fai



# Step 3: CFTR DNA Sequence Exportation
python3 CFTR-201_genomic_sequence_exporter.py CFTR-201.dna.fa && rm CFTR-201.dna.fa