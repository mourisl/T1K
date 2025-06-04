#!/bin/bash
#==========================================================================
# CFTR2 Transcript and Genomic Sequence Pipeline
#==========================================================================
#
# STEP 0: PREREQUISITE FILES
#
#   1. Download the GTF annotation file from:
#         https://www.gencodegenes.org/human/release_38.html
#      (Filename: gencode.v38.basic.annotation.gtf)
#
#   2. Download the GRCh38 DNA Chromosome 7 FASTA file from:
#         https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/
#      (Filename: Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz)
#
#   3. Package Dependencies:
#         - samtools
#
#==========================================================================
#
# STEP 1: PIPELINE EXECUTION
#
# To execute the complete pipeline, run the following commands sequentially:
#
#   {PATH}/01_parse_cftr_gtf.sh
#   {PATH}/02_cftr_sequence_extraction.sh
#
# Replace {PATH} with the appropriate directory path where the scripts reside.
#
#==========================================================================
#
# EXPLANATION OF 01_parse_cftr_gtf.sh 
#
# PURPOSE: CFTR2 TRANSCRIPT DATA EXTRACTION AND COORDINATE GENERATION
#
#   This step performs the following:
#     1.1. Filters the GTF file to include only lines that mention "CFTR".
#     1.2. Extracts and reformats the relevant information into a CSV file.
#     1.3. Generates exon-intron coordinates.
#
#==========================================================================
#
# EXPLANATION OF 02_cftr_sequence_extraction.sh
#
# PURPOSE: CFTR-201 GENOMIC SEQUENCE PREPARATION AND EXPORTATION
#
#   This step involves:
#     2.1. Standardizing the Chromosome 7 FASTA header.
#     2.2. Extracting the CFTR-201 genomic region using samtools.
#     2.3. Exporting the CFTR-201 DNA sequence as a Python variable for downstream analysis.
#
#==========================================================================


