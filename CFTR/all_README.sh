#!/bin/bash
#==========================================================================
# CFTR2 Comprehensive Variant Genotyping and Analysis Pipeline
#==========================================================================
#
# OVERVIEW:
#
# This pipeline integrates several modules to process and genotype CFTR2 variants.
# It encompasses the following major steps:
#
#   STEP 0: CFTR2 Transcript and Genomic Sequence Pipeline
#     - Preprocesses GTF and FASTA files to generate reference allele sequences.
#     - Produces key reference files:
#           • CFTR_Sequence.py
#           • CFTR-201_Exon_Intron_Complete_Coordinates.csv
#
#   STEP 1a: CFTR2 Variant Integration and Ensembl Formatter Pipeline
#     - Integrates CFTR2 variant data and processes variant cDNA names to create
#       mutant Ensembl-formatted coordinates.
#     - Command-line arguments allow customization (e.g., allele frequency threshold).
#     - Produces output files:
#           • CFTR_Mimic_Ensembl_Format.dat
#           • CFTR_cDNA_Legacy_Allele_Reference.csv
#         Optionally for debugging:
#           • Combined_Variants_to_Drop.csv
#           • Combined_Variants_to_Keep.csv
#
#   STEP 1b: CFTR2 T1K Reference Sequence Extraction and Analysis
#     - Executes the T1K pipeline to genotype CFTR2 using the generated reference
#       allele sequences.
#     - Reference sequences (CFTR2_dna.fa and/or CFTR2_rna.fa) are generated via
#       the Perl script ParseDatFile.pl.
#
#   STEP 2: RUN T1K Analysis
#     - Performs variant calling using the T1K pipeline on sequencing data.
#
# PREREQUISITE FILES AND MODULES:
#
#   1. CFTR2 Variant List:
#         • CFTR2_25September2024.xlsx (download from http://cftr2.org/mutations_history)
#
#   2. Variant Input File for Downstream Analysis:
#         • Variant_Name_Input.xlsx (prepared from CFTR2_25September2024.xlsx)
#
#   3. Files from the Transcript and Genomic Sequence Pipeline:
#         • CFTR_Sequence.py
#         • CFTR-201_Exon_Intron_Complete_Coordinates.csv
#
#   4. Required Python Modules:
#         • Genomic_Coordinate_Mapping.py
#         • Germline_Ensembl_Variant_Formatter.py
#         • CFTR_Sequence.py
#         • VariantMappingAndMutantEnsemblFormatUtils.py
#         • Codon_AA.py
#
#   5. Required Modules:
#         • T1K: Clone and install from https://github.com/mourisl/T1K.git
#
#   6. Package Dependencies:
#         • samtools
#         • Python (with pandas, re, csv, argparse)
#
# EXECUTION:
#
# To run the entire pipeline, execute the steps sequentially:
#
#   STEP 0: Transcript and Genomic Sequence Pipeline
#     - Run the following scripts to generate the reference allele sequences:
#           {PATH}/01_parse_cftr_gtf.sh
#           {PATH}/02_cftr_sequence_extraction.sh
#
#   STEP 1a: Variant Integration and Ensembl Formatter Pipeline
#     - Execute:
#           python3 Variant_Integration_Ensembl_Formatting.py
#       (Optional command-line arguments may be provided as needed.)
#
#   STEP 1b: CFTR2 T1K Reference Sequence Extraction and Analysis Pipeline
#     - Generate reference FASTA files from the Ensembl-formatted .dat file:
#
#       For DNA:
#           perl {PATH}/ParseDatFile.pl CFTR_Mimic_Ensembl_Format.dat --mode dna > CFTR2_dna.fa
#
#       For RNA:
#           perl {PATH}/ParseDatFile.pl CFTR_Mimic_Ensembl_Format.dat --mode rna > CFTR2_rna.fa
#
#       (Replace {PATH} with the directory path where ParseDatFile.pl is located.)
#
#   STEP 2: T1K Analysis
#     - Run the T1K pipeline to genotype CFTR2 using paired-end FASTQ files:
#
#       For RNA-seq data:
#           ./run-t1k -1 read_1.fq -2 read_2.fq -f CFTR2_rna.fa
#
#       For whole genome or exome sequencing data:
#           ./run-t1k -1 read_1.fq -2 read_2.fq -f CFTR2_dna.fa
#
# Each module generates output files that serve as input for the subsequent steps.
#
#==========================================================================