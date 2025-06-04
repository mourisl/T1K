#!/usr/bin/env python3
import sys
import pandas as pd
import re


def read_fasta(file_path):

    """
    Read a FASTA file and return the sequence as a single string.
    Assumes the file has a header line (starting with '>') that is skipped.
    """

    with open(file_path, 'r') as file:
        # Skip the header line (starts with '>')
        header = file.readline()
        # Read the sequence (remove newlines and join into a single string)
        sequence = ''.join(line.strip() for line in file)
    return sequence

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(1)

    # Get the FASTA file path from the command line argument
    fasta_file = sys.argv[1]

# Read the sequence
dna_sequence = read_fasta(fasta_file)

# Create the content for the Python script
script_content = f"""\
# DNA sequence for CFTR gene
CFTR_DNA = ("{dna_sequence}")
"""

# Path to the output Python script
output_script = 'CFTR_Sequence.py'

# Write the content to the Python script
with open(output_script, 'w') as file:
    file.write(script_content)
