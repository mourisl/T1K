#!/bin/bash
# Filename: extract_cftr_coordinates.sh
# Description: Extract and reformat relevant CFTR transcript (CFTR-201) information from the GTF file
#              into a CSV file.
# Usage: ./extract_cftr_coordinates.sh
# Requirements: Ensure that gencode.v38.basic.annotation.gtf is in the same directory, or adjust the path accordingly.


output_file="gencode.v38.cftr.clean.csv"
echo "chr7,type,pos1,pos2,exon_num,trans_name" > "$output_file"

# Use grep to filter lines containing "CFTR" from the basic annotation file,
# Process only exon lines with "CFTR-201".
grep -i "CFTR" gencode.v38.basic.annotation.gtf | \
awk 'BEGIN {
    FS="\t"; OFS=","
}
$3 == "exon" && /transcript_name "CFTR-201"/ {
    exon = ""
    # Extract the exon number if available.
    if (match($0, /exon_number ([0-9]+)/, arr))
        exon = arr[1]
    # Print the CSV line.
    print $1, $3, $4, $5, exon, "CFTR-201"
}' >> "$output_file"



