#!/bin/bash

# Step 0: Prerequisite files
# Download gencode.v38.basic annotation from https://www.gencodegenes.org/human/release_38.html (gencode.v38.basic.annotation.gtf)



# Step 1: Extract CFTR2 transcript data
## 1.1 Filter the original GTF file to include only lines that mention "CFTR".

## 1.2 Extract and reformat the relevant information into a CSV.
./extract_cftr_coordinates.sh



# Step 2: Generate exon-intron coordinates
python3 cftr_exon_intron_coordinate_mapper.py gencode.v38.cftr.clean.csv && rm gencode.v38.cftr.clean.csv