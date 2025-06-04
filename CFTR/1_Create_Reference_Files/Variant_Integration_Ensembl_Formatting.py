import pandas as pd
import re
import csv
import sys
import argparse

from Genomic_Coordinate_Mapping import generate_result_mapping
from Germline_Ensembl_Variant_Formatter import create_original_Ensembl_format
from CFTR_Sequence import CFTR_DNA as CFTR2
from Codon_AA import translate_codon
from VariantMappingAndMutantEnsemblFormatUtils import (
    map_cdna_to_dna,
    add_200_to_numbers,
    clean_cdna_name,
    clean_number,
    extract_numbers_with_logic,
    check_overlap_mutation_regions,
    determine_sequence,
    build_sequence,
    get_region_for_position,
    map_regions_for_dna_pos,
    final_bp_counts,
    clean_bp_with_indicators,
    assign_region_names,
    create_mutant_Ensembl_format, 
    build_cdna_sequence, 
    translate_full_sequence,
    assign_protein_family_and_allele,
    export_to_dat,
    format_dna_sequence
)



# -----------------------------------------------------------------------------
# Command-Line Argument Parsing
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Variant Integration and Ensembl Formatting Pipeline for CFTR2."
)
parser.add_argument(
    "-allelethreshold", 
    type=float, 
    default=0.01,
    help="Allele frequency threshold for creating combined alleles (default: 0.01)"
)
parser.add_argument(
    "-completepos_dir",
    type=str,
    default="./CFTR-201_Exon_Intron_Complete_Coordinates.csv",
    help="Path to the CFTR exon/intron coordinate CSV file (default: ./CFTR-201_Exon_Intron_Complete_Coordinates.csv)"
)
parser.add_argument(
    "-variantnames_dir",
    type=str,
    default="./Variant_Names_Input.xlsx",
    help="Path to the variant name list (default: ./Variant_Names_Input.xlsx)"
)
parser.add_argument(
    "-export_kept_and_dropped",
    action="store_true",
    help="Export the dropped variants DataFrame and kept variants DataFrame to CSV (default: do not export)"
)
args = parser.parse_args()

allele_freq_threshold = args.allelethreshold
complete_pos = args.completepos_dir
mutation_pos = args.variantnames_dir
export_kept_and_dropped = args.export_kept_and_dropped



# -------------------------------
# Step 1: Read Variant Name List 
# -------------------------------

# Read the variant name list
df = pd.read_excel(mutation_pos)

# Parse the "or" condition: split the "Variant cDNA name" column by "|"
df['Variant cDNA name'] = df['Variant cDNA name'].str.split('|')
df = df.explode('Variant cDNA name')
df.reset_index(drop=True, inplace=True)



# -------------------------------
# Step 2: Expand Variant Name List by Creating Combined Alleles 
# -------------------------------

# Determine top variants based on the allele frequency threshold provided from the command line.
# Build lookup from each original cDNA name
class_map = df.set_index('Variant cDNA name')['Class'].to_dict()

# Select top variants
top_variants = df.loc[
    df['Allele frequency'] >= allele_freq_threshold,
    ['Variant cDNA name', 'Variant legacy name']
].reset_index(drop=True)

# repare a list to hold all rows (original + paired)
records = []

# Add all the original alleles (preserves their original Class)
for _, row in df.iterrows():
    records.append({
        'Variant cDNA name': row['Variant cDNA name'],
        'Variant legacy name': row['Variant legacy name'],
        'Class': row['Class']
    })

# Generate every combination of a “top” variant with every *other* variant
for i in range(len(top_variants)):
    v_i   = top_variants.loc[i, 'Variant cDNA name']
    l_i   = top_variants.loc[i, 'Variant legacy name']
    d_i   = v_i.replace('c.', '')
    
    for j in range(i+1, len(df)):
        v_j = df.loc[j, 'Variant cDNA name']
        l_j = df.loc[j, 'Variant legacy name']
        d_j = v_j.replace('c.', '')
        
        # build the combined cDNA name
        if '[' in d_i:
            d_j_clean     = d_j.strip('[]')
            combined_cDNA = f'c.[{d_i.strip("[]")};{d_j_clean}]'
        elif '[' in d_j:
            d_i_clean     = d_i.strip('[]')
            combined_cDNA = f'c.[{d_j.strip("[]")};{d_i_clean}]'
        else:
            combined_cDNA = f'c.[{d_i};{d_j}]'
        
        # build the combined legacy name
        combined_legacy = f'{l_i};{l_j}'
        
        # look up the two component classes
        c1 = class_map.get(v_i, "")
        c2 = class_map.get(v_j, "")
        
        # apply your hierarchy
        if 'CF-causing' in (c1, c2):
            combined_class = 'CF-causing'
        elif 'Varying clinical consequence' in (c1, c2):
            combined_class = 'Varying clinical consequence'
        else:
            combined_class = 'Non CF-causing'
        
        # append the new record
        records.append({
            'Variant cDNA name':    combined_cDNA,
            'Variant legacy name':  combined_legacy,
            'Class':                combined_class
        })

# Build final DataFrame
df_combined = pd.DataFrame(records)
df_combined = df_combined.drop(columns=[col for col in df_combined.columns
                                        if col == 'Allele frequency'])



# -------------------------------
# Step 3: Establish Coordinate Ranges for Genomic DNA and cDNA Positions
# -------------------------------

# Read the complete CFTR exon/intron coordinates
df_pos = pd.read_csv(complete_pos)

# Generate the mapping between DNA (0-indexed) and RNA positions.
result_mapping = generate_result_mapping(df_pos)



# -------------------------------
# Step 4: Exclude Duplicate Variants and Combined Variants with Overlapping Mutation Regions
# -------------------------------

variants_to_drop_final = []
variants_to_keep_final = []
adjusted_variant_names = {}

for index, row in df_combined.iterrows():
    variant_name = row['Variant cDNA name']

    # Clean the variant name and extract mutation details.
    cleaned_variant = clean_cdna_name(variant_name)
    pos, mutation, seq_change, dna_pos = extract_numbers_with_logic(cleaned_variant, result_mapping)

    # Check for overlaps among mutation regions and adjust the variant name accordingly
    variants_to_drop_final, variants_to_keep_final, adjusted_variant_names = check_overlap_mutation_regions(
        [dna_pos], [mutation], [variant_name],
        variants_to_drop=variants_to_drop_final,
        variants_to_keep=variants_to_keep_final,
        adjusted_variant_names=adjusted_variant_names
    )


# -----------------------------------------------------------------------------
# Step 4A: Create DataFrame of Variants for Downstream Pipeline
# -----------------------------------------------------------------------------

# Create a DataFrame for variants without overlapping mutation regions.
df_keep = pd.DataFrame({'Variant cDNA name': variants_to_keep_final})
df_keep = df_keep.merge(df_combined[['Variant cDNA name', 'Variant legacy name', 'Class']], on='Variant cDNA name', how='left')

# Map the original variant names to their adjusted names.
df_keep['Adjusted Variant cDNA name'] = df_keep['Variant cDNA name'].map(adjusted_variant_names)
df_keep = df_keep.drop(columns=['Variant cDNA name'])

# Remove any duplicate entries in the adjusted variant names, keeping the first occurrence.
df_keep = df_keep.drop_duplicates(subset=['Adjusted Variant cDNA name'])


# -----------------------------------------------------------------------------
# Step 4B: Save DataFrame for Variants Included and Excluded
# -----------------------------------------------------------------------------

if export_kept_and_dropped:
    df_drop = pd.DataFrame({'Variant cDNA name': variants_to_drop_final})
    df_drop = df_drop.merge(df_combined[['Variant cDNA name', 'Variant legacy name', 'Class']],
                            on='Variant cDNA name', how='left')
    df_drop.to_csv("Combined_Variants_to_Drop.csv", index=False)
    df_keep.to_csv("Combined_Variants_to_Keep.csv", index=False)



# -------------------------------
# Step 5: Create Germline Ensembl Format
# -------------------------------

# Extract the num_bp values as a list from the coordinate DataFrame.
num_bp_values = df_pos['bp_num'].tolist()

# Prepend 200 before the first item (5' UTR) and append 200 after the last item (3' UTR).
num_of_bp_origin = [200] + num_bp_values + [200]

# Create the germline Ensembl format using the computed base pair counts
original_Ensembl = create_original_Ensembl_format(num_of_bp_origin)



# -------------------------------
# Step 6: Generate Mutant Ensembl Format for All Variants
# -------------------------------

# Initialize the output list for mutant cDNA data.
CFTR_mutant_cDNA_output = []

# -----------------------------------------------------------------------------
# Step 6A: Create DataFrame of Variants for Downstream Pipeline
# -----------------------------------------------------------------------------

# Generate the wildtype entry using the original Ensembl coordinates.
wt_list = original_Ensembl.copy()
wt_list_modified = [(region[0], region[3], region[4]) for region in wt_list]
wt_cleaned_list = assign_region_names(wt_list_modified)
wt_results = create_mutant_Ensembl_format(wt_cleaned_list)
full_genome_sequence = CFTR2
full_cdna_sequence = build_cdna_sequence(full_genome_sequence, wt_results)
protein_sequence = translate_full_sequence(full_cdna_sequence)
protein_sequence_len = len(protein_sequence)

# Add the wildtype row.
CFTR_mutant_cDNA_output.append({
    "ID": "wildtype",
    "Cleaned Variant": None,
    "Positions": None,
    "Mutation Types": None,
    "Sequence Changes": None,
    "DNA Positions (Mapped)": None,
    "Mutant DNA Sequence": full_genome_sequence,
    "Mapped Regions": None,
    "Final BP": None,
    "Final Format": wt_results,
    "Protein Sequence": protein_sequence,
    "Length of Protein Sequence": protein_sequence_len
})


# -----------------------------------------------------------------------------
# Step 6B: Process Each Variant Entry for Mutant Sequence Generation
# -----------------------------------------------------------------------------

for index, row in df_keep.iterrows():
    variant_name = row['Adjusted Variant cDNA name']
    
    # (i) Clean the variant name and extract mutation details.
    cleaned_variant = clean_cdna_name(variant_name)
    pos, mutation, seq_change, dna_pos = extract_numbers_with_logic(cleaned_variant, result_mapping)

    # (ii) Identify sequence segments from CFTR2 based on the mapped DNA positions.
    sequences = determine_sequence(dna_pos, CFTR2)

    # (iii) Construct the complete mutant DNA sequence by integrating modifications.
    full_genome_sequence = build_sequence(sequences, seq_change)

    # (iv) Mimic Ensembl Format.
    mapped_regions_all = map_regions_for_dna_pos(dna_pos, mutation, seq_change, original_Ensembl)
    result_bps = final_bp_counts(mapped_regions_all, original_Ensembl)
    cleaned_bps = clean_bp_with_indicators(result_bps)
    cleaned_list = assign_region_names(cleaned_bps)
    results = create_mutant_Ensembl_format(cleaned_list)

    # (v) Construct the protein sequence from cDNA sequence.
    full_cdna_sequence = build_cdna_sequence(full_genome_sequence, results)
    protein_sequence = translate_full_sequence(full_cdna_sequence)
    protein_sequence_len = len(protein_sequence)

    CFTR_mutant_cDNA_output.append({
        "ID": variant_name,
        "Cleaned Variant": cleaned_variant,
        "Positions": pos,
        "Mutation Types": mutation,
        "Sequence Changes": seq_change,
        "DNA Positions (Mapped)": dna_pos,
        "Mutant DNA Sequence": full_genome_sequence,
        "Mapped Regions": mapped_regions_all,
        "Final BP": cleaned_list,
        "Final Format": results,
        "Protein Sequence": protein_sequence,
        "Length of Protein Sequence": protein_sequence_len
    })


# -----------------------------------------------------------------------------
# Step 6C: Assign Protein Family and Allele IDs
# -----------------------------------------------------------------------------

CFTR_mutant_cDNA_output = pd.DataFrame(CFTR_mutant_cDNA_output)
CFTR_mutant_cDNA_output = assign_protein_family_and_allele(CFTR_mutant_cDNA_output)


# -----------------------------------------------------------------------------
# Step 6D: Export the Mimic Ensembl Format Data to a .dat File
# -----------------------------------------------------------------------------

df_CFTR_mutant_output = CFTR_mutant_cDNA_output[["ID", "DE", "allele", "Mutant DNA Sequence", "Final Format"]]

# Reorder the columns as desired.
desired_order = [
    "ID", 
    "DE", 
    "allele", 
    "Final Format",
    "Mutant DNA Sequence"
]
df_CFTR_mutant_output = df_CFTR_mutant_output[desired_order]

output_file = "CFTR_Mimic_Ensembl_Format.dat"
export_to_dat(df_CFTR_mutant_output, output_file)


# -----------------------------------------------------------------------------
# Step 6E: Save Reference Data (Variant Name, Legacy Name, and Allele)
# -----------------------------------------------------------------------------

# Create a reference DataFrame with variant IDs and allele identifiers.
df_reference = CFTR_mutant_cDNA_output[['ID', 'allele']].copy()

# Merge with the variant DataFrame to add the legacy name.
df_reference = df_reference.merge(
    df_keep[['Adjusted Variant cDNA name', 'Variant legacy name', 'Class']], 
    how='left', 
    left_on='ID', 
    right_on='Adjusted Variant cDNA name'
)

df_reference.drop(columns=['Adjusted Variant cDNA name'], inplace=True)

# Set the legacy name for the wildtype entry.
df_reference.loc[df_reference['ID'] == 'wildtype', 'Variant legacy name'] = 'wildtype'

# Export the resulting DataFrame as a CSV file
df_reference.to_csv("CFTR_cDNA_Legacy_Allele_Reference.csv", index=False)


# -----------------------------------------------------------------------------
# Step 6F: Save Reference Protein Data (Protein Family, Protein Sequence)
# -----------------------------------------------------------------------------

# Create a reference DataFrame with Protein Family ID and Protein Sequence.
df_protein_reference = CFTR_mutant_cDNA_output[['protein_family_id', 'Protein Sequence', 'Length of Protein Sequence']].copy()
df_protein_reference = df_protein_reference.drop_duplicates(subset='protein_family_id')
df_protein_reference = df_protein_reference.rename(columns={'protein_family_id': 'Protein Family ID'})
# Export the resulting DataFrame as a CSV file
df_protein_reference.to_csv("CFTR_Protein_Family_Reference.csv", index=False)


