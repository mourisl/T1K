import re
import pandas as pd
from Codon_AA import translate_codon


"""
Module: VariantMappingAndMutantEnsemblFormatUtils
Description:
    This module provides a comprehensive set of utilities for the integration and 
    standardized formatting of variant data for downstream analyses. It includes functions to:
      - Map cDNA variant numbers to genomic DNA positions.
      - Clean and extract structured numeric information from cDNA variant names.
      - Check for overlaps in mutation regions.
      - Map genomic regions and adjust base pair counts.
      - Construct Ensembl-format data for mutant variants.
      - Clean and combine base pair counts with region indicators.
      - Assign standardized region names.
      - Build mutant Ensembl coordinate formats.
      - Add allele identifier columns to a DataFrame.
      - Export a .dat file to mimic Ensembl format.
"""

# -----------------------------------------------------------------------------
# Constants+++
# -----------------------------------------------------------------------------
MUTATION_KEYWORDS = ['del', '>', 'ins', 'dup']

SPECIAL_CASES = {
    "(?_1)": "-70",
    "(?_-1)": "-70",
    "(*1_?)": "6000"
}

# -----------------------------------------------------------------------------
# Function: map_cdna_to_dna
# -----------------------------------------------------------------------------
def map_cdna_to_dna(cdna_number, mappings):
    """
    Map a cDNA position to the corresponding genomic DNA position using a provided mapping.
    
    Parameters:
        cdna_number: A numeric value (or string) representing the cDNA position.
        mappings: A DataFrame with a 'mapping' column containing dictionaries mapping DNA to RNA positions.
    
    Returns:
        The corresponding DNA position (integer), or None if no mapping is found.
    """
    if pd.isna(cdna_number) or cdna_number in ['', '-']:
        return None

    cdna_str = str(cdna_number).strip()

    # Check for offset indicated by '+' or '-' in the cDNA string.
    if '+' in cdna_str or '-' in cdna_str[1:]:
        match = re.match(r"(-?\d+)([+-]\d+)", cdna_str)
        if match:
            base_number = int(match.group(1))
            offset = int(match.group(2))
            adjusted_number = base_number + 69 if base_number >= 0 else base_number + 70

            for _, mapping_row in mappings.iterrows():
                mapping = mapping_row['mapping']
                if adjusted_number in mapping.values():
                    reversed_mapping = {v: k for k, v in mapping.items()}
                    dna_position = reversed_mapping[adjusted_number] + offset
                    return dna_position
    else:
        try:
            cdna_number = int(cdna_number)
        except ValueError:
            return None

        adjusted_number = cdna_number + 69 if cdna_number >= 0 else cdna_number + 70

        for _, mapping_row in mappings.iterrows():
            mapping = mapping_row['mapping']
            if adjusted_number in mapping.values():
                reversed_mapping = {v: k for k, v in mapping.items()}
                return reversed_mapping[adjusted_number]

    return None


# -----------------------------------------------------------------------------
# Function: add_200_to_numbers
# -----------------------------------------------------------------------------
def add_200_to_numbers(value):
    """
    Increment the given numeric value by 200 if it is valid to consider UTR.
    
    Parameters:
        value: A numeric value (int or float) representing a DNA position.
    
    Returns:
        The value plus 200, or the original value if not a valid number.
    """
    try:
        return value + 200 if pd.notna(value) and isinstance(value, (int, float)) else value
    except Exception:
        return value


# -----------------------------------------------------------------------------
# Function: clean_cdna_name
# -----------------------------------------------------------------------------
def clean_cdna_name(cdna_name):
    """
    Clean a cDNA variant name by replacing defined special cases.
    
    Parameters:
        cdna_name: A string representing the cDNA variant name.
    
    Returns:
        A cleaned cDNA variant name with special cases replaced, or None if input is NaN.
    """
    if pd.isna(cdna_name):
        return None

    clean_name = cdna_name
    for pattern, replacement in SPECIAL_CASES.items():
        clean_name = clean_name.replace(pattern, f"({replacement}_{replacement})")
    return clean_name


# -----------------------------------------------------------------------------
# Function: clean_number
# -----------------------------------------------------------------------------
def clean_number(number):
    """
    Clean a numeric string by removing non-numeric prefixes and extracting the valid portion.
    
    Parameters:
        number: A string containing a numeric value.
    
    Returns:
        A cleaned numeric string, or None if the input is empty.
    """
    if not number:
        return None

    # Remove any non-numeric prefix while preserving '+' and '-'
    number = re.sub(r'^[^\d\+\-]+', '', number)
    
    # Handle compound numbers (e.g., '273+7982')
    match = re.match(r'(\d+\+\d+)', number)
    if match:
        return match.group(1)
    
    return re.sub(r'[^\d\+\-]', '', number)


# -----------------------------------------------------------------------------
# Function: extract_numbers_with_logic
# -----------------------------------------------------------------------------
def extract_numbers_with_logic(cdna_name, result_mapping):
    """
    Extract structured numeric information from a cDNA variant name and map it to DNA positions
    with an additional +200 offset.
    
    The function handles multiple mutation groups (separated by ';') and returns lists containing:
        - Extracted numeric values from cDNA variant name (as a pair for each group)
        - Mutation types
        - Sequence changes
        - Mapped DNA positions after adjustments.
    
    Parameters:
        cdna_name: A string representing the cDNA variant name.
        result_mapping: The mapping between DNA (0-indexed) and RNA positions
    
    Returns:
        A tuple of four lists: (pos, mutation, seq_change, dna_pos)
    """
    if pd.isna(cdna_name):
        return [], [], [], []

    groups = re.split(r';', str(cdna_name))
    pos = [[] for _ in range(len(groups))]
    mutation = [[] for _ in range(len(groups))]
    seq_change = [[] for _ in range(len(groups))]
    dna_pos = [[] for _ in range(len(groups))]

    for i, group in enumerate(groups):
        num_0, num_1 = None, None
        mutation_type, sequence_change = None, ''

        # Attempt to match complex pattern with two numeric pairs
        match = re.match(r'.*\(([^_]+)_([^_]+)\)_\(([^_]+)_([^_]+)\)', group)
        if match:
            num_0 = clean_number(match.group(1))
            num_1 = clean_number(match.group(4))
        else:
            # Handle simple underscore-separated numeric parts
            ranges = re.split(r'_', group)
            if len(ranges) >= 2:
                num_0 = clean_number(ranges[0])
                num_1 = clean_number(ranges[1])
            elif len(ranges) == 1:
                num_0 = clean_number(ranges[0])
                num_1 = None
        if num_1 is None:
            num_1 = num_0

        # Determine mutation type and sequence change based on mutation keywords
        mutation_count = sum(1 for keyword in MUTATION_KEYWORDS if keyword in group)
        if mutation_count > 1 and "delins" in group:
            mutation_type = "delins"
            match = re.search(r"delins([A-Z]+)", group)
            sequence_change = match.group(1) if match else ''
        else:
            if "del" in group:
                mutation_type = "del"
                sequence_change = ''
            elif ">" in group:
                mutation_type = "mut"
                sequence_change = group.split(">")[1][0]
            elif "ins" in group:
                mutation_type = "ins"
                match = re.search(r"ins([A-Z]+)", group)
                sequence_change = match.group(1) if match else ''
            elif "dup" in group:
                mutation_type = "dup"
                sequence_change = ''
                if num_0 is not None and num_1 is not None:
                    num_0, num_1 = num_1, num_0

        # Map extracted numbers to DNA positions using the provided mapping
        mapped_num0 = map_cdna_to_dna(num_0, result_mapping)
        mapped_num1 = map_cdna_to_dna(num_1, result_mapping)

        # Adjust positions for certain mutation types
        if mutation_type in ["delins", "del", "mut"]:
            if mapped_num0 is not None:
                mapped_num0 -= 1
            if mapped_num1 is not None:
                mapped_num1 += 1

        # Apply a +200 offset to the mapped positions
        mapped_num0 = add_200_to_numbers(mapped_num0)
        mapped_num1 = add_200_to_numbers(mapped_num1)

        pos[i] = [num_0, num_1]
        mutation[i].append(mutation_type)
        seq_change[i].append(sequence_change)
        dna_pos[i] = [mapped_num0, mapped_num1]

    return pos, mutation, seq_change, dna_pos


# -----------------------------------------------------------------------------
# Function: check_overlap_mutation_regions
# -----------------------------------------------------------------------------
def check_overlap_mutation_regions(dna_pos, mutation_types, cdna_names, 
                                   variants_to_drop=None, variants_to_keep=None,
                                   adjusted_variant_names=None):
    """
    Evaluate mutation regions for overlaps within each cDNA variant and produce adjusted names.
    
    For each variant, the function:
      - Extracts numeric ranges from the cDNA variant name,
      - Adjusts the DNA position ranges based on mutation type rules,
      - Checks for overlaps among mutation regions.
      
    Variants with overlapping regions are flagged for exclusion, while non-overlapping variants
    are retained with adjusted names standardized as 'c.[mutation1;mutation2;...]'.
    
    Parameters:
        dna_pos: List of lists containing mapped DNA position ranges for each mutation group.
        mutation_types: List of lists containing mutation types for each mutation group.
        cdna_names: List of original cDNA variant names.
        variants_to_drop: List to store variants with overlapping regions.
        variants_to_keep: List to store variants without overlaps.
        adjusted_variant_names: Dictionary to store standardized variant names.
    
    Returns:
        A tuple: (variants_to_drop, variants_to_keep, adjusted_variant_names)
    """
    if variants_to_drop is None:
        variants_to_drop = []
    if variants_to_keep is None:
        variants_to_keep = []
    if adjusted_variant_names is None:
        adjusted_variant_names = {}

    for idx, positions_group in enumerate(dna_pos):
        ranges = []
        # Clean and split the cDNA variant name
        mutations = cdna_names[idx].replace('c.[', '').replace(']', '').split(';')
        mutations = [m.replace('c.', '') for m in mutations]

        # Adjust positions based on mutation type rules
        for j, positions_list in enumerate(positions_group):
            pos0, pos1 = positions_list[0], positions_list[-1]
            mt = mutation_types[idx][j] if isinstance(mutation_types[idx][j], str) else mutation_types[idx][j][0]

            if mt in ["delins", "del", "mut"]:
                pos0 += 1
                pos1 -= 1
            elif mt == "dup":
                pos0, pos1 = pos1, pos0

            start, end = sorted([pos0, pos1])
            ranges.append((start, end))

        # Check for overlaps within the variant
        overlap_found = False
        ranges_sorted = sorted(ranges, key=lambda x: x[0])
        for i in range(len(ranges_sorted) - 1):
            current_end = ranges_sorted[i][1]
            next_start = ranges_sorted[i + 1][0]
            if current_end >= next_start:
                overlap_found = True
                break

        if overlap_found:
            variants_to_drop.append(cdna_names[idx])
        else:
            variants_to_keep.append(cdna_names[idx])
            # Sort mutations by starting position
            mutations_with_pos = sorted(zip(mutations, ranges), key=lambda x: x[1][0])
            sorted_mutations = [m[0] for m in mutations_with_pos]
            adjusted_variant_names[cdna_names[idx]] = f'c.[{";".join(sorted_mutations)}]'

    return variants_to_drop, variants_to_keep, adjusted_variant_names


# -----------------------------------------------------------------------------
# Function: determine_sequence
# -----------------------------------------------------------------------------
def determine_sequence(dna_pos, CFTR2):
    """
    Identify sequence segments from the full CFTR2 DNA sequence based on specified region boundaries.
    
    This function partitions the CFTR2 sequence into segments using the provided positional ranges.
    These segments can later be used to reconstruct the complete sequence with appropriate modifications.

    Parameters:
        dna_pos (list of tuple): A list of tuples, each tuple containing the start and end 
                                 positions (e.g., [(start1, end1), (start2, end2), ...]).
        CFTR2 (str): The complete CFTR2 DNA sequence.

    Returns:
        list: A list of sequence segments extracted from CFTR2.
    """
    sequences = []

    # Starting sequence segment: from start of CFTR2 to pos1 of the first part
    if dna_pos[0][0] is not None:
        seq0 = CFTR2[:dna_pos[0][0]+1]
        sequences.append(seq0)

    # Loop through pairs in pos to define intermediate sequence segments
    for i in range(len(dna_pos) - 1):
        pos2_current = dna_pos[i][1]  # pos2 of the current part
        pos1_next = dna_pos[i + 1][0]  # pos1 of the next part

        if pos2_current is not None and pos1_next is not None:
            seq_i = CFTR2[pos2_current : pos1_next + 1]
            sequences.append(seq_i)

    # Last sequence segment: from pos2 of the last part to the end of CFTR2
    if dna_pos[-1][1] is not None:
        last_seq = CFTR2[dna_pos[-1][1]:]
        sequences.append(last_seq)
    
    return sequences


# -----------------------------------------------------------------------------
# Function: build_sequence
# -----------------------------------------------------------------------------
def build_sequence(sequences, seq_change):
    """
    Construct the final DNA sequence by integrating original sequence segments with modifications.

    This function concatenates each segment from the 'sequences' list with its corresponding
    modified sequence from 'seq_change'. If no segments are provided, the function returns the 
    complete CFTR2 sequence. In that case, the full CFTR2 sequence must be provided via the optional
    CFTR2 parameter.

    Parameters:
        sequences (list): A list of sequence segments extracted from the full CFTR2 sequence.
        seq_change (list): A list of modifications, where each element is a list whose first element
                           contains the inserted/modified sequence.

    Returns:
        str: The fully assembled DNA sequence after applying all modifications.
    """
    DNA_sequence = ""

    # Check if sequences is None or empty
    if not sequences:
        # If sequences is None, use the full CFTR2 sequence
        DNA_sequence = CFTR2
    else:
        # Loop through each segment in sequences, up to len(seq_change)
        for i in range(len(seq_change)):
            # Add the original sequence segment
            DNA_sequence += sequences[i]

            # Add the corresponding changed_seq
            DNA_sequence += seq_change[i][0]

        # Append the last segment of sequences
        DNA_sequence += sequences[-1]

    return DNA_sequence


# -----------------------------------------------------------------------------
# Function: get_region_for_position
# -----------------------------------------------------------------------------
def get_region_for_position(pos, ensembl_data):
    """
    Retrieve the region identifier from ensembl_data for a given position.
    
    Parameters:
        pos (int): The genomic position.
        ensembl_data (list of tuples): Each tuple contains 
            (region_identifier, start, end, base_pair_count, region_indicator).
    
    Returns:
        str or None: The region identifier if pos is within a region; otherwise, None.
    """
    for (region, start, end, num_bp, indicator) in ensembl_data:
        if start <= pos <= end:
            return region
    return None


# -----------------------------------------------------------------------------
# Function: map_regions_for_dna_pos
# -----------------------------------------------------------------------------
def map_regions_for_dna_pos(dna_pos, mutation_types, seq_change, ensembl_data):
    """
    Map DNA positions to Ensembl regions and calculate the net base pair change.
    
    For each mutation, adjust the start and end positions based on the mutation type,
    then compute the base pair change and identify the affected regions.
    
    Parameters:
        dna_pos (list): A list of position pairs (each a list/tuple of positions).
        mutation_types (list): A list of mutation types corresponding to each region.
        seq_change (list): A list of sequence change modifications.
        ensembl_data (list of tuples): Ensembl region data with structure 
            (region, start, end, num_bp, indicator).
    
    Returns:
        list: A list of mappings in the form:
              [affected_region0, affected_region1, net_bp_change, pos0, pos1, mutation_type]
    """
    mapped_regions_all = []

    for i, positions_list in enumerate(dna_pos):
        pos0, pos1 = positions_list[0], positions_list[-1]
        mt = mutation_types[i] if isinstance(mutation_types[i], str) else mutation_types[i][0]

        # Adjust positions based on mutation type
        if mt in ["delins", "del", "mut"]:
            pos0 += 1
            pos1 -= 1
        elif mt == "dup":
            pos0, pos1 = pos1, pos0

        # Calculate base pair change based on mutation type
        if mt == "delins":
            num_bp_changed = len(seq_change[i][0]) - (pos1 - pos0 + 1)
        elif mt == "dup":
            num_bp_changed = pos1 - pos0 + 1
        elif mt == "del":
            num_bp_changed = -(pos1 - pos0 + 1)
        elif mt == "ins":
            num_bp_changed = len(seq_change[i][0])
        else:
            num_bp_changed = 0

        affected_region0 = get_region_for_position(pos0, ensembl_data)
        affected_region1 = get_region_for_position(pos1, ensembl_data)
        mapped_regions_all.append([affected_region0, affected_region1, num_bp_changed, pos0, pos1, mt])
    
    return mapped_regions_all


# -----------------------------------------------------------------------------
# Function: final_bp_counts
# -----------------------------------------------------------------------------
def final_bp_counts(mapped_regions_all, ensembl_data):
    """
    Compute the final base pair counts for each Ensembl region after mutation adjustments.
    
    The function updates base pair counts based on the mutation mapping and handles
    cases such as deletions, insertions, duplications, and delins. If a mutation spans
    multiple regions, the counts are adjusted accordingly.
    
    Parameters:
        mapped_regions_all (list): Output list from map_regions_for_dna_pos.
        ensembl_data (list of tuples): Ensembl region data with structure 
            (region, start, end, base_pair_count, indicator).
    
    Returns:
        list: A list of tuples in the form (region, final_bp_count, region_indicator) for regions
              with non-zero base pair counts, or a message if further exploration is needed.
    """
    region_names = [r for (r, s, e, nb, ind) in ensembl_data]
    final_bp = [nb for (r, s, e, nb, ind) in ensembl_data]
    indicators = [ind for (r, s, e, nb, ind) in ensembl_data]
    region_boundaries = {r: (s, e) for (r, s, e, nb, ind) in ensembl_data}
    
    for mapping in mapped_regions_all:
        region0, region1, bp_change, pos0, pos1, mt = mapping
        
        if region0 == region1 and region0 is not None:
            if region0 in region_names:
                idx = region_names.index(region0)
                final_bp[idx] += bp_change
            else:
                return "Needs more exploratory"
        else:
            # When the mutation spans two different regions
            left_bound = region_boundaries.get(region0, (None, None))[0]
            right_bound = region_boundaries.get(region1, (None, None))[1]
            if left_bound is None or right_bound is None:
                return "Needs more exploratory"
            idx0 = region_names.index(region0)
            idx1 = region_names.index(region1)

            if pos0 == left_bound and pos1 == right_bound and mt == "del":
                # Set bp count to 0 for regions spanning the deletion
                for idx in range(idx0, idx1 + 1):
                    final_bp[idx] = 0
            elif pos0 == left_bound and pos1 == right_bound and mt == "dup":
                duplicate_entries = [
                    (region_names[idx], final_bp[idx], indicators[idx])
                    for idx in range(idx0, idx1 + 1)
                ]
                original_entries = [
                    (region_names[i], final_bp[i], indicators[i])
                    for i in range(len(final_bp))
                ]
                first_segment = original_entries[:idx1 + 1]
                second_segment = original_entries[idx1 + 1:]
                combined_entries = first_segment + duplicate_entries + second_segment
                region_names = [entry[0] for entry in combined_entries]
                final_bp = [entry[1] for entry in combined_entries]
                indicators = [entry[2] for entry in combined_entries]
            elif mt == "del":
                new_val_region0 = pos0 - left_bound
                new_val_region1 = right_bound - pos1
                final_bp[idx0] = new_val_region0
                final_bp[idx1] = new_val_region1
                for idx in range(idx0 + 1, idx1):
                    final_bp[idx] = 0
            elif mt == "ins":
                ind0 = indicators[region_names.index(region0)] if region0 in region_names else None
                ind1 = indicators[region_names.index(region1)] if region1 in region_names else None
                if ind0 == "exon":
                    idx = region_names.index(region0)
                    final_bp[idx] = final_bp[idx] + bp_change
                elif ind1 == "exon":
                    idx = region_names.index(region1)
                    final_bp[idx] = final_bp[idx] + bp_change
                else:
                    return "Needs more exploratory"
            elif mt == "delins":
                new_val_region0 = pos0 - left_bound
                new_val_region1 = right_bound - pos1
                final_bp[idx0] = new_val_region0
                final_bp[idx1] = new_val_region1
                for idx in range(idx0 + 1, idx1):
                    final_bp[idx] = 0
                ind0 = indicators[region_names.index(region0)] if region0 in region_names else None
                ind1 = indicators[region_names.index(region1)] if region1 in region_names else None
                if ind0 == "exon" and ind1 != "exon":
                    idx_target = region_names.index(region0)
                elif ind1 == "exon" and ind0 != "exon":
                    idx_target = region_names.index(region1)
                elif ind0 == "exon" and ind1 == "exon":
                    idx_target = region_names.index(region0)
                else:
                    idx_target = region_names.index(region0)
                final_bp[idx_target] += bp_change + pos1 - pos0 + 1
            elif mt == "dup":
                left_bound_region1 = region_boundaries[region1][0]
                dup_bp_region1 = pos1 - left_bound_region1 + 1
                dup_entry_region1 = (region1, dup_bp_region1, indicators[region_names.index(region1)])
                right_bound_region0 = region_boundaries[region0][1]
                dup_bp_region0 = right_bound_region0 - pos0 + 1
                dup_entry_region0 = (region0, dup_bp_region0, indicators[region_names.index(region0)])
                duplicate_entries = [dup_entry_region1, dup_entry_region0]
                original_entries = [
                    (region_names[i], final_bp[i], indicators[i])
                    for i in range(len(final_bp))
                ]
                first_segment = original_entries[:idx1]
                second_segment = original_entries[idx1 - 1:]
                combined_entries = first_segment + duplicate_entries + second_segment
                region_names = [entry[0] for entry in combined_entries]
                final_bp = [entry[1] for entry in combined_entries]
                indicators = [entry[2] for entry in combined_entries]
            else:
                return "Needs more exploratory"

    final_bp_with_indicators = [
        (region_names[i], bp, indicators[i]) 
        for i, bp in enumerate(final_bp) if bp != 0
    ]
    return final_bp_with_indicators


# -----------------------------------------------------------------------------
# Function: clean_bp_with_indicators
# -----------------------------------------------------------------------------
def clean_bp_with_indicators(final_bp_with_indicators):
    """
    Merge adjacent entries with the same region indicator to simplify the base pair counts.
    
    Parameters:
        final_bp_with_indicators (list): A list of tuples (region, bp_count, indicator).
    
    Returns:
        list: A cleaned list of tuples with merged base pair counts.
    """
    if not final_bp_with_indicators:
        return []

    cleaned = []
    current_region, current_bp, current_ind = final_bp_with_indicators[0]

    for entry in final_bp_with_indicators[1:]:
        region, bp, ind = entry
        if ind == current_ind:
            current_bp += bp
        else:
            cleaned.append((current_region, current_bp, current_ind))
            current_region, current_bp, current_ind = region, bp, ind

    cleaned.append((current_region, current_bp, current_ind))
    return cleaned


# -----------------------------------------------------------------------------
# Function: assign_region_names
# -----------------------------------------------------------------------------
def assign_region_names(cleaned):
    """
    Assign standardized region names to each cleaned entry.
    
    The first entry is labeled as UTR, intermediate entries are alternately 
    labeled as exon or intron (with sequential numbering), and the last entry is UTR.
    
    Parameters:
        cleaned (list): A list of tuples (region, bp_count, indicator).
    
    Returns:
        list: A list of tuples (region_label, bp_count) with assigned region names.
    """
    n = len(cleaned)
    if n == 0:
        return []
    if n == 1:
        # Only one entry: assign it as UTR.
        _, bp, _ = cleaned[0]
        return [("UTR", bp)]
    
    new_list = []
    # The first entry is always UTR.
    _, bp, _ = cleaned[0]
    new_list.append(("UTR", bp))
    
    # Intermediate entries: assign alternating exon/intron labels.
    for i in range(1, n - 1):
        _, bp, _ = cleaned[i]
        j = i - 1  # Zero-indexed adjustment.
        if j % 2 == 0:
            region_label = f"exon{(j // 2) + 1}"
        else:
            region_label = f"intron{(j // 2) + 1}"
        new_list.append((region_label, bp))
    
    # The last entry is always UTR.
    _, bp, _ = cleaned[-1]
    new_list.append(("UTR", bp))
    
    return new_list


# -----------------------------------------------------------------------------
# Function: create_mutant_Ensembl_format
# -----------------------------------------------------------------------------
def create_mutant_Ensembl_format(new_list):
    """
    Generate the mutant Ensembl coordinate format based on region names and base pair counts.
    
    Each region in new_list is assigned a start and end coordinate in a continuous system.
    
    Parameters:
        new_list (list): A list of tuples (region_label, base_pair_count).
    
    Returns:
        list: A list of tuples (region_label, start_coordinate, end_coordinate, base_pair_count).
    """
    pos0 = 0
    results = []
    
    for region_label, length in new_list:
        pos1 = pos0 + length - 1
        results.append((region_label, pos0, pos1, length))
        pos0 = pos1 + 1
        
    return results


# -----------------------------------------------------------------------------
# Function: build_cdna_sequence
# -----------------------------------------------------------------------------
def build_cdna_sequence(full_genome_sequence, mapped_regions):
    """
    Extracts and concatenates exonic regions from the full genome sequence.

    Parameters:
    - full_genome_sequence (str): The complete mutant genomic DNA sequence.
    - mapped_regions (list of tuples): List of tuples of the form 
      (region_label, pos0, pos1, length). Regions with 'exon' in their label will be extracted.

    Returns:
    - full_cdna_sequence (str): The concatenated cDNA sequence composed of exonic segments.
    """
    exon_segments = []
    
    for region_label, pos0, pos1, length in mapped_regions:
        # Check if the region label contains 'exon' (case-insensitive)
        if 'exon' in region_label.lower():
            # Extract the segment. Note that Python slicing is zero-indexed
            # and is [start, end) so we use pos1 + 1 to include pos1.
            exon_seq = full_genome_sequence[pos0:pos1+1]
            exon_segments.append(exon_seq)
    
    full_cdna_sequence = ''.join(exon_segments)
    full_cdna_sequence = full_cdna_sequence[70:]
    return full_cdna_sequence


# -----------------------------------------------------------------------------
# Function: translate_full_sequence
# -----------------------------------------------------------------------------

def translate_full_sequence(full_cdna_sequence: str) -> str:
    """
    Translates the full cDNA sequence into the corresponding amino acid sequence.
    
    The function groups the sequence into codons (triplets), translates each codon 
    using the translate_codon function, and stops translating when the stop codon 'X' is encountered.
    
    Parameters:
    - full_cdna_sequence (str): The complete cDNA sequence.
    
    Returns:
    - protein_sequence (str): The concatenated amino acid sequence.
    """
    protein_sequence = ""
    
    # Iterate over the sequence in steps of 3 (one codon at a time)
    for i in range(0, len(full_cdna_sequence), 3):
        codon = full_cdna_sequence[i:i+3]
        
        # If we get an incomplete codon (less than 3 bases), break out of the loop
        if len(codon) < 3:
            break
        
        aa = translate_codon(codon)
        
        # If a stop codon is encountered, stop translation immediately.
        if aa == "X":
            break
        
        protein_sequence += aa
        
    return protein_sequence


# -----------------------------------------------------------------------------
# Function: assign_protein_family_and_allele
# -----------------------------------------------------------------------------
def assign_protein_family_and_allele(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assigns protein family and allele identifiers to each record based on the 
    'Protein Sequence' column, then formats these as a standardized identifier.
    
    This function computes a unique protein family ID using pd.factorize on the 
    'Protein Sequence' column. Within each protein family, a cumulative count is 
    used to assign an allele ID. These two values are combined into a formatted 
    identifier of the form "CFTR*XXXX:YYYY", where XXXX and YYYY are zero-padded 
    to 4 digits. An additional column 'allele' is created as an alias of 'DE'.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least the 'Protein Sequence' column.
    
    Returns
    -------
    pd.DataFrame
        The original DataFrame augmented with the following columns:
          - protein_family_id: Unique ID for each protein sequence (starting from 1).
          - allele_id: The cumulative allele count within each protein family.
          - DE: Formatted identifier (e.g., "CFTR*0001:0001").
          - allele: Duplicate of DE, to serve as an alias.
    """
    # Assign a unique protein family ID based on the 'Protein Sequence' column.
    df['protein_family_id'] = pd.factorize(df['Protein Sequence'])[0] + 1

    # Within each protein family group, assign allele IDs as the cumulative occurrence count.
    df['allele_id'] = df.groupby('protein_family_id').cumcount() + 1

    # Create the formatted identifier in the 'DE' column.
    df['DE'] = df.apply(
        lambda row: f"CFTR*{row['protein_family_id']:04d}:{row['allele_id']:04d}", axis=1
    )
    df['allele'] = df['DE']

    return df


# -----------------------------------------------------------------------------
# Function: export_to_dat
# -----------------------------------------------------------------------------
def export_to_dat(df, output_file):
    """
    Export variant data to a .dat file mimicking the Ensembl format.

    For each variant in the DataFrame, this function writes header information 
    (ID, DE, and allele) followed by feature table (FT) lines for each region 
    specified in the "Final Format" column. Regions are output with 1-indexed coordinates,
    and the mutant DNA sequence is formatted and appended at the end.

    Parameters:
        df (DataFrame): DataFrame containing variant data.
        output_file (str): The output file name for the .dat export.
    """
    with open(output_file, "w") as f:
        for idx, row in df.iterrows():
            # Header lines
            f.write("ID\t" + str(row["ID"]) + "\n")
            f.write("DE\t" + str(row["DE"]) + "\n")
            f.write(f'FT\t/allele="{row["allele"]}"\n')
            
            # Process each entry in the Final Format data
            final_format_data = row["Final Format"]
            for entry in final_format_data:
                # Determine region_label, pos0, pos1 from entry.
                if isinstance(entry, dict):
                    region_label = entry.get("region_label", "")
                    pos0 = entry.get("pos0", "")
                    pos1 = entry.get("pos1", "")
                elif isinstance(entry, (list, tuple)) and len(entry) >= 3:
                    region_label, pos0, pos1 = entry[0], entry[1], entry[2]
                else:
                    region_label, pos0, pos1 = "", "", ""
                
                # Convert pos0 and pos1 to 1-indexed values.
                pos0_1 = int(pos0) + 1
                pos1_1 = int(pos1) + 1
                
                # Write FT lines based on region type.
                match = re.match(r'^(exon|intron)(\d+)$', region_label, re.IGNORECASE)
                if match:
                    feature_type = match.group(1) 
                    number = match.group(2) 
                    f.write("FT\t" + f"{feature_type:<15}" + f"{pos0_1}..{pos1_1}" + "\n")
                    f.write("FT\t" + " " * 15 + f'/number="{number}"' + "\n")
                elif region_label.upper() == "UTR":
                    f.write("FT\t" + f"{region_label:<15}" + f"{pos0_1}..{pos1_1}" + "\n")
                    last_utr_pos = pos1_1
                    
            if last_utr_pos is not None:
                mutant_seq = str(row["Mutant DNA Sequence"]).lower()
                formatted_lines = format_dna_sequence(mutant_seq, last_utr_pos, chunk_size=60, group_size=10)
                for line in formatted_lines:
                    f.write(line + "\n")
            
            f.write("//\n")


# -----------------------------------------------------------------------------
# Function: format_dna_sequence
# -----------------------------------------------------------------------------
def format_dna_sequence(seq, total_bp, chunk_size=60, group_size=10):
    """
    Format the DNA sequence to mimic the Ensembl .dat sequence block.

    The function splits the sequence into chunks (default: 60 bases per line)
    and groups bases into subgroups (default: 10 bases per group), appending a 
    running total of bases at the end of each line.

    Parameters:
        seq (str): The DNA sequence to be formatted.
        total_bp (int): Total base pair count to be displayed in the header.
        chunk_size (int, optional): Number of bases per line (default is 60).
        group_size (int, optional): Number of bases per group (default is 10).

    Returns:
        list: A list of formatted string lines representing the sequence block.
    """
    lines = []
    
    # First line: "SQ   Sequence 14738 BP;"
    lines.append(f"SQ\tSequence {total_bp} BP;")
    
    prefix = "        "
    
    i = 0
    total_count = 0
    while i < len(seq):
        chunk = seq[i : i + chunk_size]
        i += chunk_size
        total_count += len(chunk)
        
        # Group the chunk into subgroups.
        groups = [chunk[j : j + group_size] for j in range(0, len(chunk), group_size)]
        grouped_str = " ".join(groups)
        line_str = f"{prefix}{grouped_str:<65}{str(total_count).rjust(8)}"
        
        lines.append(line_str)
    
    return lines

