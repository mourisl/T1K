def create_original_Ensembl_format(num_of_bp):

    """
    Generate contiguous genomic coordinate ranges based on specified base pair counts.

    For each region length in the provided list, this function computes the start (pos0) 
    and end (pos1) coordinates for that region in a continuous genomic coordinate system.
    
    - The first region (region0) starts at pos0 = 0 and ends at pos1 = pos0 + length - 1.
    - For each subsequent region i, the start coordinate is defined as:
          pos0 = (previous region's pos1 + 1).

    Each region is assigned an indicator:
        - 'intron' for even-indexed regions,
        - 'exon' for odd-indexed regions.

    Parameters:
        num_of_bp (list of int): A list of base pair counts for each region.

    Returns:
        List[tuple]: A list of tuples in the format:
                     (region_identifier, start_coordinate, end_coordinate, base_pair_count, region_indicator)
    """

    pos0 = 0
    results = []
    
    for i, length in enumerate(num_of_bp):
        pos1 = pos0 + length - 1
        indicator = 'intron' if i % 2 == 0 else 'exon'
        results.append((f"region{i}", pos0, pos1, length, indicator))
        pos0 = pos1 + 1
        
    return results
