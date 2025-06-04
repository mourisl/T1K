import pandas as pd



def generate_result_mapping(df_pos):
    """
    Generate DNA-RNA mapping from position dataframe.

    Returns:
    - DataFrame with 'exon-intron_structure' and 'mapping' columns
    """

    # Create range for DNA
    df_pos['range'] = df_pos.apply(
        lambda row: list(range(row['pos1_CFTR_py'], row['pos2_CFTR_py'] + 1)), axis=1
    )

    # Filter only exons for RNA range
    df_exon_rna = df_pos[df_pos['type'] == 'exon'].copy()
    df_exon_rna['range_RNA'] = df_exon_rna.apply(
        lambda row: list(range(int(row['pos1_RNA']) - 1, int(row['pos2_RNA']))), axis=1
    )

    # Merge DNA and RNA ranges for mapping
    merged = pd.merge(
        df_pos[['exon-intron_structure', 'range']],
        df_exon_rna[['exon-intron_structure', 'range_RNA']],
        on='exon-intron_structure'
    )
    merged['mapping'] = merged.apply(
        lambda row: {dna: rna for dna, rna in zip(row['range'], row['range_RNA'])}, axis=1
    )

    return merged[['exon-intron_structure', 'mapping']] 