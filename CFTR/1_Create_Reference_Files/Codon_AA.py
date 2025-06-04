"""
codon_lookup.py

Simple script to map a user-supplied 3-letter codon (e.g. "TGG") 
to its corresponding single-letter amino acid code (e.g. "W").
"""

def translate_codon(codon: str) -> str:
    """
    Return the single-letter amino acid code that corresponds 
    to the given codon, or 'Stop' if it's a stop codon.
    """
    # Standard genetic code (DNA) codon table
    codon_table = {
        # Phenylalanine
        "TTT": "F", "TTC": "F",
        # Leucine
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        # Isoleucine
        "ATT": "I", "ATC": "I", "ATA": "I",
        # Methionine (Start codon)
        "ATG": "M",
        # Valine
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        # Serine
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", 
        "AGT": "S", "AGC": "S",
        # Proline
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        # Threonine
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        # Alanine
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        # Tyrosine
        "TAT": "Y", "TAC": "Y",
        # Histidine
        "CAT": "H", "CAC": "H",
        # Glutamine
        "CAA": "Q", "CAG": "Q",
        # Asparagine
        "AAT": "N", "AAC": "N",
        # Lysine
        "AAA": "K", "AAG": "K",
        # Aspartic acid
        "GAT": "D", "GAC": "D",
        # Glutamic acid
        "GAA": "E", "GAG": "E",
        # Cysteine
        "TGT": "C", "TGC": "C",
        # Tryptophan
        "TGG": "W",
        # Arginine
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGA": "R", "AGG": "R",
        # Glycine
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        # STOP codons
        "TAA": "X", "TAG": "X", "TGA": "X"
    }

    return codon_table[codon]