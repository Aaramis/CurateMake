# utils.py
import pandas as pd
from typing import List

def create_fasta_header(row: pd.Series, tax_columns: List[str], id_col: str = "seq_id") -> str:
    """
    Build a FASTA header from a DataFrame row.
    Format: >ID|kingdom=...|phylum=...|...

    Args:
        row: pandas Series representing a sequence + taxonomy.
        tax_columns: Taxonomy columns to include.
        id_col: Identifier column name (default seq_id).

    Returns:
        FASTA header string.
    """
    # Grab the base identifier
    base_id = row[id_col]
    
    taxonomy_parts = []
    for col_name in tax_columns:
        # Ensure column exists and value is not NaN
        if col_name in row and pd.notna(row[col_name]):
            value = row[col_name]
            # Append as column=value
            taxonomy_parts.append(f"{col_name}={value}")

    # Assemble taxonomy chain
    if taxonomy_parts:
        taxonomy_string = "|".join(taxonomy_parts)
        # Combine ID + taxonomy
        return f">{base_id}|{taxonomy_string}"
    else:
        # No taxonomy data available
        return f">{base_id}"


def get_fasta_ids(fasta_path: str) -> list:
    """Return the sequence headers from a FASTA file."""
    ids = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Grab the first token after '>'
                    ids.append(line[1:].split()[0].strip())
    except FileNotFoundError:
        print(f"Warning: FASTA file not found {fasta_path}")
    return ids


def analyze_msa_file(msa_path, logger):
    try:
        with open(msa_path, 'r') as f:
            sequences = []
            current_seq = ''
            for line in f:
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                    current_seq = ''
                else:
                    current_seq += line.strip()
            if current_seq:
                sequences.append(current_seq)
            
            if not sequences:
                return 0, 0
            
            num_sequences = len(sequences)
            alignment_length = len(sequences[0]) if sequences else 0
            
            for seq in sequences[1:]:
                if len(seq) != alignment_length:
                    logger.warning(f"Inconsistent alignment length detected in {msa_path}")
                    break
            
            return num_sequences, alignment_length
    except Exception as e:
        logger.error(f"Failed to analyze MSA file {msa_path}: {str(e)}")
        return 0, 0
