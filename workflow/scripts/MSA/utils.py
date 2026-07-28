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
            # Escape spaces: MAFFT/trimAl truncate FASTA headers at the first
            # whitespace, which would otherwise drop the species epithet
            # (e.g. "species=Russula azurea" -> "species=Russula"). Underscores
            # are restored to spaces when the taxonomy is read back
            # (see decode_header_tax_value / SATIVA integration).
            value = str(row[col_name]).replace(" ", "_")
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
