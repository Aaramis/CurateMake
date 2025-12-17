# utils.py
import pandas as pd
from typing import List

def create_fasta_header(row: pd.Series, tax_columns: List[str], id_col: str = "seq_id") -> str:
    """
    Crée un en-tête FASTA complet à partir d'une ligne de DataFrame.
    Format : >ID|k__kingdom;p__phylum;c__class;...

    Args:
        row: Une Series pandas représentant une séquence et sa taxonomie.
        tax_columns: Colonnes taxonomiques
        id_col: Nom de la colonne d'identifiant

    Returns:
        Un en-tête formaté pour un fichier FASTA.
    """
    # Récupère l'identifiant de base de la séquence
    base_id = row[id_col]
    
    taxonomy_parts = []
    for col_name in tax_columns:
        # Vérifie si la colonne existe et si sa valeur n'est pas nulle
        if col_name in row and pd.notna(row[col_name]):
            value = row[col_name]
            # Ajoute la partie au format 'nom_colonne=valeur'
            taxonomy_parts.append(f"{col_name}={value}")

    # Assemble la chaîne de taxonomie avec des '|'
    if taxonomy_parts:
        taxonomy_string = "|".join(taxonomy_parts)
        # Combine l'ID et la taxonomie pour l'en-tête final
        return f">{base_id}|{taxonomy_string}"
    else:
        # Si aucune information taxonomique n'est trouvée, retourne juste l'ID
        return f">{base_id}"


def get_fasta_ids(fasta_path: str) -> list:
    """Lit un fichier FASTA et retourne une liste des en-têtes de séquences."""
    ids = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extrait l'ID (le premier mot après '>')
                    ids.append(line[1:].split()[0].strip())
    except FileNotFoundError:
        print(f"Avertissement: Fichier non trouvé {fasta_path}")
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
                    logger.warning(f"Incohérence dans la longueur des séquences dans {msa_path}")
                    break
            
            return num_sequences, alignment_length
    except Exception as e:
        logger.error(f"Erreur lors de l'analyse du fichier MSA {msa_path}: {str(e)}")
        return 0, 0