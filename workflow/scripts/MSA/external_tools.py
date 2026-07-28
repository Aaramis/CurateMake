# external_tools.py
import os
import subprocess
import tempfile
import pandas as pd
import utils
from typing import List, Set
from pathlib import Path


def existing_msa_ids(msa_path: str) -> Set[str]:
    """Raw seq_ids already present in an aligned FASTA (header up to first space / '|').

    Used to make ``--add`` integrations idempotent: on a partial re-run the MSA on
    disk may already contain previously-integrated sequences (build_family is cached
    on its manifest, not on the per-MSA FASTA), so re-adding them would duplicate.
    """
    ids: Set[str] = set()
    if not msa_path or not os.path.exists(msa_path):
        return ids
    with open(msa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0].split("|")[0])
    return ids


def run_hmmsearch(sequences_df: pd.DataFrame,
                  hmm_db_path: str,
                  seq_col: str,
                  max_thread_workers: int,
                  tax_columns: List[str],
                  id_col: str = "seq_id",
                  max_evalue: float = 0.001,
                  dom_e_value: float = 0.001) -> pd.DataFrame:
    """
    Classify multiple sequences against a combined HMM profile library using
    hmmsearch (profiles = query, sequences = target). Each sequence is assigned
    its best-scoring profile by E-value. The profile library is plain text; no
    ``hmmpress`` index is required.

    Args:
        sequences_df: DataFrame containing the sequences.
        hmm_db_path: path to the combined HMM file (plain text; no hmmpress needed).
        seq_col: sequence column name.
        max_thread_workers: number of threads passed to hmmsearch.
        tax_columns: taxonomy columns to render FASTA headers.
        id_col: identifier column.
        max_evalue: maximum E-value threshold (enforced by hmmsearch -E and re-checked).
        dom_e_value: domain-level E-value (currently unused, placeholder).

    Returns:
        DataFrame with columns [id_col, 'assigned_msa', 'score', 'e_value'].
    """
    # Create temporary FASTA/tbl files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as fasta_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.tbl', delete=False) as tbl_file:

        # Build FASTA input
        for _, row in sequences_df.iterrows():
            header = utils.create_fasta_header(row, tax_columns, id_col)
            fasta_file.write(f"{header}\n{row[seq_col]}\n")
        fasta_path = fasta_file.name
        tbl_path = tbl_file.name

    try:
        # Execute hmmsearch: query = HMM profile library, target = sequence batch.
        # -E caps the reported hits at the classification threshold (keeps the
        # tblout small); the best hits sit far below it (E ~ 1e-60 in practice).
        command = [
            'hmmsearch',
            '--tblout', tbl_path,
            '--cpu', str(max_thread_workers),
            '--noali',
            '-E', str(max_evalue),        # Full-sequence E-value reporting threshold
            # '--domE', str(dom_e_value),   # Domain-level threshold (currently unused)
            hmm_db_path,
            fasta_path
        ]
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Parse hmmsearch --tblout: target (parts[0]) = seq_id, query (parts[2]) = profile_name.
        hits = []
        with open(tbl_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 18:
                    continue
                seq_id = parts[0].split('|')[0]
                hits.append({
                    id_col: seq_id,  # UUID string
                    'assigned_msa': parts[2],  # Nom du profil HMM
                    'score': float(parts[5]),
                    'e_value': float(parts[4])
                })

        if not hits:
            return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score', 'e_value'])

        df = pd.DataFrame(hits)
        if max_evalue is not None:
            df = df[df['e_value'] <= max_evalue]
        if df.empty:
            return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score', 'e_value'])

        df.sort_values(['e_value', 'score'], ascending=[True, False], inplace=True)
        return df.drop_duplicates(id_col)

    except subprocess.CalledProcessError as e:
        print(f"hmmsearch failed: {e.stderr}")
        return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score'])
    finally:
        # Cleanup temporary files
        for path in [fasta_path, tbl_path]:
            try:
                if path and os.path.exists(path):
                    os.remove(path)
            except OSError:
                pass


def run_mafft_add(new_sequences_df: pd.DataFrame,
                  existing_msa_path: str,
                  output_msa_path: str,
                 seq_col: str,
                 tax_columns: List[str],
                 max_thread_workers: int,
                 id_col: str = "seq_id") -> bool:
    """
    Adds sequences from a DataFrame to an existing MSA file using mafft --add.

    Args:
        new_sequences_df: DataFrame with the sequences to add.
        existing_msa_path: Path to the existing MSA file (in FASTA format).
        output_msa_path: Path to write the updated MSA file.
        seq_col: Sequence column name.
        tax_columns: Taxonomy columns.
        id_col: Identifier column name.
        
    Returns:
        True if successful, False otherwise.
    """
    if not os.path.exists(existing_msa_path):
        print(f"Error: existing MSA file not found: {existing_msa_path}")
        return False

    # Create a temporary file for the new sequences to be added
    with tempfile.NamedTemporaryFile(mode='w+', delete=True, suffix='.fasta', encoding='utf-8') as temp_fasta:
        # Write new sequences to the temporary FASTA file
        for _, row in new_sequences_df.iterrows():
            header = utils.create_fasta_header(row, tax_columns, id_col)
            temp_fasta.write(f"{header}\n{row[seq_col]}\n")
        temp_fasta.flush() # Ensure all data is written to disk

        cmd = [
            'mafft',
            '--add', temp_fasta.name,
            '--thread', str(max_thread_workers),
            '--reorder',
            '--quiet',
            existing_msa_path
        ]

        try:
            # Run the command, capturing the output
            process = subprocess.run(cmd, capture_output=True, text=True, check=True, encoding='utf-8')
            
            # Write the new, combined alignment to the final output path
            with open(output_msa_path, 'w', encoding='utf-8') as f_out:
                f_out.write(process.stdout)
            
            return True

        except subprocess.CalledProcessError as e:
            print(f"mafft --add failed for {os.path.basename(existing_msa_path)}.")
            print(f"stderr: {e.stderr}")
            return False
        except Exception as e:
            print(f"Unexpected error during profile alignment: {str(e)}")
            return False


def run_hmmalign_add(new_sequences_df: pd.DataFrame,
                     existing_msa_path: str,
                     output_msa_path: str,
                     profile_path: str,
                     seq_col: str,
                     tax_columns: List[str],
                     id_col: str = "seq_id") -> bool:
    """
    Add sequences to an existing clade MSA using ``hmmalign --mapali`` — a much
    faster alternative to ``mafft --add`` for large reference alignments.

    The clade's HMM profile aligns the new (orphan) sequences to the profile
    match states, and ``--mapali`` maps the existing reference alignment back in
    WITHOUT re-aligning it, producing one consistent MSA at a cost close to
    aligning the orphans alone (benchmark 16.mafft_vs_hmmalign: ~1.7s vs ~195s
    for mafft --add on a 13k-sequence clade). The afa output is normalised to the
    project convention (uppercase, ``.`` insert-gaps -> ``-``) so it stays a plain
    aligned FASTA usable by SATIVA.

    NOTE: ``--mapali`` expects the reference alignment the HMM was built from.
    A full pipeline run rebuilds the MSAs from scratch, so this holds; the caller
    falls back to ``mafft --add`` if hmmalign fails, which keeps re-runs on an
    already-enriched MSA safe.

    Returns True on success, False otherwise (caller may fall back to MAFFT).
    """
    if not os.path.exists(existing_msa_path):
        print(f"Error: existing MSA file not found: {existing_msa_path}")
        return False
    if not profile_path or not os.path.exists(profile_path) or os.path.getsize(profile_path) == 0:
        print(f"Error: HMM profile missing/empty: {profile_path}")
        return False

    # Write orphan sequences (raw, unaligned) with the standard header format.
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as orph:
        for _, row in new_sequences_df.iterrows():
            header = utils.create_fasta_header(row, tax_columns, id_col)
            orph.write(f"{header}\n{row[seq_col]}\n")
        orph_path = orph.name

    raw_out = existing_msa_path + ".hmmalign.tmp"
    clean_out = output_msa_path + ".clean.tmp"
    try:
        command = [
            'hmmalign', '--dna', '--outformat', 'afa',
            '--mapali', existing_msa_path,   # map reference back in (not re-aligned)
            '-o', raw_out,
            profile_path, orph_path
        ]
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Normalise afa -> plain aligned FASTA (uppercase; '.' insert-gaps -> '-').
        with open(raw_out) as fin, open(clean_out, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    fout.write(line)
                else:
                    fout.write(line.rstrip('\n').upper().replace('.', '-') + '\n')
        os.replace(clean_out, output_msa_path)
        return True

    except subprocess.CalledProcessError as e:
        print(f"hmmalign --mapali failed for {os.path.basename(existing_msa_path)}: {e.stderr}")
        return False
    except Exception as e:
        print(f"Unexpected error during hmmalign integration: {str(e)}")
        return False
    finally:
        for p in [orph_path, raw_out, clean_out]:
            try:
                if p and os.path.exists(p):
                    os.remove(p)
            except Exception:
                pass
