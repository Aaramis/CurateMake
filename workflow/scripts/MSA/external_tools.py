# external_tools.py
import os
import glob
import subprocess
import tempfile
import pandas as pd
import os
import utils
from typing import Optional, List
from pathlib import Path
    

def run_mafft(sequences_df: pd.DataFrame,
              barcode: str,
              clade_name: str,
              level: str,
              output_folder: str,
              seq_col: str,
              tax_columns: List[str],
              id_col: str = "seq_id",
              deterministic: bool = False) -> Optional[str]:
    """
    Run MAFFT on a set of sequences, save the result, and return the FASTA path.

    Args:
        sequences_df: DataFrame containing the sequences to align.
        barcode: Barcode label.
        clade_name: Clade name.
        level: Taxonomic level (family, genus, etc.).
        output_folder: Root output directory.
        seq_col: Column holding sequences.
        tax_columns: Taxonomy columns for FASTA headers.
        id_col: Identifier column.
        deterministic: If True, use a fully reproducible MAFFT invocation
            (FFT-NS-2, single thread, input order preserved) and sort sequences
            deterministically before writing. If False (default), use the
            production behaviour (``--auto --thread -1``).

    Returns:
        Path to the MAFFT alignment file.
    """
    output_dir = Path(output_folder) / barcode / f"alignments_{level}"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = f"{output_dir}/{clade_name}_msa.fasta"

    # Write temporary input file for MAFFT
    with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta", delete=False) as infile:
        try:
            # Deterministic mode: sort sequences by stable content keys (not UUID)
            # so the MAFFT input order is reproducible across runs.
            if deterministic:
                sort_cols = [c for c in ['source_id', 'original_id', 'sequence']
                             if c in sequences_df.columns]
                rows_iter = (sequences_df.sort_values(sort_cols, kind="mergesort")
                             if sort_cols else sequences_df).iterrows()
            else:
                rows_iter = sequences_df.iterrows()
            # Stream sequences into the temporary input
            for _, row in rows_iter:
                header = utils.create_fasta_header(row, tax_columns, id_col)
                infile.write(f"{header}\n{row[seq_col]}\n")
            infile_path = infile.name
        except Exception as e:
            print(f"Error while writing MAFFT input file: {str(e)}")
            return None

    # Execute MAFFT
    if deterministic:
        # FFT-NS-2 (--retree 2 --maxiterate 0), single-threaded, input order preserved:
        # fully reproducible alignment.
        command = ["mafft", "--retree", "2", "--maxiterate", "0", "--quiet",
                   "--thread", "1", "--threadit", "0", "--inputorder", infile_path]
    else:
        command = ["mafft", "--auto", "--quiet", "--thread", "-1", infile_path]
    try:
        # Redirect output directly to the output file
        with open(output_file, 'w') as outfile:
            result = subprocess.run(
                command, 
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        # Cleanup temporary file
        os.unlink(infile_path)
        
        return output_file
        
    except subprocess.CalledProcessError as e:
        print(f"MAFFT failed: {e.stderr}")
        # Cleanup on failure
        if os.path.exists(infile_path):
            os.unlink(infile_path)
        return None

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


def create_hmm_database(hmm_dir: str, combined_hmm_path: str) -> bool:
    """
    Concatenate all .hmm files from a directory into a single combined HMM.

    No hmmpress: the combined HMM is consumed by hmmsearch (run_hmmsearch),
    which reads the plain-text profile file directly. Any stale hmmpress index
    (.h3{fimp}) next to the output is purged to avoid confusion / wasted disk.

    Args:
        hmm_dir (str): Directory containing .hmm files (recursively searched).
        combined_hmm_path (str): Destination path for the combined HMM.

    Returns:
        bool: True on success, False otherwise.
    """
    try:
        # Recursively gather all .hmm files, excluding empty marker files
        all_hmm_files = sorted(glob.glob(os.path.join(hmm_dir, '**', '*.hmm'), recursive=True))
        hmm_files = [f for f in all_hmm_files if os.path.getsize(f) > 0]

        skipped = len(all_hmm_files) - len(hmm_files)
        if skipped > 0:
            print(f"Info: skipping {skipped} empty HMM file(s)")

        if not hmm_files:
            print(f"Warning: no valid .hmm files found in {hmm_dir}")
            return False

        # Purge any stale hmmpress index (no longer produced nor read)
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            if os.path.exists(combined_hmm_path + ext):
                os.remove(combined_hmm_path + ext)

        with open(combined_hmm_path, 'w') as outfile:
            for hmm_file in hmm_files:
                with open(hmm_file, 'r') as infile:
                    outfile.write(infile.read())

        return True

    except Exception as e:
        print(f"Unexpected error while creating HMM database: {str(e)}")
        return False

def combine_hmm_databases(hmm_paths: List[str], output_path: str) -> bool:
    """
    Merge multiple HMM databases into one.

    No hmmpress: the merged HMM is consumed by hmmsearch (run_hmmsearch),
    which reads the plain-text profile file directly. Any stale hmmpress index
    (.h3{fimp}) next to the output is purged.

    Args:
        hmm_paths: List of paths to existing HMM databases.
        output_path: Destination for the combined database.

    Returns:
        True on success, False otherwise.
    """
    try:
        # Purge any stale hmmpress index (no longer produced nor read)
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            if os.path.exists(output_path + ext):
                os.remove(output_path + ext)

        # Concatenate HMM files, skipping empty marker files
        skipped = 0
        with open(output_path, 'w') as outfile:
            for hmm_path in sorted(hmm_paths):
                if not os.path.exists(hmm_path):
                    print(f"Warning: missing HMM file {hmm_path}")
                    continue
                if os.path.getsize(hmm_path) == 0:
                    skipped += 1
                    continue
                with open(hmm_path, 'r') as infile:
                    outfile.write(infile.read())

        if skipped > 0:
            print(f"Info: skipping {skipped} empty HMM file(s)")

        return True
    except Exception as e:
        print(f"Unexpected error while combining HMM databases: {str(e)}")
        return False


def run_hmmbuild(msa_path: str, hmm_path: str, barcode: str = None) -> bool:
    """
    Execute hmmbuild to generate a profile from a single MSA.

    Args:
        msa_path (str): path to the MSA FASTA.
        hmm_path (str): destination HMM file.
        barcode (str, optional): optional barcode prefix for the profile name.

    Returns:
        bool: True on success, False otherwise.
    """
    try:
        # Extract clade name from filename
        clade_name = Path(msa_path).stem.replace('_msa', '')

        # Prefix with barcode if provided
        if barcode:
            profile_name = f"{barcode}_{clade_name}"
        else:
            profile_name = clade_name

        subprocess.run([
            'hmmbuild',
            '--dna',
            '-n', profile_name,  # Nom explicite du profile
            hmm_path,
            msa_path
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"hmmbuild failed for {msa_path}: {e.stderr}")
        return False
    except Exception as e:
        print(f"Unexpected error during hmmbuild for {msa_path}: {str(e)}")
        return False



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
    A full pipeline run rebuilds the trimmed MSAs from scratch, so this holds;
    the caller (integrate_clade_task) falls back to ``mafft --add`` if hmmalign
    fails, which keeps re-runs on an already-enriched MSA safe.

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


def run_trimal_and_identify_removed(input_msa_path: str, seq_overlap_threshold: int, seq_resoverlap_threshold: int) -> set:
    """
    Run trimal on an MSA, overwrite the original file with the cleaned result,
    and return the IDs of sequences that were removed.
    """
    if not os.path.exists(input_msa_path) or os.path.getsize(input_msa_path) == 0:
        return set()

    output_dir = os.path.dirname(input_msa_path)
    trimmed_msa_path = "" # Initialisation

    try:
        # 1. Read original sequence IDs before trimming
        original_ids = utils.get_fasta_ids(input_msa_path)
        if not original_ids:
            return set()

        # 2. Create a temporary output file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", dir=output_dir) as tmp_out:
            trimmed_msa_path = tmp_out.name

        # 3. Run trimal with the configured thresholds
        command = [
            'trimal',
            '-in', input_msa_path,
            '-out', trimmed_msa_path,
            '-seqoverlap', str(seq_overlap_threshold),
            '-resoverlap', str(seq_resoverlap_threshold),
            '-noallgaps'  # Drop columns containing only gaps
        ]
        subprocess.run(command, capture_output=True, text=True, check=True)

        # 4. Read IDs from the cleaned file
        trimmed_ids = utils.get_fasta_ids(trimmed_msa_path)

        # 5. Identify removed sequences and replace the original file
        removed_ids = set(original_ids) - set(trimmed_ids)
        os.replace(trimmed_msa_path, input_msa_path)

        return removed_ids

    except subprocess.CalledProcessError as e:
        print(f"trimal failed on {os.path.basename(input_msa_path)}: {e.stderr}")
        if os.path.exists(trimmed_msa_path):
            os.remove(trimmed_msa_path)
        return set()
    except Exception as e:
        print(f"Unexpected error during Trimal cleanup: {e}")
        if os.path.exists(trimmed_msa_path):
            os.remove(trimmed_msa_path)
        return set()
