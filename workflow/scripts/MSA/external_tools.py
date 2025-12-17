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
              id_col: str = "seq_id") -> Optional[str]:
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

    Returns:
        Path to the MAFFT alignment file.
    """
    output_dir = Path(output_folder) / barcode / f"alignments_{level}"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = f"{output_dir}/{clade_name}_msa.fasta"
    # print(f"DEBUG3: {output_file}")
    
    # Write temporary input file for MAFFT
    with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta", delete=False) as infile:
        try:
            # Stream sequences into the temporary input
            for _, row in sequences_df.iterrows():
                header = utils.create_fasta_header(row, tax_columns, id_col)
                infile.write(f"{header}\n{row[seq_col]}\n")
            infile_path = infile.name
        except Exception as e:
            print(f"Error while writing MAFFT input file: {str(e)}")
            return None

    # Execute MAFFT
    command = ["mafft", "--auto", "--quiet", "--thread", "-1", infile_path]
    try:
        # Rediriger la sortie directement vers le fichier de sortie
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

def run_hmmscan(sequences_df: pd.DataFrame,
                hmm_db_path: str,
                seq_col: str,
                max_thread_workers: int,
                tax_columns: List[str],
                id_col: str = "seq_id",
                max_evalue: float = 0.001,
                dom_e_value: float = 0.001) -> pd.DataFrame:
    """
    Run hmmscan against a pressed HMM database for multiple sequences.

    Args:
        sequences_df: DataFrame containing the sequences.
        hmm_db_path: path to the pressed HMM database (hmmpress outputs present).
        seq_col: sequence column name.
        max_thread_workers: number of threads passed to hmmscan.
        tax_columns: taxonomy columns to render FASTA headers.
        id_col: identifier column.
        max_evalue: maximum E-value threshold.
        dom_e_value: domain-level E-value (currently unused, placeholder).

    Returns:
        DataFrame with columns [id_col, 'assigned_msa', 'score', 'e_value'].
    """
    # Ensure hmmpress artifacts exist; otherwise run hmmpress
    index_files = ['.h3f', '.h3i', '.h3m', '.h3p']
    needs_press = any(not os.path.exists(hmm_db_path + ext) for ext in index_files)

    if needs_press:
        try:
            print(f"Running hmmpress on {hmm_db_path}...")
            subprocess.run(
                ['hmmpress', '-f', hmm_db_path],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"hmmpress failed: {e.stderr}")
            return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score'])

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
        # Execute hmmscan
        command = [
            'hmmscan',
            '--tblout', tbl_path,
            '--cpu', str(max_thread_workers),
            '--noali',
            # '-E', str(max_evalue),        # Uncomment to enforce an E-value threshold
            # '--domE', str(dom_e_value),   # Domain-level threshold (currently unused)
            hmm_db_path,
            fasta_path
        ]
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Parse hmmscan output (note: hmmscan swaps profile/sequence columns vs hmmsearch)
        hits = []
        with open(tbl_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 18:
                    continue
                # hmmscan: query (parts[2]) = seq_id, target (parts[0]) = profile_name
                seq_id = parts[2].split('|')[0]
                hits.append({
                    id_col: seq_id,  # UUID string
                    'assigned_msa': parts[0],  # Nom du profil HMM
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
        print(f"hmmscan failed: {e.stderr}")
        return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score'])
    finally:
        # Cleanup temporary files
        for path in [fasta_path, tbl_path]:
            try:
                if path and os.path.exists(path):
                    os.remove(path)
            except:
                pass


def create_hmm_database(hmm_dir: str, combined_hmm_path: str) -> bool:
    """
    Concatenate all .hmm files from a directory and run hmmpress on the result.

    Args:
        hmm_dir (str): Directory containing .hmm files (recursively searched).
        combined_hmm_path (str): Destination path for the combined HMM.

    Returns:
        bool: True on success, False otherwise.
    """
    try:
        # Recursively gather all .hmm files
        hmm_files = glob.glob(os.path.join(hmm_dir, '**', '*.hmm'), recursive=True)
        if not hmm_files:
            print(f"Warning: no .hmm files found in {hmm_dir}")
            return False

        # Remove previous hmmpress artifacts
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            if os.path.exists(combined_hmm_path + ext):
                os.remove(combined_hmm_path + ext)

        with open(combined_hmm_path, 'w') as outfile:
            for hmm_file in hmm_files:
                with open(hmm_file, 'r') as infile:
                    outfile.write(infile.read())

        # Run hmmpress
        subprocess.run(
            ['hmmpress', combined_hmm_path],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        return True

    except subprocess.CalledProcessError as e:
        print(f"hmmpress failed: {e.stderr}")
        return False
    except Exception as e:
        print(f"Unexpected error while creating HMM database: {str(e)}")
        return False

def combine_hmm_databases(hmm_paths: List[str], output_path: str) -> bool:
    """
    Merge multiple HMM databases into one.

    Args:
        hmm_paths: List of paths to existing HMM databases.
        output_path: Destination for the combined database.

    Returns:
        True on success, False otherwise.
    """
    try:
        # Remove previous hmmpress artifacts
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            if os.path.exists(output_path + ext):
                os.remove(output_path + ext)

        # Concatenate HMM files
        with open(output_path, 'w') as outfile:
            for hmm_path in hmm_paths:
                if os.path.exists(hmm_path):
                    with open(hmm_path, 'r') as infile:
                        outfile.write(infile.read())
                else:
                    print(f"Warning: missing HMM file {hmm_path}")

        # Rebuild indexes with hmmpress
        subprocess.run(
            ['hmmpress', '-f', output_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        return True
    except subprocess.CalledProcessError as e:
        print(f"hmmpress failed: {e.stderr}")
        return False
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
