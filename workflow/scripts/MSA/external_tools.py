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
    Exécute MAFFT sur un ensemble de séquences, sauvegarde le résultat dans le dossier de sortie
    et retourne l'alignement au format FASTA.
    
    Args:
        sequences_df: DataFrame contenant les séquences à aligner
        barcode: Code-barres des séquences
        clade_name: Nom du clade
        level: Niveau taxonomique (genre, famille, etc.)
        output_folder: Dossier de sortie
        seq_col: Nom de la colonne de séquences
        tax_columns: Colonnes taxonomiques
        id_col: Nom de la colonne d'identifiant

        
    Returns:
        str: Path du fichier de sortie
    """
    # Créer le répertoire de sortie

    # output_dir = f"{output_folder}/{barcode}/alignments_{level}"
    output_dir = Path(output_folder) / barcode / f"alignments_{level}"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = f"{output_dir}/{clade_name}_msa.fasta"
    # print(f"DEBUG3: {output_file}")
    
    # Créer un fichier d'entrée temporaire
    with tempfile.NamedTemporaryFile(mode='w+', suffix=".fasta", delete=False) as infile:
        try:
            # Écrire les séquences dans le fichier temporaire
            for _, row in sequences_df.iterrows():
                header = utils.create_fasta_header(row, tax_columns, id_col)
                infile.write(f"{header}\n{row[seq_col]}\n")
            infile_path = infile.name
        except Exception as e:
            print(f"Erreur lors de l'écriture du fichier d'entrée: {str(e)}")
            return None
    
    # Exécuter MAFFT
    #  mafft --genafpair --maxiterate 1000 input_file > output_file   Futur version 
    # https://www.drive5.com/muscle5/ -> Muscle > mafft
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
        
        # Supprimer le fichier temporaire
        os.unlink(infile_path)
        
        return output_file
        
    except subprocess.CalledProcessError as e:
        print(f"Erreur lors de l'exécution de MAFFT: {e.stderr}")
        # Supprimer le fichier temporaire en cas d'erreur
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
    Exécute hmmscan sur un DataFrame de séquences contre une base HMM.

    hmmscan est l'outil approprié pour scanner plusieurs séquences contre une base de profils HMM.
    Il est optimisé pour ce cas d'usage (plus rapide que hmmsearch pour ce scénario).

    Args:
        sequences_df: DataFrame contenant les colonnes ID_COL et SEQUENCE_COL
        hmm_db_path: Chemin vers le fichier HMM pressé (avec hmmpress)
        seq_col: Nom de la colonne de séquences
        max_thread_workers: Nombre max de threads par worker
        tax_columns: Colonnes taxonomiques
        id_col: Nom de la colonne d'identifiant
        max_evalue: Seuil E-value maximum pour conserver un hit (défaut: 0.001)
        dom_e_value: Seuil E-value pour les domaines (défaut: 0.001)

    Returns:
        DataFrame avec les colonnes [id_col, 'assigned_msa', 'score']
    """
    # Vérifier si la base HMM a été pressée, sinon exécuter hmmpress
    index_files = ['.h3f', '.h3i', '.h3m', '.h3p']
    needs_press = any(not os.path.exists(hmm_db_path + ext) for ext in index_files)

    if needs_press:
        try:
            print(f"Exécution de hmmpress sur {hmm_db_path}...")
            subprocess.run(
                ['hmmpress', '-f', hmm_db_path],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Erreur lors de l'exécution de hmmpress: {e.stderr}")
            return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score'])

    # Créer les fichiers temporaires
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as fasta_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.tbl', delete=False) as tbl_file:

        # Écrire le fichier FASTA
        for _, row in sequences_df.iterrows():
            header = utils.create_fasta_header(row, tax_columns, id_col)
            fasta_file.write(f"{header}\n{row[seq_col]}\n")
        fasta_path = fasta_file.name
        tbl_path = tbl_file.name

    try:
        # Exécuter hmmscan (optimisé pour scanner plusieurs séquences)
        command = [
            'hmmscan',
            '--tblout', tbl_path,
            '--cpu', str(max_thread_workers),
            '--noali',
            # '-E', str(max_evalue),        # Seuil E-value pour accélération
            # '--domE', str(dom_e_value), # Seuil E-value pour domaines
            hmm_db_path,
            fasta_path
        ]
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Parser les résultats
        # IMPORTANT: hmmscan inverse les colonnes par rapport à hmmsearch (profil/seq).
        hits = []
        with open(tbl_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 18:
                    continue
                # Pour hmmscan: query (parts[2]) = seq_id, target (parts[0]) = profile_name
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
        print(f"Erreur hmmscan: {e.stderr}")
        return pd.DataFrame(columns=[id_col, 'assigned_msa', 'score'])
    finally:
        # Nettoyer les fichiers temporaires
        for path in [fasta_path, tbl_path]:
            try:
                if path and os.path.exists(path):
                    os.remove(path)
            except:
                pass


def create_hmm_database(hmm_dir: str, combined_hmm_path: str) -> bool:
    """
    Concatène tous les fichiers .hmm d'un dossier et crée une base de données avec hmmpress.

    Args:
        hmm_dir (str): Chemin du dossier contenant les fichiers .hmm
        combined_hmm_path (str): Chemin du fichier .hmm combiné à créer

    Returns:
        bool: True si succès, False sinon
    """
    try:
        # Concaténer les fichiers .hmm (recherche récursive pour supporter des sous-répertoires par niveau)
        hmm_files = glob.glob(os.path.join(hmm_dir, '**', '*.hmm'), recursive=True)
        if not hmm_files:
            print(f"Avertissement: Aucun fichier .hmm trouvé dans {hmm_dir}")
            return False

        # Supprimer les anciens fichiers d'index s'ils existent
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            if os.path.exists(combined_hmm_path + ext):
                os.remove(combined_hmm_path + ext)

        with open(combined_hmm_path, 'w') as outfile:
            for hmm_file in hmm_files:
                with open(hmm_file, 'r') as infile:
                    outfile.write(infile.read())

        # Lancer hmmpress
        subprocess.run(
            ['hmmpress', combined_hmm_path],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        return True

    except subprocess.CalledProcessError as e:
        print(f"Erreur lors de l'exécution de hmmpress: {e.stderr}")
        return False
    except Exception as e:
        print(f"Erreur inattendue lors de la création de la base HMM: {str(e)}")
        return False

def combine_hmm_databases(hmm_paths: List[str], output_path: str) -> bool:
    """
    Combine plusieurs bases HMM en une seule.

    Args:
        hmm_paths: Liste des chemins vers les bases HMM à combiner
        output_path: Chemin de sortie de la base combinée

    Returns:
        True si succès, False sinon
    """
    try:
        # Supprimer les anciens index
        for ext in ['.h3f', '.h3i', '.h3m', '.h3p']:
            if os.path.exists(output_path + ext):
                os.remove(output_path + ext)

        # Concaténer les fichiers HMM
        with open(output_path, 'w') as outfile:
            for hmm_path in hmm_paths:
                if os.path.exists(hmm_path):
                    with open(hmm_path, 'r') as infile:
                        outfile.write(infile.read())
                else:
                    print(f"Avertissement: Fichier HMM introuvable: {hmm_path}")

        # Indexer avec hmmpress
        subprocess.run(
            ['hmmpress', '-f', output_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        return True
    except subprocess.CalledProcessError as e:
        print(f"Erreur lors de l'exécution de hmmpress: {e.stderr}")
        return False
    except Exception as e:
        print(f"Erreur inattendue lors de la combinaison des bases HMM: {str(e)}")
        return False


def run_hmmbuild(msa_path: str, hmm_path: str, barcode: str = None) -> bool:
    """
    Exécute hmmbuild pour générer un fichier HMM à partir d'un alignement MSA avec un nom unique.

    Args:
        msa_path (str): Chemin du fichier d'alignement MSA (FASTA)
        hmm_path (str): Chemin de sortie du fichier HMM
        barcode (str, optional): Barcode pour préfixer le nom du profile (éviter les conflits)

    Returns:
        bool: True si succès, False sinon
    """
    try:
        # Extraire le nom du clade depuis le path
        clade_name = Path(msa_path).stem.replace('_msa', '')

        # Préfixer avec le barcode si fourni
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
        print(f"Erreur lors de l'exécution de hmmbuild pour {msa_path}: {e.stderr}")
        return False
    except Exception as e:
        print(f"Erreur inattendue lors de hmmbuild pour {msa_path}: {str(e)}")
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
        seq_col: Nom de la colonne de séquences
        tax_columns: Colonnes taxonomiques
        id_col: Nom de la colonne d'identifiant
        
    Returns:
        True if successful, False otherwise.
    """
    if not os.path.exists(existing_msa_path):
        print(f"Erreur : Le fichier MSA existant n'a pas été trouvé : {existing_msa_path}")
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
            print(f"Erreur lors de l'exécution de mafft --add pour {os.path.basename(existing_msa_path)}.")
            print(f"Stderr: {e.stderr}")
            return False
        except Exception as e:
            print(f"Une erreur inattendue est survenue lors de l'alignement de profil : {str(e)}")
            return False


def run_trimal_and_identify_removed(input_msa_path: str, seq_overlap_threshold: int, seq_resoverlap_threshold: int) -> set:
    """
    Exécute trimal sur un MSA, écrase le fichier d'entrée avec le résultat nettoyé,
    et retourne les IDs des séquences qui ont été supprimées.
    """
    if not os.path.exists(input_msa_path) or os.path.getsize(input_msa_path) == 0:
        return set()

    output_dir = os.path.dirname(input_msa_path)
    trimmed_msa_path = "" # Initialisation

    try:
        # 1. Lire les IDs de séquences originales avant le nettoyage
        original_ids = utils.get_fasta_ids(input_msa_path)
        if not original_ids:
            return set()

        # 2. Créer un fichier de sortie temporaire
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", dir=output_dir) as tmp_out:
            trimmed_msa_path = tmp_out.name

        # 3. Construire et exécuter la commande trimal
        command = [
            'trimal',
            '-in', input_msa_path,
            '-out', trimmed_msa_path,
            '-seqoverlap', str(seq_overlap_threshold),
            '-resoverlap', str(seq_resoverlap_threshold),
            '-noallgaps' # Retire les colonnes avec uniquement des gaps
        ]
        subprocess.run(command, capture_output=True, text=True, check=True)

        # 4. Lire les IDs du fichier nettoyé
        trimmed_ids = utils.get_fasta_ids(trimmed_msa_path)

        # 5. Identifier les séquences supprimées et mettre à jour le fichier MSA
        removed_ids = set(original_ids) - set(trimmed_ids)
        os.replace(trimmed_msa_path, input_msa_path)

        return removed_ids

    except subprocess.CalledProcessError as e:
        print(f"Erreur Trimal sur {os.path.basename(input_msa_path)}: {e.stderr}")
        if os.path.exists(trimmed_msa_path):
            os.remove(trimmed_msa_path)
        return set()
    except Exception as e:
        print(f"Erreur inattendue durant le nettoyage Trimal: {e}")
        if os.path.exists(trimmed_msa_path):
            os.remove(trimmed_msa_path)
        return set()
