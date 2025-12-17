# processor.py
import pandas as pd
import os
import logging
import json
import time
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import concurrent.futures
import threading
import external_tools
import utils
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Set


def get_next_tax_level(current_level, tax_columns):
    try:
        current_index = tax_columns.index(current_level)
        if current_index < len(tax_columns) - 1:
            return tax_columns[current_index + 1]
    except (ValueError, IndexError):
        pass
    return None


def generate_reassignment_statistics(assignments_df: pd.DataFrame, id_col: str = 'seq_id') -> dict:
    """
    Génère les statistiques de réaffectation inter-barcodes.

    Args:
        assignments_df: DataFrame avec colonnes 'source_barcode', 'assigned_msa', 'seq_id'
        id_col: Nom de la colonne d'identifiant

    Returns:
        Dictionnaire avec les statistiques détaillées
    """
    if assignments_df.empty:
        return {}

    stats = {}

    # Créer une copie pour éviter les modifications
    df = assignments_df.copy()

    # Extraire le barcode cible depuis le profile name
    # Exemple: "ITS2_Agaricus_msa" → "ITS2"
    df['target_barcode'] = df['assigned_msa'].str.extract(r'^(ITS\d?)_')[0]

    # Gérer les cas où le barcode cible ne peut pas être extrait (ancien format)
    df['target_barcode'].fillna('Unknown', inplace=True)

    # Calculer les pourcentages
    for source in ['ITS', 'ITS1', 'ITS2']:
        source_total = (df['source_barcode'] == source).sum()

        for target in ['ITS', 'ITS1', 'ITS2']:
            count = ((df['source_barcode'] == source) &
                    (df['target_barcode'] == target)).sum()

            percentage = (count / source_total * 100) if source_total > 0 else 0

            stats[f'{source}_to_{target}'] = {
                'count': int(count),
                'percentage': float(percentage)
            }

    return stats


class HierarchicalAligner:
    """
    Gère le processus d'alignement hiérarchique des séquences.
    """
    def __init__(self,
                 data_df: pd.DataFrame,
                 min_seq_threshold: int = 20,
                 max_seq_threshold: int = 1000,
                 large_task_threshold: int = 2000,
                 min_domain: str = "species",
                 trimal_seq_overlap: float = 75,
                 trimal_res_overlap: float = 0.5,
                 domain: str = "genus",
                 barcodes: List[str] = None,
                 hierarchy_levels: List[str] = None,
                 tax_columns: List[str] = None,
                 id_col: str = "seq_id",
                 seq_col: str = "sequences",
                 barcode_col: str = "barcode",
                 output_folder: str = "output",
                 max_barcode_workers: int = 3,
                 max_thread_workers: int = 4,
                 classification_max_workers: int = 8,
                 classification_threads_per_worker: int = 10,
                 classification_batch_size: int = 150000,
                logger: Optional[logging.Logger] = None,
                logger_folder: str = "logs/MSA",
                debug_mode: bool = False,
                test_limit: Optional[int] = None):
        """
        Initialise l'aligneur hiérarchique avec tous les paramètres nécessaires.

        Args:
            data_df: DataFrame contenant les données de séquences
            min_seq_threshold: Seuil minimum de séquences par groupe
            trimal_seq_overlap: Paramètre de chevauchement de séquence pour trimal
            trimal_res_overlap: Paramètre de chevauchement de résidu pour trimal
            domain: Niveau taxonomique principal pour le groupement
            barcodes: Liste des barcodes à traiter
            hierarchy_levels: Niveaux hiérarchiques taxonomiques
            tax_columns: Colonnes taxonomiques
            id_col: Nom de la colonne d'identifiant
            seq_col: Nom de la colonne de séquences
            barcode_col: Nom de la colonne de barcode
            output_folder: Dossier de sortie
            max_barcode_workers: Nombre max de workers pour les barcodes
            max_thread_workers: Nombre max de threads par worker
            classification_max_workers: Nombre de processus hmmscan en parallèle
            classification_threads_per_worker: Threads par processus hmmscan
            classification_batch_size: Séquences par batch pour la classification
            logger: Un objet logger pré-configuré. Si None, un nouveau sera créé.
            logger_folder: Dossier de sortie pour les logs
            debug_mode: Mode debug pour logs détaillés
            test_limit: Limitation du nombre de séquences pour les tests
        """
        # Appliquer la limitation si test_limit est défini
        if test_limit is not None and test_limit > 0:
            self.data_df = data_df.head(test_limit).copy()
            print(f"MODE TEST: Limitation à {test_limit} séquences (sur {len(data_df)} totales)")
        else:
            self.data_df = data_df.copy()

        self.min_seq_threshold = min_seq_threshold
        self.max_seq_threshold = max_seq_threshold
        self.large_task_threshold = large_task_threshold
        self.trimal_seq_overlap = trimal_seq_overlap
        self.trimal_res_overlap = trimal_res_overlap
        self.domain = domain
        self.min_domain = min_domain
        self.barcodes = barcodes or ["ITS", "ITS1", "ITS2"]
        self.hierarchy_levels = hierarchy_levels or ["genus", "family", "order", "class", "phylum"]
        self.tax_columns = tax_columns or ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        self.id_col = id_col
        self.seq_col = seq_col
        self.barcode_col = barcode_col
        self.output_folder = output_folder
        self.max_barcode_workers = max_barcode_workers
        self.max_thread_workers = max_thread_workers
        self.logger_folder = logger_folder
        self.debug_mode = debug_mode

        # Configuration de la classification par lots avec hmmscan (depuis config.yaml)
        self.classification_max_workers = classification_max_workers
        self.classification_threads_per_worker = classification_threads_per_worker
        self.classification_batch_size = classification_batch_size

        # Créer un verrou pour protéger l'accès aux sets partagés
        self.lock = threading.Lock()

        # Créer le dossier de sortie
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
    
        if logger:
            self.logger = logger
        else:
            self._setup_logging()

    def _setup_logging(self):
        """Configure le système de logging par défaut si aucun logger n'est fourni."""
        log_dir = os.path.join(self.output_folder, "logs")
        Path(log_dir).mkdir(parents=True, exist_ok=True)
        
        # Logger principal
        self.logger = logging.getLogger(f"HierarchicalAligner")
        self.logger.setLevel(logging.INFO)
        
        # Évite d'ajouter plusieurs handlers si la méthode est appelée plusieurs fois
        if not self.logger.handlers:
            handler = logging.FileHandler(os.path.join(log_dir, "main.log"))
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def _get_barcode_logger(self, barcode: str) -> logging.Logger:
        """Obtient ou crée un logger spécifique pour un barcode."""
        logger_name = f"HierarchicalAligner.{barcode}"
        logger = logging.getLogger(logger_name)
        
        if not logger.handlers:
            log_dir = os.path.join(self.logger_folder, "logs")
            Path(log_dir).mkdir(parents=True, exist_ok=True)
            
            handler = logging.FileHandler(os.path.join(log_dir, f"{barcode}.log"))
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        
        return logger

    def run(self) -> Dict[str, Any]:
        """
        Lance le workflow MSA avec l'architecture unifiée en 2 phases.

        Phase 1: Préparation des MSA et HMM par barcode
        Phase 2: Classification unifiée inter-barcodes

        Returns:
            Dictionnaire avec les résultats complets
        """
        self.logger.info("Démarrage du processus d'alignement hiérarchique (architecture unifiée).")

        # PHASE 1: Préparation par barcode (inchangé)
        barcode_results = self._phase1_prepare_msa_and_hmm()

        # PHASE 2: Classification unifiée (NOUVEAU)
        global_assignments = self._phase2_unified_classification(barcode_results)

        # Construction du rapport final
        serializable_results = {}
        for barcode, data in barcode_results.items():
            serializable_data = {}
            for key, value in data.items():
                if key == "hmm_db_path" and value is not None:
                    serializable_data[key] = value
                elif key == "msa_paths" and value is not None:
                    serializable_data[key] = f"Number of paths: {len(value)}"
                elif isinstance(value, int):
                    serializable_data[key] = value
            serializable_results[barcode] = serializable_data

        self.logger.info("Workflow MSA terminé avec succès.")

        return {
            "barcode_results": serializable_results,
            "global_assignments": global_assignments,
            "total_barcodes_processed": len(barcode_results),
            "phases_completed": 2
        }


    def _split_group_into_tasks(self, group: pd.DataFrame, level: str, 
                                clade_name: str, max_sequences: int) -> List[Dict]:
        """
        Split un groupe trop grand en plusieurs tâches/MSA de taille max_sequences.
        
        Args:
            group: DataFrame du groupe à diviser
            level: Niveau taxonomique
            clade_name: Nom du clade
            max_sequences: Taille maximale par sous-groupe
        
        Returns:
            Liste de tâches pour les sous-groupes
        """
        tasks = []
        num_splits = (len(group) + max_sequences - 1) // max_sequences  # Arrondi supérieur
        
        for i in range(num_splits):
            start_idx = i * max_sequences
            end_idx = min((i + 1) * max_sequences, len(group))
            sub_group = group.iloc[start_idx:end_idx].copy()
            
            # Créer un nom unique pour chaque sous-MSA
            sub_clade_name = f"{clade_name}_part{i+1}"
            
            tasks.append({
                'level': level,
                'clade_name': sub_clade_name,
                'original_clade_name': clade_name,  # Garder trace du nom original
                'part_number': i + 1,
                'total_parts': num_splits,
                'dataframe': sub_group,
                'sequence_count': len(sub_group),
                'is_split': True  # Marquer que c'est un split artificiel
            })
        
        self.logger.info(f"  - {clade_name} divisé en {num_splits} parties de ~{max_sequences} séquences")
        
        return tasks

    def _create_tasks_recursive(self, df, level, barcode, max_sequences, processed_indices):
        """Méthode de classe maintenant"""
        tasks = []
        grouped = df.groupby(level)

        for clade_name, group in grouped:
            if isinstance(clade_name, str):
                clade_name = clade_name.replace(" ", "_")

            # Vérifier le seuil minimum de séquences pour tous les niveaux
            if len(group) < self.min_seq_threshold:
                continue

            # Ignorer les clades marqués comme "incertae sedis"
            if isinstance(clade_name, str) and "incertae sedis" in clade_name.lower():
                continue

            # CAS 1: Le groupe est de taille acceptable
            if len(group) <= max_sequences:
                tasks.append({
                    'level': level,
                    'clade_name': clade_name,
                    'dataframe': group.copy(),
                    'sequence_count': len(group)
                })
                processed_indices.update(group.index)
            
            # CAS 2: Le groupe est trop grand
            else:
                # Si on n'est PAS encore au niveau minimum, on peut descendre
                if level != self.min_domain:
                    next_level = get_next_tax_level(level, self.tax_columns)
                    if next_level is not None and next_level in group.columns:
                        if not group[next_level].dropna().empty:
                            # self.logger.info(f"{barcode} - Subdivision de {clade_name} ({len(group)} séquences) au niveau {next_level}")
                            sub_tasks = self._create_tasks_recursive(group, next_level, barcode, max_sequences, processed_indices)
                            tasks.extend(sub_tasks)
                        else:
                            # Pas de données au niveau suivant, on split artificiellement
                            self.logger.warning(f"{barcode} - Impossible de subdiviser {clade_name} au niveau {next_level} (données manquantes). Split artificiel.")
                            tasks.extend(self._split_group_into_tasks(group, level, clade_name, max_sequences))
                            processed_indices.update(group.index)
                    else:
                        # Pas de niveau suivant, on split artificiellement
                        self.logger.warning(f"{barcode} - Pas de niveau suivant pour {clade_name}. Split artificiel.")
                        tasks.extend(self._split_group_into_tasks(group, level, clade_name, max_sequences))
                        processed_indices.update(group.index)
                
                # Si on EST au niveau minimum (species), on split en plusieurs MSA
                else:
                    # self.logger.info(f"{barcode} - Niveau minimum atteint pour {clade_name} ({len(group)} séquences). Création de plusieurs MSA.")
                    tasks.extend(self._split_group_into_tasks(group, level, clade_name, max_sequences))
                    processed_indices.update(group.index)
        return tasks


    def create_tasks(self, barcode_df: pd.DataFrame, barcode: str):
        """Crée les tâches de traitement pour un barcode donné.
        Args:
            barcode_df: DataFrame contenant les données du barcode à traiter
        Returns:
            tuple: (tasks, unprocessed_sequences_df) où:
            - tasks: Liste des tâches à exécuter
            - unprocessed_sequences_df: DataFrame des séquences non traitées
        """
        # Compter les occurrences de chaque domaine taxonomique
        self.logger.info(f"{barcode} - nombre de séquence total {len(barcode_df)}")
        domain_counts = barcode_df[self.domain].dropna().value_counts()
        valid_domains = domain_counts[domain_counts >= self.min_seq_threshold].index
        
        # Filtrage des données valides
        valid_df = barcode_df[barcode_df[self.domain].isin(valid_domains)]

        # Set pour tracker les indices des séquences traitées
        processed_indices = set()

        # Création des tâches avec subdivision automatique
        all_tasks = self._create_tasks_recursive(valid_df, self.domain, barcode, self.max_seq_threshold, processed_indices)
        
        # Séparer les tâches normales des grandes tâches
        normal_tasks = []
        large_tasks = []
        
        for task in all_tasks:
            if task['sequence_count'] >= self.large_task_threshold:
                large_tasks.append(task)
            else:
                normal_tasks.append(task)


        # Log des tâches créées
        tot_sequence_count = 0
        for i, task in enumerate(all_tasks, 1):
            tot_sequence_count += task['sequence_count']
            self.logger.info(f" {barcode}. {task['clade_name']} ({task['level']}) - {task['sequence_count']} séquences")
        self.logger.info(f" {barcode} - {tot_sequence_count} séquences traitées dans les tasks")
        self.logger.info(f" {barcode} - {len(normal_tasks)} tâches normales, {len(large_tasks)} grandes tâches")

        # Récupération des séquences non traitées : toutes celles qui ne sont pas dans une tâche
        unprocessed_df = barcode_df[~barcode_df.index.isin(processed_indices)].copy()
        self.logger.info(f" {barcode} - {len(unprocessed_df)} séquences non traitées dans les tasks")

        return normal_tasks, large_tasks, unprocessed_df

    def _add_sequences_to_msa(self, assignments_df: pd.DataFrame, unprocessed_df: pd.DataFrame, 
                             barcode: str, msa_output_paths: Dict):
        """
        Ajoute les séquences nouvellement assignées à leurs MSA respectifs.
        """
        if assignments_df.empty:
            return

        logger = self._get_barcode_logger(barcode)
        logger.info(f"{barcode} - Ajout de {len(assignments_df)} séquences assignées aux MSA existants")
        
        # Le nom du MSA assigné est dans la colonne 'assigned_msa' (ex: 'ITS2_Cupressus_msa')
        # On retire le préfixe barcode (ITS, ITS1, ITS2) et le suffixe '_msa'
        assignments_df = assignments_df.copy()
        assignments_df['clade_name'] = (assignments_df['assigned_msa']
                                        .str.replace(r'^(ITS\d?)_', '', regex=True)
                                        .str.replace('_msa', ''))
        
        grouped_assignments = assignments_df.groupby('clade_name')

        for clade_name, group_df in tqdm(grouped_assignments, desc=f"Mise à jour des MSA pour {barcode}"):
            # Récupérer le chemin du MSA original ET son niveau taxonomique
            msa_info = next(((k, v) for k, v in msa_output_paths.items() if k[0] == clade_name), None)

            if not msa_info:
                logger.warning(f"{barcode} - Impossible de trouver le chemin du MSA pour le clade '{clade_name}'. Séquences ignorées.")
                continue

            (clade_key, level), original_msa_path = msa_info

            # FILTRE : Ne traiter que les MSA de niveau "species" à commenter
            # if level != "species":
            # if level == "species":
            #     logger.debug(f"{barcode} - MSA {clade_name} (niveau {level}) ignoré")
            #     continue

            # Obtenir les IDs des séquences à ajouter
            sequence_ids_to_add = group_df[self.id_col].tolist()
            
            # Filtrer le DataFrame non traité pour obtenir les séquences complètes
            sequences_to_add_df = unprocessed_df[unprocessed_df[self.id_col].isin(sequence_ids_to_add)]

            if sequences_to_add_df.empty:
                logger.warning(f"{barcode} - Aucune séquence correspondante trouvée pour le clade '{clade_name}'.")
                continue

            # Lancer mafft --add
            try:
                success = external_tools.run_mafft_add(
                    new_sequences_df=sequences_to_add_df,
                    existing_msa_path=original_msa_path,
                    output_msa_path=original_msa_path,
                    seq_col=self.seq_col,
                    tax_columns=self.tax_columns,
                    max_thread_workers=self.max_thread_workers,
                    id_col=self.id_col
                )

                if not success:
                    logger.error(f"{barcode} - Échec de la mise à jour du MSA pour {clade_name}.")
                # else:
                #     logger.info(f"{barcode} - MSA mis à jour avec succès pour {clade_name}.")
                    
            except Exception as e:
                logger.error(f"{barcode} - Erreur lors de la mise à jour du MSA pour {clade_name}: {e}")



    def _prepare_barcode_data(self, barcode: str, logger) -> Tuple[pd.DataFrame, pd.DataFrame, List[Dict]]:
        """Prépare les données et crée les tâches pour un barcode."""
        barcode_df = self.data_df[self.data_df[self.barcode_col] == barcode]
        
        if barcode_df.empty:
            logger.warning(f" {barcode} - Aucune donnée trouvée")
            return barcode_df, pd.DataFrame(), []
        
        normal_tasks, large_tasks, unprocessed_df = self.create_tasks(barcode_df, barcode)
        all_tasks = normal_tasks + large_tasks
        
        if not unprocessed_df.empty:
            logger.info(f"{barcode} - {len(barcode_df)} séquences")
            logger.info(f"{barcode} - {len(unprocessed_df)} séquences non traitées (pas assez de représentants par {self.domain})")
        logger.info(f"{barcode} - {len(all_tasks)} tâches créées")
        
        return barcode_df, unprocessed_df, all_tasks


    def _generate_msa(self, barcode: str, all_tasks: List[Dict], logger) -> Dict[Tuple[str, str], str]:
        """Génère les alignements multiples (MSA) pour toutes les tâches."""
        msa_output_paths = {}
        
        if not all_tasks:
            return msa_output_paths
        
        with ThreadPoolExecutor(max_workers=self.max_thread_workers) as executor:
            futures = {
                executor.submit(
                    external_tools.run_mafft,
                    task['dataframe'],
                    barcode,
                    task['clade_name'],
                    task['level'],
                    self.output_folder,
                    self.seq_col,
                    self.tax_columns,
                    self.id_col
                ): (task['clade_name'], task['level']) for task in all_tasks
            }
            
            for future in tqdm(
                concurrent.futures.as_completed(futures),
                total=len(all_tasks),
                desc=f"Exécution des alignements pour {barcode}"
            ):
                clade_name, level = futures[future]
                try:
                    output_path = future.result()
                    if output_path is None:
                        logger.warning(f"{barcode} - Aucun résultat pour {clade_name}")
                    else:
                        msa_output_paths[(clade_name, level)] = output_path
                except Exception as e:
                    logger.error(f"{barcode} - Erreur lors de l'alignement pour {clade_name}: {str(e)}")
        
        return msa_output_paths


    def _clean_msa_with_trimal(
        self, barcode: str, barcode_df: pd.DataFrame, 
        barcode_output_dir: str, msa_output_paths: Dict, 
        all_tasks: List[Dict], unprocessed_df: pd.DataFrame, logger
    ) -> pd.DataFrame:
        """Nettoie les MSA avec Trimal et identifie les séquences mal alignées."""
        normal_tasks = [task for task in all_tasks if not task.get('is_large', False)]
        
        if not msa_output_paths or not normal_tasks:
            return unprocessed_df
        
        logger.info(f"{barcode} - Lancement du nettoyage des MSA avec Trimal")
        
        # Filtrer les MSA des tâches normales
        normal_task_clades = {task['clade_name'] for task in normal_tasks}
        normal_msa_paths = {
            (clade_name, level): msa_path 
            for (clade_name, level), msa_path in msa_output_paths.items() 
            if clade_name in normal_task_clades
        }
        
        # Exécuter Trimal en parallèle
        all_removed_ids = self._run_trimal_parallel(barcode, normal_msa_paths, logger)
        
        # Mettre à jour le pool de séquences non traitées
        if all_removed_ids:
            unprocessed_df = self._update_unprocessed_pool(
                barcode, barcode_df, all_removed_ids, unprocessed_df, logger
            )
        
        # Sauvegarder le pool de séquences non traitées
        if not unprocessed_df.empty:
            unprocessed_path = os.path.join(barcode_output_dir, 'unprocessed_after_trimal.csv')
            unprocessed_df.to_csv(unprocessed_path, index=False)
            logger.info(f"{barcode} - Pool de séquences non traitées sauvegardé dans {unprocessed_path}")
        
        # Logger les informations
        self._log_trimal_summary(barcode, all_tasks, msa_output_paths, logger)
        
        return unprocessed_df


    def _run_trimal_parallel(
        self, barcode: str, normal_msa_paths: Dict, logger
    ) -> Set[str]:
        """Exécute Trimal en parallèle sur tous les MSA."""
        all_removed_ids = set()
        
        with ThreadPoolExecutor(max_workers=self.max_thread_workers) as executor:
            futures_to_clade = {
                executor.submit(
                    external_tools.run_trimal_and_identify_removed, 
                    msa_path,
                    self.trimal_seq_overlap,
                    self.trimal_res_overlap
                ): clade_name for (clade_name, _), msa_path in normal_msa_paths.items()
            }

            for future in tqdm(
                concurrent.futures.as_completed(futures_to_clade),
                total=len(futures_to_clade),
                desc=f"Nettoyage des MSA (Trimal) pour {barcode}"
            ):
                clade_name = futures_to_clade[future]
                try:
                    removed_ids = future.result()
                    if removed_ids:
                        all_removed_ids.update(removed_ids)
                except Exception as e:
                    logger.error(f"{barcode} - Erreur lors du nettoyage du MSA pour {clade_name} avec Trimal: {str(e)}")
        
        return all_removed_ids


    def _update_unprocessed_pool(
        self, barcode: str, barcode_df: pd.DataFrame,
        removed_ids: Set[str], unprocessed_df: pd.DataFrame, logger
    ) -> pd.DataFrame:
        """Met à jour le pool de séquences non traitées avec les séquences retirées par Trimal."""
        logger.info(f"{barcode} - Total de {len(removed_ids)} séquences mal alignées à ajouter au pool non traité.")

        # DEBUG: Analyser l'état initial
        initial_unprocessed_count = len(unprocessed_df)
        logger.info(f"{barcode} - État initial: {initial_unprocessed_count} séquences non traitées")

        if self.debug_mode:
            logger.info(f"{barcode} - Colonnes dans unprocessed_df: {list(unprocessed_df.columns)}")
            logger.info(f"{barcode} - Colonnes dans barcode_df: {list(barcode_df.columns)}")
            logger.info(f"{barcode} - Type de {self.id_col} dans unprocessed_df: {unprocessed_df[self.id_col].dtype if len(unprocessed_df) > 0 else 'DF vide'}")
            logger.info(f"{barcode} - Type de {self.id_col} dans barcode_df: {barcode_df[self.id_col].dtype}")

        # Parser les IDs pour extraire uniquement la partie ID (avant le premier |)
        # Les IDs venant de FASTA peuvent contenir des métadonnées: "622932|kingdom=Fungi|..."
        parsed_ids = set()
        for id_val in removed_ids:
            id_str = str(id_val)
            # Extraire la partie avant le premier "|" si présent
            if '|' in id_str:
                id_clean = id_str.split('|')[0]
            else:
                id_clean = id_str
            parsed_ids.add(id_clean)

        logger.info(f"{barcode} - IDs parsés: {len(parsed_ids)} IDs extraits")

        # Convertir les IDs parsés en strings (UUIDs sont toujours des strings)
        # Note: Ancien code supportait integer IDs, maintenant nous utilisons exclusivement UUID strings
        removed_ids_converted = {str(id_val) for id_val in parsed_ids}
        logger.info(f"{barcode} - IDs préparés: {len(removed_ids_converted)} IDs (UUID strings)")

        # DEBUG: Vérifier quelques exemples
        if self.debug_mode and removed_ids_converted:
            example_id = next(iter(removed_ids_converted))
            example_df_id = barcode_df[self.id_col].iloc[0]
            logger.info(f"{barcode} - Example removed ID: {example_id} (type: {type(example_id)})")
            logger.info(f"{barcode} - Example DF ID: {example_df_id} (type: {type(example_df_id)})")
            logger.info(f"{barcode} - ID match test: {example_id in barcode_df[self.id_col].values}")

        # Filtrer pour obtenir les séquences retirées
        removed_sequences_df = barcode_df[barcode_df[self.id_col].isin(removed_ids_converted)].copy()
        logger.info(f"{barcode} - Séquences trouvées correspondant aux IDs retirés: {len(removed_sequences_df)}")

        if len(removed_sequences_df) == 0:
            logger.warning(
                f"{barcode} - ATTENTION: Aucune séquence trouvée malgré {len(removed_ids)} IDs retirés! "
                f"Problème probable de correspondance de types."
            )
            if self.debug_mode:
                # Afficher quelques IDs pour déboguer
                sample_removed = list(removed_ids_converted)[:5]
                sample_df = barcode_df[self.id_col].head(5).tolist()
                logger.info(f"{barcode} - Sample removed IDs: {sample_removed}")
                logger.info(f"{barcode} - Sample DF IDs: {sample_df}")
            return unprocessed_df

        # DEBUG: Vérifier les colonnes avant concat
        if self.debug_mode:
            logger.info(f"{barcode} - Forme de unprocessed_df avant concat: {unprocessed_df.shape}")
            logger.info(f"{barcode} - Forme de removed_sequences_df: {removed_sequences_df.shape}")

        # Concatener les DataFrames
        unprocessed_df = pd.concat([unprocessed_df, removed_sequences_df], ignore_index=True)
        logger.info(f"{barcode} - Après concat: {len(unprocessed_df)} séquences (ajout de {len(removed_sequences_df)})")

        # Supprimer les duplicats
        before_dedup = len(unprocessed_df)
        unprocessed_df.drop_duplicates(subset=[self.id_col], keep='first', inplace=True)
        after_dedup = len(unprocessed_df)

        if before_dedup != after_dedup:
            logger.info(f"{barcode} - {before_dedup - after_dedup} duplicats supprimés")

        final_count = len(unprocessed_df)
        added_count = final_count - initial_unprocessed_count
        logger.info(f"{barcode} - Taille finale du pool: {final_count} ({added_count:+d} séquences)")

        return unprocessed_df

    def _log_trimal_summary(
        self, barcode: str, all_tasks: List[Dict], 
        msa_output_paths: Dict, logger
    ):
        """Log le résumé du traitement Trimal."""
        large_tasks = [task for task in all_tasks if task.get('is_large', False)]
        
        if large_tasks:
            logger.info(f"{barcode} - {len(large_tasks)} grandes tâches non traitées par Trimal:")
            for task in large_tasks:
                logger.info(f"  - {task['clade_name']}: {task['sequence_count']} séquences")
        
        logger.info(f"{barcode} - Résumé des MSA générés:")
        for (clade_name, _), msa_path in msa_output_paths.items():
            try:
                num_seqs, aln_length = utils.analyze_msa_file(msa_path, self.logger)
                # logger.info(f"{barcode} - {clade_name}: {num_seqs} séquences, {aln_length} positions d'alignement")
            except Exception as e:
                logger.warning(f"Impossible d'analyser le MSA pour {clade_name}: {e}")


    def _generate_hmm_database(
        self, barcode: str, barcode_output_dir: str, 
        msa_output_paths: Dict, logger
    ) -> Optional[str]:
        """Génère les HMM et crée la base de données combinée."""
        if not msa_output_paths:
            return None
        
        hmm_dir = os.path.join(barcode_output_dir, 'hmm')
        os.makedirs(hmm_dir, exist_ok=True)
        
        # Générer les HMM individuels
        self._generate_individual_hmm(barcode, msa_output_paths, hmm_dir, logger)
        
        # Créer la base de données combinée
        combined_hmm_path = os.path.join(barcode_output_dir, 'all_profiles.hmm')
        logger.info(f"Création de la base de données HMM pour {barcode}...")
        
        success = external_tools.create_hmm_database(hmm_dir, combined_hmm_path)
        if not success:
            logger.error("Échec de la création de la base de données HMM")
            return None
        
        return combined_hmm_path


    def _generate_individual_hmm(
        self, barcode: str, msa_output_paths: Dict, hmm_dir: str, logger
    ):
        """Génère les fichiers HMM individuels pour chaque MSA avec préfixe barcode."""
        with ThreadPoolExecutor(max_workers=self.max_thread_workers) as executor:
            futures = {}
            for (clade_name, _), msa_path in msa_output_paths.items():
                hmm_path = os.path.join(hmm_dir, f"{clade_name}.hmm")
                # Passer le barcode pour nommage unique des profiles HMM
                future = executor.submit(external_tools.run_hmmbuild, msa_path, hmm_path, barcode)
                futures[future] = clade_name

            for future in tqdm(
                concurrent.futures.as_completed(futures),
                total=len(futures),
                desc=f"Génération des HMM pour {barcode}"
            ):
                clade_name = futures[future]
                try:
                    success = future.result()
                    if not success:
                        logger.warning(f"La génération du HMM a échoué pour {clade_name}")
                except Exception as e:
                    logger.error(f"Erreur lors de la génération du HMM pour {clade_name}: {str(e)}")


    def _create_batches(self, df: pd.DataFrame, batch_size: int) -> List[pd.DataFrame]:
        """Crée des lots de DataFrames pour le traitement parallèle.

        Args:
            df: DataFrame à découper
            batch_size: Nombre de séquences par lot

        Returns:
            Liste de DataFrames (batches)
        """
        return [
            df.iloc[i:i + batch_size].copy()
            for i in range(0, len(df), batch_size)
        ]

    def _init_checkpoint(self, checkpoint_path: str, num_batches: int, barcode: str) -> Dict:
        """Initialise un nouveau fichier checkpoint.

        Args:
            checkpoint_path: Chemin du fichier checkpoint
            num_batches: Nombre total de batches
            barcode: Code-barres en cours de traitement

        Returns:
            Dictionnaire du checkpoint
        """
        checkpoint = {
            "barcode": barcode,
            "total_batches": num_batches,
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "last_updated": time.strftime("%Y-%m-%d %H:%M:%S"),
            "batches": {}
        }

        for i in range(num_batches):
            checkpoint["batches"][str(i)] = {
                "status": "pending",
                "batch_number": i,
                "result_file": None,
                "num_sequences": 0,
                "num_assignments": 0,
                "started_at": None,
                "completed_at": None,
                "error": None
            }

        # Sauvegarder le checkpoint initial
        with open(checkpoint_path, 'w') as f:
            json.dump(checkpoint, f, indent=2)

        return checkpoint

    def _load_checkpoint(self, checkpoint_path: str) -> Optional[Dict]:
        """Charge un fichier checkpoint existant.

        Args:
            checkpoint_path: Chemin du fichier checkpoint

        Returns:
            Dictionnaire du checkpoint ou None si le fichier n'existe pas
        """
        if not os.path.exists(checkpoint_path):
            return None

        try:
            with open(checkpoint_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            self.logger.error(f"Erreur lors du chargement du checkpoint: {e}")
            return None

    def _update_checkpoint(self, checkpoint_path: str, checkpoint: Dict):
        """Met à jour le fichier checkpoint.

        Args:
            checkpoint_path: Chemin du fichier checkpoint
            checkpoint: Dictionnaire du checkpoint
        """
        checkpoint["last_updated"] = time.strftime("%Y-%m-%d %H:%M:%S")

        with open(checkpoint_path, 'w') as f:
            json.dump(checkpoint, f, indent=2)

    def _save_batch_result(self, batch_df: pd.DataFrame, batch_dir: str, batch_number: int) -> str:
        """Sauvegarde le résultat d'un batch dans un fichier CSV.

        Args:
            batch_df: DataFrame avec les résultats du batch
            batch_dir: Dossier pour sauvegarder les résultats
            batch_number: Numéro du batch

        Returns:
            Chemin du fichier sauvegardé
        """
        os.makedirs(batch_dir, exist_ok=True)
        batch_file = os.path.join(batch_dir, f"batch_{batch_number:04d}_results.csv")
        batch_df.to_csv(batch_file, index=False)
        return batch_file

    def _load_batch_results(self, batch_dir: str, checkpoint: Dict) -> List[pd.DataFrame]:
        """Charge tous les résultats de batches complétés depuis les fichiers.

        Args:
            batch_dir: Dossier contenant les résultats des batches
            checkpoint: Dictionnaire du checkpoint

        Returns:
            Liste de DataFrames avec les résultats
        """
        results = []

        for batch_id, batch_info in checkpoint["batches"].items():
            if batch_info["status"] == "completed" and batch_info["result_file"]:
                result_file = batch_info["result_file"]
                if os.path.exists(result_file):
                    try:
                        df = pd.read_csv(result_file)
                        results.append(df)
                        self.logger.info(f"Batch {batch_id} chargé: {len(df)} assignations")
                    except Exception as e:
                        self.logger.error(f"Erreur lors du chargement du batch {batch_id}: {e}")

        return results

    def _get_checkpoint_status(self, checkpoint: Dict) -> Dict[str, Any]:
        """Obtient un résumé du statut du checkpoint.

        Args:
            checkpoint: Dictionnaire du checkpoint

        Returns:
            Dictionnaire avec les statistiques du checkpoint
        """
        statuses = {
            "pending": 0,
            "in_progress": 0,
            "completed": 0,
            "failed": 0
        }

        total_sequences = 0
        total_assignments = 0

        for batch_info in checkpoint["batches"].values():
            status = batch_info["status"]
            statuses[status] = statuses.get(status, 0) + 1
            total_sequences += batch_info.get("num_sequences", 0)
            total_assignments += batch_info.get("num_assignments", 0)

        return {
            "total_batches": checkpoint["total_batches"],
            "status_counts": statuses,
            "total_sequences": total_sequences,
            "total_assignments": total_assignments,
            "created_at": checkpoint.get("created_at"),
            "last_updated": checkpoint.get("last_updated")
        }

    def _cleanup_checkpoint(self, batch_dir: str):
        """Nettoie les fichiers de checkpoint et de batches après succès complet.

        Args:
            batch_dir: Dossier contenant les résultats des batches
        """
        checkpoint_path = os.path.join(batch_dir, 'checkpoint.json')

        if not os.path.exists(checkpoint_path):
            return

        checkpoint = self._load_checkpoint(checkpoint_path)

        if checkpoint is None:
            return

        # Vérifier si tous les batches sont complétés
        all_completed = all(
            batch_info["status"] == "completed"
            for batch_info in checkpoint["batches"].values()
        )

        if all_completed:
            self.logger.info(f"Tous les batches sont complétés. Archivage du checkpoint...")

            # Archiver le checkpoint au lieu de le supprimer
            archive_path = checkpoint_path.replace('.json', '_completed.json')
            os.rename(checkpoint_path, archive_path)

            self.logger.info(f"Checkpoint archivé: {archive_path}")
        else:
            self.logger.info("Certains batches ne sont pas complétés, conservation du checkpoint")

    def _reset_failed_batches(self, checkpoint_path: str):
        """Réinitialise les batches échoués pour permettre un nouveau traitement.

        Args:
            checkpoint_path: Chemin du fichier checkpoint
        """
        checkpoint = self._load_checkpoint(checkpoint_path)

        if checkpoint is None:
            return

        reset_count = 0
        for batch_info in checkpoint["batches"].values():
            if batch_info["status"] == "failed":
                batch_info["status"] = "pending"
                batch_info["error"] = None
                batch_info["started_at"] = None
                batch_info["completed_at"] = None
                reset_count += 1

        if reset_count > 0:
            self._update_checkpoint(checkpoint_path, checkpoint)
            self.logger.info(f"{reset_count} batches échoués réinitialisés à 'pending'")
        else:
            self.logger.info("Aucun batch échoué à réinitialiser")

    def _process_classification_batches(
        self, barcode: str, batches: List[pd.DataFrame],
        combined_hmm_path: str, logger, batch_dir: Optional[str] = None
    ) -> List[pd.DataFrame]:
        """Traite les lots de classification en parallèle avec hmmscan et checkpoint.

        Stratégie optimisée avec gestion d'erreurs robuste et système de checkpoint:
        - Plusieurs processus Python en parallèle (classification_max_workers)
        - Chaque processus lance 1 hmmscan avec classification_threads_per_worker threads
        - Sauvegarde le résultat de chaque batch dans un fichier séparé
        - Utilise un checkpoint JSON pour suivre l'état de chaque batch
        - Permet de reprendre le traitement en cas d'interruption

        Args:
            barcode: Code-barres en cours de traitement
            batches: Liste de DataFrames (batches de séquences)
            combined_hmm_path: Chemin vers la base HMM combinée
            logger: Logger pour les messages
            batch_dir: Dossier pour sauvegarder les résultats des batches

        Returns:
            Liste de DataFrames contenant les assignations de chaque batch
        """
        total_batches = len(batches)

        # Définir le dossier de batch par défaut
        if batch_dir is None:
            batch_dir = os.path.join(self.output_folder, barcode, 'batches')

        os.makedirs(batch_dir, exist_ok=True)

        # Chemin du fichier checkpoint
        checkpoint_path = os.path.join(batch_dir, 'checkpoint.json')

        # Charger ou initialiser le checkpoint
        checkpoint = self._load_checkpoint(checkpoint_path)

        if checkpoint is None:
            logger.info(f"{barcode} - Initialisation d'un nouveau checkpoint pour {total_batches} batches")
            checkpoint = self._init_checkpoint(checkpoint_path, total_batches, barcode)
        else:
            logger.info(f"{barcode} - Checkpoint existant trouvé, reprise du traitement")
            # Charger les résultats déjà complétés
            completed_results = self._load_batch_results(batch_dir, checkpoint)
            logger.info(f"{barcode} - {len(completed_results)} batches déjà complétés chargés")

        logger.info(f"{barcode} - Configuration parallélisation:")
        logger.info(f"  - {self.classification_max_workers} processus en parallèle")
        logger.info(f"  - {self.classification_threads_per_worker} threads par processus")
        logger.info(f"  - Total: {self.classification_max_workers * self.classification_threads_per_worker} threads effectifs")

        # Déterminer quels batches doivent être traités
        batches_to_process = []
        for i, batch_df in enumerate(batches):
            batch_status = checkpoint["batches"][str(i)]["status"]
            if batch_status != "completed":
                batches_to_process.append((i, batch_df))

        logger.info(f"{barcode} - {len(batches_to_process)} batches à traiter (sur {total_batches})")

        if not batches_to_process:
            logger.info(f"{barcode} - Tous les batches sont déjà complétés")
            return self._load_batch_results(batch_dir, checkpoint)

        # Traiter les batches en parallèle
        with ProcessPoolExecutor(max_workers=self.classification_max_workers) as executor:
            futures = {}
            for batch_idx, batch_df in batches_to_process:
                # Marquer le batch comme en cours
                checkpoint["batches"][str(batch_idx)]["status"] = "in_progress"
                checkpoint["batches"][str(batch_idx)]["started_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
                checkpoint["batches"][str(batch_idx)]["num_sequences"] = len(batch_df)

                future = executor.submit(
                    external_tools.run_hmmscan,
                    batch_df,
                    combined_hmm_path,
                    self.seq_col,
                    self.classification_threads_per_worker,
                    self.tax_columns,
                    self.id_col
                )
                futures[future] = batch_idx

            # Mettre à jour le checkpoint initial
            self._update_checkpoint(checkpoint_path, checkpoint)

            # Traiter les résultats au fur et à mesure
            for future in tqdm(
                concurrent.futures.as_completed(futures),
                total=len(futures),
                desc=f"Classification des lots pour {barcode}",
                unit="lot"
            ):
                batch_idx = futures[future]
                batch_info = checkpoint["batches"][str(batch_idx)]

                try:
                    batch_assignments = future.result()

                    if not batch_assignments.empty:
                        # Sauvegarder le résultat du batch
                        result_file = self._save_batch_result(batch_assignments, batch_dir, batch_idx)

                        # Mettre à jour le checkpoint
                        batch_info["status"] = "completed"
                        batch_info["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
                        batch_info["result_file"] = result_file
                        batch_info["num_assignments"] = len(batch_assignments)

                        logger.info(f"{barcode} - Batch {batch_idx+1}/{total_batches}: {len(batch_assignments)} séquences assignées")
                    else:
                        # Aucune assignation, marquer comme complété quand même
                        batch_info["status"] = "completed"
                        batch_info["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
                        batch_info["num_assignments"] = 0

                        logger.info(f"{barcode} - Batch {batch_idx+1}/{total_batches}: aucune assignation")

                except Exception as e:
                    # Marquer le batch comme échoué
                    batch_info["status"] = "failed"
                    batch_info["error"] = str(e)
                    batch_info["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")

                    logger.error(f"{barcode} - Erreur lors du traitement du batch {batch_idx+1}/{total_batches}: {e}")

                # Mettre à jour le checkpoint après chaque batch
                self._update_checkpoint(checkpoint_path, checkpoint)

        # Charger tous les résultats (y compris ceux précédemment complétés)
        all_assignments = self._load_batch_results(batch_dir, checkpoint)

        # Afficher un résumé
        completed = sum(1 for b in checkpoint["batches"].values() if b["status"] == "completed")
        failed = sum(1 for b in checkpoint["batches"].values() if b["status"] == "failed")
        logger.info(f"{barcode} - Résumé: {completed} complétés, {failed} échoués sur {total_batches} batches")

        return all_assignments


    def _empty_result(self, barcode: str) -> Dict[str, Any]:
        """Retourne un résultat vide pour un barcode sans données."""
        return {
            "barcode": barcode,
            "unprocessed_df": pd.DataFrame(),
            "hmm_db_path": None,
            "msa_paths": {},
            "tasks_processed": 0,
            "sequences_processed": 0,
            "sequences_unprocessed": 0
        }


    def _build_result(
        self, barcode: str, unassigned_df: pd.DataFrame, 
        combined_hmm_path: Optional[str], msa_output_paths: Dict,
        all_tasks: List[Dict], barcode_df: pd.DataFrame,
        classification_sequences: Optional[int] = None 
    ) -> Dict[str, Any]:
        return {
            "barcode": barcode,
            "unprocessed_df": unassigned_df,
            "hmm_db_path": combined_hmm_path,
            "msa_paths": msa_output_paths,
            "tasks_processed": len(all_tasks),
            "sequences_processed": len(barcode_df),
            "sequences_unprocessed": len(unassigned_df),
            "classification_sequences": int(classification_sequences or 0),
        }

    def _load_existing_data(
        self, barcode: str, barcode_output_dir: str, start_step: int, logger
    ) -> Tuple[pd.DataFrame, pd.DataFrame, List[Dict], Dict, Optional[str]]:
        """Charge les données existantes depuis les fichiers pour reprendre le traitement.
        
        Args:
            barcode: Le barcode à traiter
            barcode_output_dir: Le dossier de sortie du barcode
            start_step: L'étape de départ
            logger: Le logger
        
        Returns:
            Tuple contenant (barcode_df, unprocessed_df, all_tasks, msa_output_paths, combined_hmm_path)
        """
        barcode_df = pd.DataFrame()
        unprocessed_df = pd.DataFrame()
        all_tasks = []
        msa_output_paths = {}
        combined_hmm_path = None
        
        # Charger barcode_df depuis les données originales si nécessaire (étapes 2+)
        if start_step >= 2:
            barcode_df = self.data_df[self.data_df[self.barcode_col] == barcode]
            if barcode_df.empty:
                logger.warning(f"{barcode} - Aucune donnée trouvée dans le dataset original")
        
        # Charger unprocessed_df (nécessaire pour étapes 2+)
        if start_step >= 2:
            unprocessed_path = os.path.join(barcode_output_dir, 'unprocessed_after_trimal.csv')
            if os.path.exists(unprocessed_path):
                unprocessed_df = pd.read_csv(unprocessed_path)
                logger.info(f"{barcode} - {len(unprocessed_df)} séquences non traitées chargées")
            else:
                logger.warning(f"{barcode} - Fichier unprocessed_after_trimal.csv non trouvé")
        
        # Charger les MSA paths (nécessaire pour étapes 2+)
        if start_step >= 2:
            msa_output_paths = self._load_msa_paths(barcode, barcode_output_dir, logger)
        
        # Charger combined_hmm_path (nécessaire pour étapes 3+)
        # Note: Avec la nouvelle architecture, l'étape 3 génère le HMM
        if start_step >= 3:
            combined_hmm_path = os.path.join(barcode_output_dir, 'all_profiles.hmm')
            if os.path.exists(combined_hmm_path):
                logger.info(f"{barcode} - Base de données HMM chargée: {combined_hmm_path}")
            else:
                logger.warning(f"{barcode} - Fichier all_profiles.hmm non trouvé")
                combined_hmm_path = None
        
        return barcode_df, unprocessed_df, all_tasks, msa_output_paths, combined_hmm_path


    def _load_msa_paths(self, barcode: str, barcode_output_dir: str, logger) -> Dict:
        """Charge les chemins des fichiers MSA existants.

        Returns:
            Dictionnaire avec les clés (clade_name, level) et les chemins des MSA
        """
        msa_output_paths = {}
        msa_dir = os.path.join(self.output_folder, barcode)

        if not os.path.exists(msa_dir):
            logger.warning(f"{barcode} - Dossier MSA introuvable : {msa_dir}")
            return msa_output_paths

        logger.info(f"{barcode} - Recherche des fichiers MSA dans {msa_dir}")

        for root, _, files in os.walk(msa_dir):
            folder_name = os.path.basename(root).lower()

            # Identifier le niveau taxonomique à partir du nom du dossier
            if "species" in folder_name:
                level = "species"
            elif "genus" in folder_name:
                level = "genus"
            elif "family" in folder_name:
                level = "family"
            elif "order" in folder_name:
                level = "order"
            elif "class" in folder_name:
                level = "class"
            elif "phylum" in folder_name:
                level = "phylum"
            else:
                continue  # ignorer les autres dossiers

            for filename in files:
                if not filename.endswith(".fasta"):
                    continue

                msa_path = os.path.join(root, filename)

                # Extraire le clade_name avant le suffixe "_msa.fasta"
                if filename.endswith("_msa.fasta"):
                    clade_name = filename[:-len("_msa.fasta")]
                else:
                    clade_name = filename.replace(".fasta", "")

                msa_output_paths[(clade_name, level)] = msa_path

        logger.info(f"{barcode} - {len(msa_output_paths)} fichiers MSA détectés")
        return msa_output_paths

    def _collect_all_unprocessed(self, barcode_results: Dict) -> pd.DataFrame:
        """
        Collecte toutes les séquences non traitées de tous les barcodes.

        Args:
            barcode_results: Résultats de la phase 1 pour chaque barcode

        Returns:
            DataFrame avec toutes les séquences non traitées et la colonne 'source_barcode'
        """
        self.logger.info("Collecte de toutes les séquences non traitées...")

        unprocessed_sequences = []
        for barcode, data in barcode_results.items():
            unprocessed_df = data.get('unprocessed_df', pd.DataFrame())

            if not unprocessed_df.empty:
                # Ajouter la colonne source_barcode pour traçabilité
                unprocessed_df = unprocessed_df.copy()
                unprocessed_df['source_barcode'] = barcode
                unprocessed_sequences.append(unprocessed_df)
                self.logger.info(f"{barcode} - {len(unprocessed_df)} séquences non traitées")

        if not unprocessed_sequences:
            self.logger.warning("Aucune séquence non traitée trouvée")
            return pd.DataFrame()

        all_unprocessed = pd.concat(unprocessed_sequences, ignore_index=True)
        self.logger.info(f"Total: {len(all_unprocessed)} séquences non traitées à classifier")

        # Sauvegarder pour traçabilité
        combined_dir = os.path.join(self.output_folder, 'combined')
        os.makedirs(combined_dir, exist_ok=True)
        all_unprocessed.to_csv(
            os.path.join(combined_dir, 'all_unprocessed_sequences.csv'),
            index=False
        )

        return all_unprocessed

    def _combine_all_hmm_databases(self, barcode_results: Dict) -> Optional[str]:
        """
        Combine toutes les bases HMM en une seule base unifiée.

        Args:
            barcode_results: Résultats de la phase 1 pour chaque barcode

        Returns:
            Chemin vers la base HMM combinée, ou None si échec
        """
        self.logger.info("Combinaison de toutes les bases HMM...")

        hmm_paths = []
        for barcode, data in barcode_results.items():
            hmm_db_path = data.get('hmm_db_path')
            if hmm_db_path and os.path.exists(hmm_db_path):
                hmm_paths.append(hmm_db_path)
                self.logger.info(f"{barcode} - Base HMM: {hmm_db_path}")
            else:
                self.logger.warning(f"{barcode} - Aucune base HMM trouvée")

        if not hmm_paths:
            self.logger.error("Aucune base HMM disponible pour la combinaison")
            return None

        # Créer le dossier combined
        combined_dir = os.path.join(self.output_folder, 'combined')
        os.makedirs(combined_dir, exist_ok=True)

        combined_hmm_path = os.path.join(combined_dir, 'all_profiles_combined.hmm')

        # Combiner les bases HMM
        success = external_tools.combine_hmm_databases(hmm_paths, combined_hmm_path)

        if success:
            self.logger.info(f"Base HMM combinée créée: {combined_hmm_path}")
            return combined_hmm_path
        else:
            self.logger.error("Échec de la combinaison des bases HMM")
            return None

    def _classify_unified(self, all_unprocessed: pd.DataFrame,
                         combined_hmm_path: str) -> pd.DataFrame:
        """
        Classification unifiée de toutes les séquences non traitées avec checkpoint.

        Args:
            all_unprocessed: DataFrame avec toutes les séquences non traitées
            combined_hmm_path: Chemin vers la base HMM combinée

        Returns:
            DataFrame avec les assignations
        """
        if all_unprocessed.empty or not combined_hmm_path:
            self.logger.warning("Aucune séquence à classifier ou base HMM manquante")
            return pd.DataFrame()

        self.logger.info(f"Classification unifiée de {len(all_unprocessed)} séquences...")

        # Définir le dossier pour les batches unifiés
        combined_dir = os.path.join(self.output_folder, 'combined')
        batch_dir = os.path.join(combined_dir, 'batches')
        os.makedirs(batch_dir, exist_ok=True)

        # Utiliser le même système de batching que dans la version originale
        batches = self._create_batches(all_unprocessed, self.classification_batch_size)
        total_batches = len(batches)

        self.logger.info(f"Classification en {total_batches} lots de ~{self.classification_batch_size} séquences")

        # Traiter les lots en parallèle avec checkpoint
        all_assignments = self._process_classification_batches(
            "UNIFIED", batches, combined_hmm_path, self.logger, batch_dir
        )

        # Consolider les résultats
        if not all_assignments:
            self.logger.warning("Aucune assignation trouvée")
            return pd.DataFrame()

        self.logger.info(f"Concaténation de {len(all_assignments)} résultats de batches...")
        assignments_df = pd.concat(all_assignments, ignore_index=True)
        assignments_df.drop_duplicates(subset=[self.id_col], keep='first', inplace=True)

        # Ajouter la colonne source_barcode depuis all_unprocessed
        assignments_df = assignments_df.merge(
            all_unprocessed[[self.id_col, 'source_barcode']],
            on=self.id_col,
            how='left'
        )

        self.logger.info(f"Total: {len(assignments_df)} séquences assignées après déduplication")

        # Sauvegarder les assignations globales consolidées
        global_assignments_path = os.path.join(combined_dir, 'global_assignments.csv')
        assignments_df.to_csv(global_assignments_path, index=False)
        self.logger.info(f"Assignations globales sauvegardées: {global_assignments_path}")

        # Nettoyer le checkpoint après succès
        self._cleanup_checkpoint(batch_dir)

        return assignments_df

    def _distribute_assignments_by_barcode(self, assignments_df: pd.DataFrame,
                                          all_unprocessed: pd.DataFrame,
                                          barcode_results: Dict):
        """
        Distribue les assignations par barcode source et intègre dans les MSA.

        Args:
            assignments_df: DataFrame avec toutes les assignations
            all_unprocessed: DataFrame avec toutes les séquences non traitées
            barcode_results: Résultats de la phase 1
        """
        if assignments_df.empty:
            self.logger.info("Aucune assignation à distribuer")
            return

        self.logger.info("Distribution des assignations par barcode...")

        for barcode in self.barcodes:
            # Filtrer les assignations pour ce barcode
            barcode_assignments = assignments_df[
                assignments_df['source_barcode'] == barcode
            ].copy()

            if barcode_assignments.empty:
                self.logger.info(f"{barcode} - Aucune assignation")
                continue

            self.logger.info(f"{barcode} - {len(barcode_assignments)} séquences assignées")

            # Sauvegarder les assignations par barcode
            barcode_output_dir = os.path.join(self.output_folder, barcode)
            os.makedirs(barcode_output_dir, exist_ok=True)

            barcode_assignments.to_csv(
                os.path.join(barcode_output_dir, 'global_assignments.csv'),
                index=False
            )

            # Extraire le barcode cible pour déterminer à quels MSA ajouter les séquences
            barcode_assignments['target_barcode'] = barcode_assignments['assigned_msa'].str.extract(
                r'^(ITS\d?)_'
            )[0]

            # Grouper par barcode cible
            for target_barcode, group_df in barcode_assignments.groupby('target_barcode'):
                if pd.isna(target_barcode):
                    continue

                target_data = barcode_results.get(target_barcode)
                if not target_data:
                    continue

                msa_paths = target_data.get('msa_paths', {})
                if not msa_paths:
                    continue

                # Obtenir les séquences complètes depuis all_unprocessed
                seq_ids = group_df[self.id_col].tolist()
                sequences_to_add = all_unprocessed[
                    all_unprocessed[self.id_col].isin(seq_ids)
                ]

                self.logger.info(
                    f"{barcode} → {target_barcode}: {len(sequences_to_add)} séquences à intégrer"
                )

                # Intégrer dans les MSA du barcode cible
                self._add_sequences_to_msa(
                    assignments_df=group_df,
                    unprocessed_df=sequences_to_add,
                    barcode=target_barcode,
                    msa_output_paths=msa_paths
                )

    def _phase1_prepare_msa_and_hmm(self) -> Dict[str, Any]:
        """
        Phase 1 : Création des MSA et HMM par barcode (code existant).

        Returns:
            Dictionnaire des résultats par barcode
        """
        self.logger.info("PHASE 1: Préparation des MSA et HMM par barcode")

        barcode_results = {}

        with ThreadPoolExecutor(max_workers=self.max_barcode_workers) as executor:
            futures = [
                executor.submit(self._process_barcode, barcode)
                for barcode in self.barcodes
            ]

            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result and result.get("barcode"):
                        barcode_results[result["barcode"]] = result
                        self.logger.info(
                            f"{result['barcode']} - Phase 1 terminée"
                        )
                except Exception as e:
                    self.logger.error(
                        f"Erreur majeure lors du traitement d'un barcode: {str(e)}",
                        exc_info=True
                    )

        return barcode_results

    def _phase2_unified_classification(self, barcode_results: Dict) -> Dict[str, Any]:
        """
        Phase 2 : Classification unifiée sur tous les barcodes.

        Args:
            barcode_results: Résultats de la phase 1

        Returns:
            Dictionnaire avec les assignations et statistiques
        """
        self.logger.info("PHASE 2: Classification unifiée inter-barcodes")

        # 2.1. Collecter toutes les séquences non traitées
        all_unprocessed = self._collect_all_unprocessed(barcode_results)

        if all_unprocessed.empty:
            self.logger.info("Aucune séquence non traitée, phase 2 terminée")
            return {
                'total_assignments': 0,
                'reassignment_stats': {}
            }

        # 2.2. Combiner toutes les bases HMM
        combined_hmm_path = self._combine_all_hmm_databases(barcode_results)

        if not combined_hmm_path:
            self.logger.error("Impossible de continuer sans base HMM combinée")
            return {
                'total_assignments': 0,
                'reassignment_stats': {}
            }

        # 2.3. Classification unifiée
        assignments = self._classify_unified(all_unprocessed, combined_hmm_path)

        # 2.4. Distribuer les assignations par barcode
        self._distribute_assignments_by_barcode(
            assignments, all_unprocessed, barcode_results
        )

        # 2.5. Générer les statistiques de réaffectation
        reassignment_stats = generate_reassignment_statistics(
            assignments, self.id_col
        )

        self.logger.info("PHASE 2: Classification unifiée terminée")

        return {
            'total_assignments': len(assignments),
            'reassignment_stats': reassignment_stats
        }



    def _process_barcode(self, barcode: str, start_step: int = 1) -> Dict[str, Any]:
        """
        Traite toutes les séquences pour un barcode spécifique (Phase 1 uniquement).

        Avec la nouvelle architecture unifiée, cette méthode s'arrête à l'étape 3.
        La classification et l'intégration sont maintenant gérées par la Phase 2.

        Args:
            barcode: Le barcode à traiter
            start_step: L'étape de départ (1-3). Par défaut 1 pour tout exécuter.
                    1: Préparation et MSA
                    2: Nettoyage Trimal
                    3: Génération HMM
        """
        logger = self._get_barcode_logger(barcode)
        logger.info(f" {barcode} - Début du traitement Phase 1 (à partir de l'étape {start_step})")

        # Création du dossier de sortie
        barcode_output_dir = os.path.join(self.output_folder, barcode)
        os.makedirs(barcode_output_dir, exist_ok=True)

        # Initialisation des variables
        barcode_df = pd.DataFrame()
        unprocessed_df = pd.DataFrame()
        all_tasks = []
        msa_output_paths = {}
        combined_hmm_path = None

        # Récupération des données existantes si on ne démarre pas au début
        if start_step > 1:
            logger.info(f"{barcode} - Récupération des données existantes...")
            barcode_df, unprocessed_df, all_tasks, msa_output_paths, combined_hmm_path = \
                self._load_existing_data(barcode, barcode_output_dir, start_step, logger)

        # Étape 1: Préparation des données et génération des MSA
        if start_step <= 1:
            barcode_df, unprocessed_df, all_tasks = self._prepare_barcode_data(barcode, logger)
            if barcode_df.empty:
                return self._empty_result(barcode)

            msa_output_paths = self._generate_msa(barcode, all_tasks, logger)

        # Étape 2: Nettoyage des MSA avec Trimal
        if start_step <= 2:
            unprocessed_df = self._clean_msa_with_trimal(
                barcode, barcode_df, barcode_output_dir, msa_output_paths, all_tasks, unprocessed_df, logger
            )

        # Étape 3: Génération des HMM
        if start_step <= 3:
            combined_hmm_path = self._generate_hmm_database(
                barcode, barcode_output_dir, msa_output_paths, logger
            )

        logger.info(f"{barcode} - Phase 1 terminée")

        return self._build_result(
            barcode, unprocessed_df, combined_hmm_path, msa_output_paths,
            all_tasks, barcode_df, classification_sequences=0
        )