#!/usr/bin/env python3
"""
Script de déduplication des taxonomies
Détecte et supprime les doublons taxonomiques, met à jour les relations seq_tax.
"""

import os
import sys
import time
import logging
import sqlite3
import json
from datetime import datetime
from typing import List, Tuple

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TaxonomyDeduplicator:
    """Gestionnaire de déduplication des taxonomies."""
    
    def __init__(self, db_path: str):
        self.db_path = db_path
        self.stats = {
            'initial_taxonomy_count': 0,
            'final_taxonomy_count': 0,
            'duplicate_groups': 0,
            'duplicates_removed': 0,
            'sequences_updated': 0,
            'processing_time': 0.0
        }

    def get_connection(self) -> sqlite3.Connection:
        """Obtient une connexion à la base de données avec des PRAGMAs."""
        conn = sqlite3.connect(self.db_path, timeout=30)
        # Améliore les performances d'écriture et la concurrence
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        conn.execute("PRAGMA foreign_keys = ON")
        # Augmente la taille du cache pour réduire les I/O disque
        conn.execute("PRAGMA cache_size = -64000") # 64MB
        return conn

    def setup_database(self):
        """Crée les index nécessaires pour la performance si ils n'existent pas."""
        logger.info("Vérification et création des index pour l'optimisation...")
        with self.get_connection() as conn:
            cursor = conn.cursor()
            # Index crucial pour accélérer les UPDATE sur seq_tax
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_seq_tax_tax_id ON seq_tax(tax_id)")
            # Index composite pour accélérer le GROUP BY sur la table taxonomy
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_taxonomy_full ON taxonomy
                (kingdom, phylum, class_, order_, family, genus, species)
            """)
            conn.commit()
        logger.info("Index prêts.")

    def find_duplicate_taxonomies(self) -> List[Tuple[int, List[int]]]:
        """
        Trouve tous les doublons taxonomiques en utilisant une seule requête SQL.
        Retourne une liste de tuples: (tax_id_a_garder, [liste_des_ids_a_supprimer])
        """
        logger.info("Recherche des doublons taxonomiques via SQL...")
        
        duplicate_groups = []
        
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # OPTIMISATION 1: Utiliser GROUP BY et HAVING pour trouver les doublons directement en SQL.
            # GROUP_CONCAT assemble tous les tax_id d'un même groupe.
            query = """
                SELECT
                    MIN(tax_id) as keeper_id,
                    GROUP_CONCAT(tax_id) as all_ids,
                    COUNT(*) as count
                FROM taxonomy
                GROUP BY kingdom, phylum, class_, order_, family, genus, species
                HAVING COUNT(*) > 1
            """
            
            cursor.execute(query)
            
            for row in cursor.fetchall():
                keeper_id = row[0]
                # Convertir la chaîne "1,2,3" en une liste d'entiers
                all_ids = [int(id_str) for id_str in row[1].split(',')]
                # Les ID à supprimer sont tous les ID du groupe sauf celui qu'on garde
                ids_to_remove = [id_val for id_val in all_ids if id_val != keeper_id]
                
                if ids_to_remove:
                    duplicate_groups.append((keeper_id, ids_to_remove))

        self.stats['duplicate_groups'] = len(duplicate_groups)
        self.stats['duplicates_removed'] = sum(len(ids) for _, ids in duplicate_groups)
        
        logger.info(f"Trouvé {self.stats['duplicate_groups']} groupes de doublons "
                   f"({self.stats['duplicates_removed']} enregistrements à supprimer)")
        
        return duplicate_groups


    def merge_duplicate_taxonomies(self, duplicate_groups: List[Tuple[int, List[int]]]):
        """
        Fusionne les taxonomies dupliquées en utilisant la méthode haute performance
        "Créer-Remplir-Échanger" (Create-Populate-Swap), idéale pour les très grandes tables seq_tax.
        """
        logger.info("Fusion (méthode haute performance) des taxonomies dupliquées...")
        
        if not duplicate_groups:
            return

        # 1. Créer une table de correspondance complète de chaque ID à supprimer vers son ID à garder
        id_remap = {
            id_to_remove: keeper_id
            for keeper_id, ids_to_remove in duplicate_groups
            for id_to_remove in ids_to_remove
        }
        
        # IDs qui seront finalement supprimés de la table `taxonomy`
        all_ids_to_remove = list(id_remap.keys())
        
        total_sequences_updated = 0

        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            try:
                # Démarrage d'une transaction globale
                cursor.execute("BEGIN")

                # 2. Créer une nouvelle table vide avec la même structure
                logger.info("Création de la nouvelle table 'seq_tax_new'...")
                cursor.execute("DROP TABLE IF EXISTS seq_tax_new")
                create_sql = cursor.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='seq_tax'").fetchone()[0]
                cursor.execute(create_sql.replace("seq_tax", "seq_tax_new"))

                # 3. Lire l'ancienne table et insérer dans la nouvelle par lots
                logger.info("Migration des données vers la nouvelle table...")
                
                read_cursor = conn.cursor()
                read_cursor.execute("SELECT seq_id, tax_id FROM seq_tax")
                
                batch_size = 100000
                batch = []
                
                for seq_id, tax_id in read_cursor:
                    new_tax_id = id_remap.get(tax_id, tax_id)
                    if new_tax_id != tax_id:
                        total_sequences_updated += 1
                    
                    batch.append((seq_id, new_tax_id))
                    
                    if len(batch) >= batch_size:
                        cursor.executemany("INSERT INTO seq_tax_new (seq_id, tax_id) VALUES (?, ?)", batch)
                        batch = []
                
                if batch:
                    cursor.executemany("INSERT INTO seq_tax_new (seq_id, tax_id) VALUES (?, ?)", batch)

                # 4. Échanger les tables
                logger.info("Remplacement de l'ancienne table...")
                cursor.execute("DROP TABLE seq_tax")
                cursor.execute("ALTER TABLE seq_tax_new RENAME TO seq_tax")

                # 5. Recréer l'index crucial pour les performances futures
                logger.info("Création de l'index sur la nouvelle table...")
                cursor.execute("CREATE INDEX idx_seq_tax_tax_id ON seq_tax(tax_id)")
                
                # 6. Supprimer les taxonomies dupliquées de la table `taxonomy`
                logger.info(f"Suppression de {len(all_ids_to_remove)} taxonomies obsolètes en plusieurs lots...")

                delete_batch_size = 10000  # Une taille de lot sûre, bien en dessous de la limite de SQLite

                for i in range(0, len(all_ids_to_remove), delete_batch_size):
                    # On prend une "tranche" de la liste des IDs à supprimer
                    batch_to_delete = all_ids_to_remove[i:i + delete_batch_size]
                    
                    # On crée les placeholders uniquement pour ce lot
                    placeholders = ', '.join(['?'] * len(batch_to_delete))
                    delete_query = f"DELETE FROM taxonomy WHERE tax_id IN ({placeholders})"
                    
                    # On exécute la suppression pour le lot actuel
                    cursor.execute(delete_query, tuple(batch_to_delete))

                logger.info("Toutes les taxonomies obsolètes ont été supprimées.")
                
                # Validation de toutes les opérations
                conn.commit()
                logger.info("Fusion et reconstruction terminées avec succès.")

            except Exception as e:
                conn.rollback()
                logger.error(f"Erreur lors de la fusion : {e}", exc_info=True)
                raise

            self.stats['sequences_updated'] = total_sequences_updated
    
    def verify_integrity(self) -> bool:
        """Vérifie l'intégrité après la déduplication."""
        logger.info("Vérification de l'intégrité finale...")
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # Vérifier qu'il n'y a plus de doublons
            cursor.execute("""
                SELECT COUNT(*) FROM (
                    SELECT 1 FROM taxonomy
                    GROUP BY kingdom, phylum, class_, order_, family, genus, species
                    HAVING COUNT(*) > 1
                )
            """)
            remaining_duplicates = cursor.fetchone()[0]

            # Vérifier les relations orphelines
            cursor.execute("""
                SELECT COUNT(*) FROM seq_tax st
                LEFT JOIN taxonomy t ON st.tax_id = t.tax_id
                WHERE t.tax_id IS NULL
            """)
            orphaned_relations = cursor.fetchone()[0]

            if remaining_duplicates == 0 and orphaned_relations == 0:
                logger.info("✓ Intégrité vérifiée avec succès.")
                return True
            else:
                logger.warning(f"⚠ Problèmes d'intégrité détectés: "
                               f"{remaining_duplicates} doublons restants, "
                               f"{orphaned_relations} relations orphelines.")
                return False

    def optimize_database(self):
        """Optimise la base de données après la déduplication."""
        logger.info("Optimisation de la base de données (VACUUM)...")
        with self.get_connection() as conn:
            conn.execute("VACUUM")
        logger.info("Optimisation terminée.")
    
    def deduplicate_taxonomies(self):
        """Processus complet de déduplication."""
        start_time = time.time()
        logger.info("Démarrage de la déduplication taxonomique")
        self.stats['final_unique_tax_in_seq'] = 0
        
        try:
            with self.get_connection() as conn:
                self.stats['initial_taxonomy_count'] = conn.execute("SELECT COUNT(*) FROM taxonomy").fetchone()[0]
            logger.info(f"Nombre initial de taxonomies : {self.stats['initial_taxonomy_count']}")

            # 0. Préparation : création des index
            self.setup_database()

            # 1. Trouver les doublons (version SQL rapide)
            duplicate_groups = self.find_duplicate_taxonomies()
            
            if not duplicate_groups:
                logger.info("Aucun doublon trouvé.")
                self.stats['processing_time'] = time.time() - start_time
                return
            
            # 2. Fusionner les doublons (version en masse)
            self.merge_duplicate_taxonomies(duplicate_groups)
            
            # 3. Vérifier l'intégrité
            integrity_ok = self.verify_integrity()
            
            # 4. Mettre à jour les stats finales
            with self.get_connection() as conn:
                self.stats['final_taxonomy_count'] = conn.execute("SELECT COUNT(*) FROM taxonomy").fetchone()[0]
                self.stats['final_unique_tax_in_seq'] = conn.execute("SELECT COUNT(DISTINCT tax_id) FROM seq_tax").fetchone()[0]
                logger.info(f"Nombre final de tax_id uniques dans seq_tax : {self.stats['final_unique_tax_in_seq']}")

            # 5. Optimiser l'espace disque
            self.optimize_database()
            
        except Exception as e:
            logger.error(f"Erreur durant la déduplication: {e}", exc_info=True)
            raise
        finally:
            self.stats['processing_time'] = time.time() - start_time
            logger.info(f"Déduplication terminée en {self.stats['processing_time']:.2f}s")
            logger.info(f"Taxonomies : {self.stats['initial_taxonomy_count']} -> {self.stats['final_taxonomy_count']}")
            logger.info(f"Supprimé {self.stats['duplicates_removed']} doublons, "
                       f"mis à jour {self.stats['sequences_updated']} relations.")

def write_report(stats: dict, path: str, db_path: str):
    """
    Écrit un rapport JSON récapitulant la déduplication.

    Parameters
    ----------
    stats : dict
        Les statistiques retournées par TaxonomyDeduplicator.
    path : str
        Chemin du fichier de sortie.
    db_path : str
        Chemin de la base de données traitée (pour le contexte).
    """
    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "database": os.path.abspath(db_path),
        "summary": {
            "initial_taxonomy_count": stats.get("initial_taxonomy_count"),
            "final_taxonomy_count":   stats.get("final_taxonomy_count"),
            "final_unique_tax_id_in_seq_tax":   stats.get("final_unique_tax_in_seq"),
            "duplicate_groups":       stats.get("duplicate_groups"),
            "duplicates_removed":     stats.get("duplicates_removed"),
            "sequences_updated":      stats.get("sequences_updated"),
            "processing_time_sec":    round(stats.get("processing_time", 0), 3),
        }
    }

    # On crée le dossier parent si nécessaire
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    logger.info(f"Rapport de déduplication écrit dans: {path}")




def main():
    """Point d'entrée principal."""
    try:
        if 'snakemake' in globals():
            db_path = snakemake.params.db
            output_report = snakemake.output.report
            log_file = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(file_handler)
        else:
            logger.warning("Mode test - utiliser avec Snakemake en production")
            # Créez une fausse base de données pour le test si nécessaire
            if not os.path.exists("test_database.sqlite"):
                 # code pour générer une db de test ici...
                 pass
            db_path = "test_database.sqlite"
            output_report = "taxonomy_deduplication_report.json"
        
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Base de données introuvable: {db_path}")
        
        deduplicator = TaxonomyDeduplicator(db_path)
        deduplicator.deduplicate_taxonomies()
        
        print(f"\n{'='*60}")
        print("RÉSUMÉ DE LA DÉDUPLICATION TAXONOMIQUE")
        print(f"{'='*60}")
        print(f"Taxonomies initiales: {deduplicator.stats['initial_taxonomy_count']}")
        print(f"Taxonomies finales:   {deduplicator.stats['final_taxonomy_count']}")
        print(f"Taxonomies uniques dans seq_tax: {deduplicator.stats['final_unique_tax_in_seq']}")
        print(f"Groupes de doublons:  {deduplicator.stats['duplicate_groups']}")
        print(f"Doublons supprimés:   {deduplicator.stats['duplicates_removed']}")
        print(f"Relations mises à jour: {deduplicator.stats['sequences_updated']}")
        print(f"Temps de traitement:  {deduplicator.stats['processing_time']:.2f}s")
        print(f"{'='*60}")
    
        write_report(deduplicator.stats, output_report, db_path)
        
    except Exception as e:
        logger.critical(f"Erreur fatale: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()