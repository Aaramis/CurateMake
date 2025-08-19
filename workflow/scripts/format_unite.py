#!/usr/bin/env python3
"""
Script pour formatter les données UNITE avec intégration base de données
workflow/scripts/format_unite.py

Professional UNITE database formatter with SQLite integration,
FASTA parsing, taxonomy extraction, and data validation.
"""

import os
import sys
import gzip
import re
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterator
import logging

# Import du gestionnaire de base de données
sys.path.append(str(Path(__file__).parent / "utils"))
from db_manager import DatabaseManager, get_database_manager

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class UNITEParser:
    """Parser professionnel pour les données UNITE avec intégration base de données."""
    
    def __init__(self, db_manager: DatabaseManager, config_validation: Dict):
        self.db_manager = db_manager
        self.config_validation = config_validation
        self.stats = {
            'total_sequences': 0,
            'processed_sequences': 0,
            'duplicate_sequences': 0,
            'invalid_sequences': 0,
            'parsing_errors': 0
        }
    
    def parse_unite_header(self, header: str) -> Tuple[str, Dict[str, str]]:
        """
        Parse UNITE FASTA header to extract sequence ID and taxonomy.
        
        Exemple de header UNITE:
        >SH1544431.08FU_UDB044419_refs_singleton|reps|k__Fungi|p__Basidiomycota|c__Agaricomycetes|o__Agaricales|f__Agaricaceae|g__Agaricus|s__Agaricus_bisporus
        
        Returns:
            tuple: (sequence_id, taxonomy_dict)
        """
        try:
            header = header.lstrip('>')
            
            # Séparer l'ID principal des annotations taxonomiques
            if '|' in header:
                parts = header.split('|')
                seq_id = parts[0]
                taxo_unite = parts[1]
                
                # Extraire la taxonomie des parties suivantes
                taxonomy = {
                    'kingdom': None, 'phylum': None, 'class_': None,
                    'order': None, 'family': None, 'genus': None, 'species': None
                }
                
                # Parser les rangs taxonomiques
                for part in taxo_unite.split(';'):
                    if part.startswith('k__'):
                        taxonomy['kingdom'] = part[3:].replace('_', ' ')
                    elif part.startswith('p__'):
                        taxonomy['phylum'] = part[3:].replace('_', ' ')
                    elif part.startswith('c__'):
                        taxonomy['class_'] = part[3:].replace('_', ' ')
                    elif part.startswith('o__'):
                        taxonomy['order'] = part[3:].replace('_', ' ')
                    elif part.startswith('f__'):
                        taxonomy['family'] = part[3:].replace('_', ' ')
                    elif part.startswith('g__'):
                        taxonomy['genus'] = part[3:].replace('_', ' ')
                    elif part.startswith('s__'):
                        species_full = part[3:].replace('_', ' ')
                        taxonomy['species'] = species_full
                
                # Construire le nom scientifique
                if taxonomy['genus'] and taxonomy['species']:
                    # Extraire seulement l'épithète spécifique si le genre est répété
                    species_parts = taxonomy['species'].split()
                    if len(species_parts) > 1 and species_parts[0] == taxonomy['genus']:
                        taxonomy['scientific_name'] = f"{taxonomy['genus']} {species_parts[1]}"
                    else:
                        taxonomy['scientific_name'] = f"{taxonomy['genus']} {taxonomy['species']}"
                elif taxonomy['genus']:
                    taxonomy['scientific_name'] = taxonomy['genus']
                else:
                    taxonomy['scientific_name'] = taxo_unite  # Fallback
            
            else:
                # Header simple sans taxonomie
                seq_id = header
                taxonomy = {'scientific_name': header}
            
            return seq_id, taxonomy
            
        except Exception as e:
            logger.warning(f"Failed to parse header '{header}': {e}")
            return header, {'scientific_name': header}
    
    def validate_sequence(self, sequence: str) -> Tuple[bool, List[str]]:
        """
        Valide une séquence d'ADN/RNA.
        
        Returns:
            tuple: (is_valid, list_of_issues)
        """
        issues = []
        sequence = sequence.upper().strip()
        
        # Vérifier la longueur minimale
        clean_seq = sequence.replace('-', '').replace('N', '')
        if len(clean_seq) < self.config_validation['min_length']:
            issues.append(f"Sequence too short: {len(clean_seq)} bp")
        
        # Vérifier les caractères valides (nucléotides + ambigus + gaps)
        valid_chars = set(self.config_validation['valid_chars'])    
        invalid_chars = set(sequence) - valid_chars
        if invalid_chars:
            issues.append(f"Invalid characters: {invalid_chars}")
        
        # Vérifier le contenu en N
        n_content = sequence.count('N') / len(sequence) if sequence else 0
        if n_content > self.config_validation['max_n_fraction']:
            issues.append(f"Too many ambiguous bases: {n_content:.1%}")
        
        # # Vérifier les homopolymères excessifs
        # for base in 'ACGT':
        #     if base * 20 in sequence:  # 20+ bases identiques consécutives
        #         issues.append(f"Excessive homopolymer: {base}")
        
        return len(issues) == 0, issues
    
    def parse_fasta_file(self, file_path: str) -> Iterator[Tuple[str, str, str]]:
        """
        Parse un fichier FASTA (gzippé ou non).
        
        Yields:
            tuple: (original_id, sequence, header)
        """
        def open_file(path):
            if path.endswith('.gz'):
                return gzip.open(path, 'rt', encoding='utf-8')
            return open(path, 'r', encoding='utf-8')
        
        try:
            with open_file(file_path) as f:
                header = None
                sequence_lines = []
                
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # Traiter la séquence précédente si elle existe
                        if header and sequence_lines:
                            sequence = ''.join(sequence_lines)
                            original_id, _ = self.parse_unite_header(header)
                            yield original_id, sequence, header
                        
                        # Nouvelle séquence
                        header = line
                        sequence_lines = []
                    
                    else:
                        # Ligne de séquence
                        sequence_lines.append(line)
                
                # Traiter la dernière séquence
                if header and sequence_lines:
                    sequence = ''.join(sequence_lines)
                    original_id, _ = self.parse_unite_header(header)
                    yield original_id, sequence, header
                    
        except Exception as e:
            logger.error(f"Error reading FASTA file {file_path}: {e}")
            raise
    
    def process_unite_database(self, input_file: str, output_file: str, 
                              source_name: str = "UNITE", source_version: str = "unknown") -> Dict:
        """
        Traite complètement une base UNITE : parsing + base de données + export.
        
        Returns:
            dict: Statistiques de traitement
        """
        start_time = time.time()
        logger.info(f"Starting UNITE database processing: {input_file}")
        
        # Créer le dossier de sortie
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Ajouter la source à la base de données
        source_id = self.db_manager.add_source(
            name=source_name,
            version=source_version,
            file_path=input_file,
            description=f"UNITE ITS database v{source_version}"
        )
        
        processed_sequences = []
        
        try:
            # Parser le fichier FASTA
            for original_id, sequence, header in self.parse_fasta_file(input_file):
                self.stats['total_sequences'] += 1
                
                try:
                    # Parser la taxonomie depuis le header
                    seq_id, taxonomy_info = self.parse_unite_header(header)
                    
                    # Valider la séquence
                    is_valid, issues = self.validate_sequence(sequence)
                    
                    if not is_valid:
                        self.stats['invalid_sequences'] += 1
                        logger.warning(f"Invalid sequence {original_id}: {', '.join(issues)}")
                        continue
                    
                    # Ajouter à la base de données
                    db_seq_id = self.db_manager.add_sequence(
                        source_id=source_id,
                        original_id=original_id,
                        sequence=sequence,
                        taxonomy_info=taxonomy_info
                    )
                    
                    # Log de la transformation
                    self.db_manager.log_transformation(
                        seq_id=db_seq_id,
                        step="format_unite",
                        tool_name="UNITEParser",
                        parameters={"source_file": input_file},
                        success=True
                    )
                    
                    # Stocker pour l'export FASTA
                    processed_sequences.append({
                        'db_seq_id': db_seq_id,
                        'original_id': original_id,
                        'sequence': sequence,
                        'scientific_name': taxonomy_info.get('scientific_name', 'Unknown')
                    })
                    
                    self.stats['processed_sequences'] += 1
                    
                    # Log progress
                    if self.stats['processed_sequences'] % 1000 == 0:
                        logger.info(f"Processed {self.stats['processed_sequences']} sequences...")
                
                except Exception as e:
                    self.stats['parsing_errors'] += 1
                    logger.error(f"Error processing sequence {original_id}: {e}")
                    continue
            
            # # Exporter le FASTA formaté
            # self.export_formatted_fasta(processed_sequences, output_file)
            
            # Mettre à jour les statistiques de la source
            with self.db_manager.get_connection() as conn:
                conn.execute("""
                    UPDATE sources SET record_count = ? WHERE source_id = ?
                """, (self.stats['processed_sequences'], source_id))
            
            end_time = time.time()
            self.stats['processing_time'] = end_time - start_time
            
            logger.info(f"UNITE processing completed in {self.stats['processing_time']:.2f}s")
            logger.info(f"Total: {self.stats['total_sequences']}, "
                       f"Processed: {self.stats['processed_sequences']}, "
                       f"Invalid: {self.stats['invalid_sequences']}, "
                       f"Errors: {self.stats['parsing_errors']}")
            
            return self.stats
            
        except Exception as e:
            logger.error(f"Critical error during UNITE processing: {e}")
            raise
    
    # def export_formatted_fasta(self, sequences: List[Dict], output_file: str) -> None:
    #     """Exporte les séquences formatées vers un fichier FASTA."""
        
    #     logger.info(f"Exporting {len(sequences)} sequences to {output_file}")
        
    #     try:
    #         with open(output_file, 'w') as f:
    #             for seq_data in sequences:
    #                 # Header formaté avec ID de base et nom scientifique
    #                 header = f">{seq_data['db_seq_id']}|{seq_data['original_id']}|{seq_data['scientific_name']}"
    #                 f.write(f"{header}\n")
                    
    #                 # Séquence avec retour à la ligne tous les 80 caractères
    #                 sequence = seq_data['sequence']
    #                 for i in range(0, len(sequence), 80):
    #                     f.write(f"{sequence[i:i+80]}\n")
            
    #         logger.info(f"Successfully exported formatted FASTA to {output_file}")
            
    #     except Exception as e:
    #         logger.error(f"Failed to export FASTA: {e}")
    #         raise


def main():
    """Point d'entrée principal pour Snakemake."""
    
    try:
        # Récupération des paramètres Snakemake
        input_db = snakemake.input[0]
        input_file = snakemake.input[1]
        config_validation = snakemake.params.sequence_validation
        output_file = snakemake.output[0]
        config = snakemake.config

        
        # Log des paramètres
        logger.info(f"Input file: {input_file}")
        logger.info(f"Input db: {input_db}")
        logger.info(f"Output file: {output_file}")
        
        # Initialisation du gestionnaire de base de données
        db_manager = get_database_manager(config)
        db_manager.create_schema()  # Créer le schéma si nécessaire
        
        # Traitement UNITE
        parser = UNITEParser(db_manager, config_validation)
        
        # Extraire version depuis le nom de fichier si possible
        version = "unknown"
        if "unite" in input_file.lower():
            # Exemple: sh_general_release_dynamic_29.11.2022.fasta.gz
            version_match = re.search(r'(\d{1,2}\.\d{1,2}\.\d{4})', input_file)
            if version_match:
                version = version_match.group(1)
        
        # Lancer le traitement
        stats = parser.process_unite_database(
            input_file=input_file,
            output_file=output_file,
            source_name="UNITE",
            source_version=version
        )
        
        # Validation de l'intégrité de la base
        integrity_report = db_manager.validate_integrity()
        if not integrity_report['is_valid']:
            logger.warning(f"Database integrity issues detected: {integrity_report}")
        
        # Optimiser la base de données
        db_manager.vacuum_database()
        
        # Rapport final
        logger.info("=== UNITE Processing Report ===")
        logger.info(f"Total sequences found: {stats['total_sequences']}")
        logger.info(f"Successfully processed: {stats['processed_sequences']}")
        logger.info(f"Invalid sequences: {stats['invalid_sequences']}")
        logger.info(f"Parsing errors: {stats['parsing_errors']}")
        logger.info(f"Processing time: {stats['processing_time']:.2f} seconds")
        logger.info(f"Database integrity: {'OK' if integrity_report['is_valid'] else 'ISSUES'}")
        
        # Créer un rapport JSON pour Snakemake
        import json
        report_data = {
            'source': 'UNITE',
            'version': version,
            'input_file': input_file,
            'output_file': output_file,
            'stats': stats,
            'integrity_report': integrity_report,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        report_file = output_file.replace('.fasta', '_report.json')
        with open(report_file, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        logger.info(f"Processing report saved to: {report_file}")
        
    except Exception as e:
        logger.error(f"UNITE processing failed: {e}")
        raise


if __name__ == "__main__":
    main()
