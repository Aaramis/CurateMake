# MVP - Alignement Hiérarchique Multi-BarCodes

Ce projet implémente un système d'alignement hiérarchique de séquences génétiques pour l'analyse taxonomique, avec support de plusieurs marqueurs génétiques (barcodes).

## Structure du Projet

```
1.mvp_MSA/
├── output/                # Dossier de sortie (créé automatiquement)
│   ├── logs/              # Fichiers de log
│   └── barcode_*/         # Dossiers par barcode (ITS, ITS1, ITS2...)
│       ├── alignments_*/    # Alignements par niveau taxonomique
│       ├── hmm/           # Profils HMM générés
│       └── inter_barcode/ # Résultats de la classification croisée
├── data/                 
│   └── plants.csv        # Fichier d'entrée des séquences
├── config.py             # Configuration du système
├── logger_config.py      # Configuration des logs
├── utils.py              # Fonctions utilitaires
├── external_tools.py     # Intégration des outils externes
├── processor.py          # Logique principale d'alignement
├── main.py               # Point d'entrée du programme
└── requirements.txt      # Dépendances Python
```

## Fonctionnalités Principales

- **Traitement en deux phases** :
  - Phase 1 : Alignement intra-barcode
  - Phase 2 : Classification croisée inter-barcodes
- Support de multiples barcodes (ITS, ITS1, ITS2, etc.)
- Traitement hiérarchique par niveau taxonomique
- Parallélisation multi-niveaux pour des performances optimales
- Gestion robuste des erreurs et journalisation détaillée

## Configuration

Les paramètres principaux sont définis dans `config.py` :

```python
# Paramètres du workflow
BARCODES = ["ITS", "ITS1", "ITS2"]  # Barcodes à traiter
MIN_SEQUENCES_THRESHOLD = 20         # Nombre minimum de séquences par groupe
MAX_LOSS_PERCENTAGE = 20.0          # Pourcentage maximum de perte autorisé
taxonomic_levels = ["genus", "family", "order", "class", "phylum"]

# Paramètres de parallélisation
MAX_BARCODE_WORKERS = len(BARCODES)  # Nombre de barcodes traités en parallèle
MAX_GENUS_WORKERS = 45              # Nombre de genres traités en parallèle par barcode
```

## Structure des Résultats

Pour chaque barcode, une structure de dossiers est créée :

```
output/
├── ITS/
│   ├── alignments_family/ # Alignements par niveau taxonomique
│   ├── alignments_genus/ # Alignements par niveau taxonomique
│   ├── alignments_order/ # Alignements par niveau taxonomique
│   ├── alignments_species/ # Alignements par niveau taxonomique
│   ├── hmm/              # Profils HMM générés
│   └── logs/             # Logs détaillés
├── ITS1/
│   ├── alignments_genus/ # Alignements par niveau taxonomique
│   ├── alignments_genus/ # Alignements par niveau taxonomique
│   ├── alignments_order/ # Alignements par niveau taxonomique
│   ├── alignments_species/ # Alignements par niveau taxonomique
│   ├── hmm/              # Profils HMM générés
│   └── logs/             # Logs détaillés
├── ITS2/
│   ├── alignments_genus/ # Alignements par niveau taxonomique
│   ├── alignments_genus/ # Alignements par niveau taxonomique
│   ├── alignments_order/ # Alignements par niveau taxonomique
│   ├── alignments_species/ # Alignements par niveau taxonomique
│   ├── hmm/              # Profils HMM générés
│   └── logs/             # Logs détailléscd ..
└── main/                 # Fichiers de log détaillés
    ├── failures.logs     # Logs des échecs
    └── success.logs      # Logs des succès
```

## Logique du Pipeline

### Phase 1 : Traitement Intra-Barcode
1. Pour chaque barcode :
   - Regroupement des séquences par niveau taxonomique
   - Filtrage des groupes avec trop peu de séquences
   - Alignement multiple avec MAFFT
   - Nettoyage avec trimAl
   - Génération de profils HMM
   - Classification des séquences non traitées contre la base de profile hmm
   - Rajout des séquences classées aux MSA initiaux

### Phase 2 : Classification Croisée Inter-Barcodes
   - Récupération des séquences non classées et des profils hmm des autres barcodes
   - Classification des séquences non classées contre les profils hmm des autres barcodes
   - Rajout des séquences classées aux MSA initiaux
   - Sauvegarde des séquences non classées


## Utilisation

Créer un fichier csv temporaire avec ces colonnes
seq_id,sequence,kingdom,phylum,class,order,family,genus,species,barcode

logique : 

Dans la database il faut récupérer toutes les séquences et les taxonomies associées qui correspondent à ces critères :

* Récupération Des taxonomies Correct + celle qui n'ont pas fonctionné avec CoL pour espérer pouvoir les récupérer avec Sativa
* Après ITSx on a identifié toutes les séquences doublons il faut donc éviter de les prendre pour ne pas influencer les MSA

    Seq_Tax:
        backbone = 'CoL' and (status = 'Correct' or status = 'check') 
    Sequences:
        barcode_status != 'Duplicate' 


1. Préparer les données dans `data/data.csv`
2. Configurer les paramètres dans `config.py` si nécessaire
3. Exécuter le script :

```bash
python main.py
```

## Dépendances

- MAFFT (pour les alignements multiples)
- trimAl (pour le nettoyage des alignements)
- HMMER (pour la création et l'utilisation des profils HMM)

Les dépendances Python sont listées dans `requirements.txt` et peuvent être installées via :

```bash
pip install -r requirements.txt
```

## Notes

- Les outils externes (MAFFT, trimAl, HMMER) doivent être installés et accessibles dans le PATH
- Les logs détaillés sont disponibles dans le dossier `output/logs/`
- Les performances peuvent être ajustées via les paramètres de parallélisation dans `config.py`
