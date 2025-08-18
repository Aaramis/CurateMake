# Outils du pipeline

Ce document décrit les principaux outils bioinformatiques utilisés dans le pipeline, leurs versions recommandées et leurs paramètres de configuration.

## Table des matières

- [ITSx](#itsx)
- [MAFFT](#mafft)
- [SATIVA](#sativa)
- [TRICHE](#triche)
- [Autres dépendances](#autres-dépendances)

## ITSx

### Description
ITSx permet d'identifier et d'extraire les régions ITS (Internal Transcribed Spacer) des séquences d'ADN ribosomique.

### Version recommandée
- Version : 1.1.3
- Site web : [https://microbiology.se/software/itsx/](https://microbiology.se/software/itsx/)

### Paramètres principaux

```yaml
itsx:
  # Options de base
  cpu: 4
  save_regions: all  # ITS1, 5.8S, ITS2
  complement: no     # Pas de séquences complémentaires
  detailed_results: yes
  
  # Filtrage
  minlen: 50         # Longueur minimale des séquences
  allow_reorder: yes  # Réordonnancement des séquences
  
  # Fichiers de sortie
  outdir: "results/itsx"
  basename: "itsx_results"
```

### Exemple de commande

```bash
ITSx \
  -i input.fasta \
  -o results/itsx/itsx_results \
  --cpu 4 \
  --save_regions all \
  --complement no \
  --detailed_results T \
  --minlen 50
```

---

## MAFFT

### Description
Outil d'alignement multiple de séquences.

### Version recommandée
- Version : 7.520
- Site web : [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)

### Paramètres principaux

```yaml
mafft:
  algorithm: "--auto"  # Choix automatique de la méthode
  threads: 8
  adjustdirection: yes  # Détection automatique du brin
  maxiterate: 1000      # Nombre maximal d'itérations
  op: 1.53             # Gap opening penalty
  ep: 0.123            # Offset (comme recommandé)
  reorder: yes         # Réorganisation des séquences
```

### Exemple de commande

```bash
mafft \
  --auto \
  --thread 8 \
  --adjustdirection \
  --maxiterate 1000 \
  --op 1.53 \
  --ep 0.123 \
  --reorder \
  input.fasta > aligned.fasta
```

---

## SATIVA

### Description
Outil de validation taxonomique des séquences basé sur la méthode de référence.

### Version recommandée
- Version : 1.10
- Site web : [https://github.com/amkozlov/sativa](https://github.com/amkozlov/sativa)

### Paramètres principaux

```yaml
sativa:
  # Configuration de base
  threads: 8
  tmpdir: "tmp/sativa"
  outdir: "results/sativa"
  
  # Paramètres d'analyse
  min_support: 0.8      # Support minimal pour les erreurs
  min_identity: 0.98    # Identité minimale
  max_placements: 10    # Nombre maximal de placements
  
  # Filtrage
  min_seq_len: 100      # Longueur minimale des séquences
  max_ambig: 0.1        # Pourcentage maximal de bases ambiguës
```

### Exemple de commande

```bash
python sativa.py \
  -s input.fasta \
  -t reference.tax \
  -x ZXJB \
  -T 8 \
  --min-support 0.8 \
  --min-identity 0.98 \
  -o results/sativa
```

---

## TRICHE

### Description
Outil d'augmentation de données pour les séquences biologiques.

### Version recommandée
- Version : 2.0.1
- Dépôt : [Lien vers le dépôt TRICHE]

### Paramètres principaux

```yaml
triche:
  # Configuration de base
  mode: "standard"  # standard, conservative, aggressive
  num_augment: 5    # Nombre d'augmentations par séquence
  
  # Paramètres de mutation
  mutation_rate: 0.01
  indel_rate: 0.001
  max_indel_length: 3
  
  # Filtrage
  min_quality: 20
  min_coverage: 10
```

### Exemple de commande

```bash
triche \
  -i input.fasta \
  -o augmented_sequences \
  --mode standard \
  --num-augment 5 \
  --mutation-rate 0.01 \
  --indel-rate 0.001
```

---

## Autres dépendances

### Python
- Python ≥ 3.8
- Bibliothèques principales :
  - Biopython
  - pandas
  - NumPy
  - SQLAlchemy
  - psycopg2-binary (pour PostgreSQL)

### R (si nécessaire)
- R ≥ 4.0
- Paquets recommandés :
  - phyloseq
  - dada2
  - DECIPHER

### Environnement
- Conda/Mamba pour la gestion des environnements
- PostgreSQL ≥ 13 pour la base de données
- Snakemake ≥ 7.0 pour l'orchestration du workflow
