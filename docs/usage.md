# Guide d'utilisation du pipeline

Ce document explique comment exécuter le pipeline bioinformatique, configurer les paramètres et interpréter les résultats.

## Table des matières

- [Démarrage rapide](#démarrage-rapide)
- [Configuration](#configuration)
- [Exécution du pipeline](#exécution-du-pipeline)
- [Résultats](#résultats)
- [Surveillance et débogage](#surveillance-et-débogage)
- [Exemples de cas d'utilisation](#exemples-de-cas-dutilisation)

## Démarrage rapide

### 1. Préparation

```bash
# Activer l'environnement Conda
conda activate bio-pipeline

# Se placer dans le répertoire du projet
cd /chemin/vers/bio-pipeline

# Copier le fichier de configuration par défaut
cp config/config.default.yaml config/config.yaml
```

### 2. Exécution complète

```bash
# Lancer le pipeline complet (mode sec pour tester)
snakemake --use-conda -n

# Lancer réellement le pipeline (utilise 8 cœurs)
snakemake --use-conda --cores 8
```

## Configuration

### Fichier de configuration principal

Le fichier `config/config.yaml` contient tous les paramètres modifiables :

```yaml
# Paramètres généraux
project_name: "bio-pipeline"
output_dir: "results"

# Base de données
database:
  host: "localhost"
  port: 5432
  name: "bio_pipeline"
  user: "pipeline_user"
  password: "votre_mot_de_passe"

# Paramètres des règles spécifiques
itsx:
  cpu: 4
  minlen: 50

mafft:
  algorithm: "--auto"
  threads: 8

sativa:
  min_support: 0.8
  min_identity: 0.98

triche:
  mode: "standard"
  num_augment: 5
```

### Fichier des échantillons

Créez un fichier `samples.tsv` dans le dossier `config/` :

```tsv
sample_id	path
sample1	data/raw/sample1.fastq
sample2	data/raw/sample2.fastq
```

## Exécution du pipeline

### Exécution complète

```bash
# Exécuter tout le pipeline
snakemake --use-conda --cores 8
```

### Exécution partielle

```bash
# N'exécuter que jusqu'à l'extraction ITS
snakemake --use-conda --cores 8 itsx

# N'exécuter que la validation taxonomique
snakemake --use-conda --cores 8 sativa
```

### Options d'exécution courantes

```bash
# Exécuter en mode sec (affiche ce qui serait fait)
snakemake -n

# Forcer la réexécution de certaines règles
snakemake --forceall --touch

# Nettoyer les fichiers temporaires
snakemake --delete-all-output

# Exécuter sur un cluster SLURM
snakemake --cluster "sbatch -t {resources.time_min} --mem={resources.mem_mb}" -j 10
```

## Résultats

### Structure des répertoires

```
results/
├── raw/               # Données brutes
├── processed/         # Données prétraitées
│   ├── trimmed/       # Données nettoyées
│   └── normalized/    # Données normalisées
├── itsx/              # Résultats d'extraction ITS
│   ├── its1/          # Régions ITS1
│   ├── 5_8s/          # Régions 5.8S
│   └── its2/          # Régions ITS2
├── taxonomy/          # Résultats taxonomiques
│   ├── assigned/      # Assignations taxonomiques
│   └── validated/     # Résultats de validation
├── augmented/         # Données augmentées
│   ├── train/         # Ensemble d'entraînement
│   └── test/          # Ensemble de test
└── reports/           # Rapports et visualisations
    ├── stats/         # Statistiques
    └── plots/         # Graphiques
```

### Fichiers de sortie principaux

1. **ITSx**
   - `results/itsx/itsx_summary.txt` : Résumé des régions ITS extraites
   - `results/itsx/its1.fasta` : Séquences ITS1 extraites

2. **MAFFT**
   - `results/alignments/aligned.fasta` : Alignement multiple
   - `results/alignments/alignment.log` : Fichier de log

3. **SATIVA**
   - `results/sativa/misclassified.txt` : Séquences potentiellement mal classées
   - `results/sativa/summary.pdf` : Rapport de validation

4. **TRICHE**
   - `results/augmented/train/augmented.fasta` : Données d'entraînement augmentées
   - `results/augmented/test/test.fasta` : Données de test

## Surveillance et débogage

### Journaux d'exécution

Les journaux sont stockés dans le répertoire `.snakemake/log/` :

```bash
# Suivre les journaux en temps réel
tail -f .snakemake/log/*.snakemake.log
```

### Règles en échec

En cas d'échec, utilisez `--reason` pour plus d'informations :

```bash
snakemake --use-conda --cores 8 --reason
```

### Nettoyage

```bash
# Supprimer tous les fichiers générés
snakemake --delete-all-output

# Nettoyer uniquement les fichiers temporaires
snakemake cleanup
```

## Exemples de cas d'utilisation

### 1. Traitement d'un nouvel échantillon

1. Ajoutez les chemins des nouveaux fichiers FASTQ dans `config/samples.tsv`
2. Exécutez le pipeline complet :
   ```bash
   snakemake --use-conda --cores 8
   ```

### 2. Régénération des données augmentées

```bash
# Supprimer les anciennes données augmentées
rm -rf results/augmented/*

# Régénérer uniquement les données augmentées
snakemake --use-conda --cores 8 augmented
```

### 3. Exécution sur un cluster

Créez un profil de cluster dans `profiles/cluster/` et exécutez :

```bash
snakemake --profile profiles/cluster
```

## Assistance

Pour toute question ou problème, consultez la [documentation complète](overview.md) ou ouvrez une issue sur le dépôt du projet.
