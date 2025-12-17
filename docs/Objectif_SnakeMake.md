# Cahier des charges – Pipeline **SnakeMake**

## Objectif

Mettre en place un pipeline bioinformatique complet et reproductible avec **SnakeMake**, permettant le prétraitement, le nettoyage et l’augmentation des données pour la constitution de jeux d’entraînement/test destinés à l’apprentissage automatique.

Le pipeline doit garantir :

* **Reproductibilité** : chaque étape est versionnée et paramétrable.
* **Traçabilité** : les transformations sont historisées dans une base de données optimisée.
* **Cohérence des identifiants** : utilisation d’ID uniques persistants, évitant les doublons.
* **Efficacité mémoire et CPU** : pas de duplication inutile de séquences, usage optimal du parallélisme.
* **Interopérabilité** : formats standardisés (CSV/FASTA/JSON) pour faciliter l’analyse et l’intégration.

---

## 1. Règles principales du pipeline

### 1.1. Récupération et formatage des bases de données

* **Sources** :

  * UNITE
  * BOLD
  * PLANiTS
  * CALeDNA

* **Tâches** :

  * Normaliser les métadonnées (taxonomie, séquence, identifiants, source).
  * Assigner des **ID uniques globaux** (UUID ou hash stable).
  * Vérifier l’intégrité et la complétude des données.
  * Eventuellement mettre en place des filtres de qualités : taille, %N

### 1.2. Harmonisation taxonomique

* Aligner toutes les bases de données contre le **Catalogue of Life (CoL)**.
* Résoudre les synonymes taxonomiques et les conflits de rang.
* Conserver une hiérarchie taxonomique standardisée (kingdom → species).

### 1.3. Extraction ITS

* Appliquer **ITSx** sur toutes les séquences afin d’isoler les régions ITS1 et ITS2 pertinentes.
* Annoter les résultats dans un format standardisé et adapter les nouveaux ID si besoin.

### 1.4. Fusion des bases de données

* Combiner les différentes bases et les séquences ITS extraites.
* Supprimer les doublons basés sur séquences + taxonomie.
* Adapter les métadonnées et garantir la cohérence des ID.

### 1.5. Alignement multiple (MSA)

* Réaliser un **alignement multiple optimisé** (par ex. MAFFT).
* Stocker les alignements.

### 1.6. Validation taxonomique

* Appliquer **SATIVA** sur les MSA pour identifier et corriger les erreurs d’assignation taxonomique.

### 1.7. Augmentation de données

* Appliquer **TRICHE** (data augmentation artificielle).
* Générer des séquences augmentées avec métadonnées cohérentes.

### 1.8. Génération des datasets d’apprentissage et de test

* Construire des **ensembles stratifiés** (train/test/validation).
* Garantir un équilibre taxonomique et éviter les fuites de données.
* Exporter en formats standards (FASTA/CSV/Parquet).

---

## 2. Contraintes techniques

* **Langage & Workflow** : SnakeMake, Python
* **Reproductibilité** :

  * Gestion des environnements avec **Conda/Mamba**.
  * Fichiers `environment.yaml` pour chaque rule.
* **Bases de données** : SQLite.
* **Scalabilité** : exécution locale ou cluster HPC (support SLURM).
* **Optimisation** : parallélisme SnakeMake (`--cores`), vectorisation si applicable.
* **Documentation & Logs** : génération automatique de rapports HTML (Snakemake `--report`).

---

## 3. Livrables

* `Snakefile` documenté.
* Dossiers `rules/`, `scripts/`, `config/` (paramètres YAML).
* Fichier `config.yaml` définissant chemins, versions DB, paramètres outils.
* Base de données centralisée contenant l’historique des transformations.
* Jeux de données finaux (`train/test/validation`).
* Rapport de pipeline (HTML + PDF) retraçant toutes les étapes.
