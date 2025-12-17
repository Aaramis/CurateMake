# Pipeline Bioinformatique - Documentation

Bienvenue dans la documentation du pipeline bioinformatique basé sur Snakemake. Ce pipeline est conçu pour le prétraitement, le nettoyage et l'augmentation des données pour la constitution de jeux d'entraînement/test destinés à l'apprentissage automatique.

## Fonctionnalités principales

- **Reproductibilité** : Chaque étape est versionnée et paramétrable
- **Traçabilité** : Historique complet des transformations dans une base de données PostgreSQL
- **Efficacité** : Optimisé pour les performances avec gestion intelligente de la mémoire et du CPU
- **Interopérabilité** : Formats standardisés (CSV/FASTA/JSON) pour une intégration aisée

## Structure du projet

```
Source/
│── Snakefile
│── config/
│   └── config.yaml
│── rules/
│   ├── download.smk
│   ├── taxonomy.smk
│   ├── itsx.smk
│   ├── msa.smk
│   ├── sativa.smk
│   ├── triche.smk
│   ├── report.smk
│   └── datasets.smk
│── scripts/   # Scripts Python/R
│── envs/      # Fichiers .yaml pour conda
│── logs/      # Journaux d'exécution
│── db/        # Schéma et scripts de base de données
│── docs/      # Documentation
└── examples/  # Exemples de jeux de données
```

## Premiers pas

1. [Installation](installation.md)
2. [Guide d'utilisation rapide](usage.md)
3. [Configuration de la base de données](database.md)

<!-- ## Ressources

- [Documentation complète](overview.md)
- [Guide du développeur](dev.md)
- [Comment contribuer](contributing.md)
- [Journal des modifications](changelog.md) -->
