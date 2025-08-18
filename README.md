
# Bio-pipeline

Pipeline bioinformatique pour le traitement et l'analyse de séquences ADN.

## 🏗️ Structure du projet

Cette structure suit les [meilleures pratiques recommandées par Snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility) pour assurer la reproductibilité et la maintenabilité du workflow.

```
.
├── .git/               # Dossier de versionnement Git
├── .github/            # Fichiers de configuration GitHub
│   └── workflows/      # GitHub Actions workflows
│
├── workflow/           # Code source du workflow Snakemake
│   ├── rules/         # Règles modulaires du workflow
│   │   ├── download.smk
│   │   ├── itsx.smk
│   │   ├── msa.smk
│   │   └── sativa.smk
│   │
│   ├── envs/          # Environnements Conda par outil
│   │   ├── itsx.yaml
│   │   ├── mafft.yaml
│   │   └── sativa.yaml
│   │
│   ├── scripts/       # Scripts personnalisés
│   │   └── utils/     # Utilitaires
│   │       └── db_schema.py
│   │
│   ├── notebooks/     # Notebooks Jupyter pour l'analyse exploratoire
│   ├── report/        # Fichiers pour la génération de rapports
│   └── Snakefile      # Point d'entrée du workflow
│
├── config/            # Fichiers de configuration
│   ├── config.yaml    # Configuration principale
│   └── .env.example  # Variables d'environnement (à copier en .env)
│
├── data/              # Données brutes (versionnés si petits, sinon .gitignorés)
│   ├── raw/          # Données brutes originales
│   └── processed/    # Données prétraitées
│
├── results/           # Résultats des analyses (générés par le workflow)
│   ├── itsx/         # Résultats d'extraction ITS
│   ├── alignments/   # Alignements multiples
│   ├── sativa/       # Résultats de validation taxonomique
│   └── augmented/    # Données augmentées
│
├── resources/         # Fichiers de référence et ressources statiques
├── tests/             # Tests automatisés
│   ├── unit/         # Tests unitaires
│   └── integration/  # Tests d'intégration
│
├── docs/              # Documentation complète
│   ├── overview.md
│   ├── installation.md
│   ├── usage.md
│   ├── database.md
│   ├── tools.md
│   ├── dev.md
│   ├── contributing.md
│   └── changelog.md
│
├── .gitignore        # Fichiers ignorés par Git
├── LICENSE.md        # Licence du projet
├── README.md         # Ce fichier
└── setup_env.sh      # Script de configuration de l'environnement
```

### Explications de la structure

- **workflow/** : Contient tout le code exécutable du workflow
  - `rules/` : Règles modulaires organisées par fonctionnalité
  - `envs/` : Environnements Conda spécifiques à chaque outil
  - `scripts/` : Scripts personnalisés réutilisables
  - `notebooks/` : Analyses exploratoires et prototypes
  - `report/` : Modèles et ressources pour les rapports

- **config/** : Fichiers de configuration
  - `config.yaml` : Configuration principale du workflow
  - `.env.example` : Modèle pour les variables d'environnement

- **data/** : Données d'entrée (versionnées si petites)
  - `raw/` : Données brutes originales
  - `processed/` : Données prétraitées

- **results/** : Résultats générés par le workflow (non versionnés)

- **resources/** : Fichiers de référence et ressources statiques

## Installation rapide

1. **Cloner le dépôt**
   ```bash
   git clone [URL_DU_DEPOT]
   cd bio-pipeline
   ```

2. **Configurer l'environnement**
   ```bash
   # Rendre le script exécutable
   chmod +x setup_env.sh
   
   # Lancer la configuration
   ./setup_env.sh
   ```

3. **Configurer les variables d'environnement**
   ```bash
   cp config/.env.example config/.env
   # Éditer le fichier .env avec vos paramètres
   ```

4. **Activer l'environnement**
   ```bash
   conda activate bio-pipeline
   # ou pour le développement
   conda activate bio-pipeline-dev
   ```

## Documentation

La documentation complète est disponible dans le dossier `docs/`. Pour la consulter en local :

```bash
# Installer les dépendances de documentation
pip install mkdocs mkdocs-material

# Lancer le serveur de documentation
mkdocs serve
```

Puis ouvrez http://127.0.0.1:8000 dans votre navigateur.

## Développement avec Git

### Workflow de développement standard

1. **Mettre à jour votre branche locale**
   ```bash
   git pull origin main
   ```

2. **Créer une nouvelle branche pour votre fonctionnalité**
   ```bash
   git checkout -b feature/nom-de-la-fonctionnalite
   ```

3. **Faire vos modifications et les valider**

```bash
# Voir les fichiers modifiés
git status

# Ajouter les fichiers au staging
git add .  # Tous les fichiers
# ou
git add chemin/vers/fichier  # Fichier spécifique

# Valider les changements avec un message clair
git commit -m "Description claire des modifications"
```

4. **Pousser vos changements**

```bash
# Premier push d'une nouvelle branche
git push -u origin feature/nom-de-la-fonctionnalite

# Pour les pushes suivants sur la même branche
git push
```

### Gestion des branches

- `main` : Branche principale, toujours stable
- `develop` : Branche d'intégration pour les fonctionnalités
- `feature/*` : Nouvelles fonctionnalités
- `bugfix/*` : Corrections de bugs
- `hotfix/*` : Corrections critiques pour la production

### .gitignore

Le fichier `.gitignore` est déjà configuré pour ignorer :
- Fichiers sensibles (.env, mots de passe, clés API)
- Fichiers temporaires et de cache
- Fichiers de données volumineux
- Fichiers spécifiques aux éditeurs de code

## Documentation **MkDocs** :

```bash
pip install mkdocs mkdocs-material
mkdocs new .
mkdocs serve
```
Optionnal

```bash
mkdocs gh-deploy
```