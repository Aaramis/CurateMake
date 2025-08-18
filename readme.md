
# Bio-pipeline

Pipeline bioinformatique pour le traitement et l'analyse de séquences ADN.

## 🏗️ Structure du projet

```
.
├── config/               # Fichiers de configuration
│   ├── config.yaml       # Configuration principale
│   └── .env.example     # Variables d'environnement (modifier en .env)
│
├── data/                 # Données brutes et traitées
│   ├── raw/             # Données brutes
│   └── processed/       # Données traitées
│
├── docs/                # Documentation complète
│   ├── overview.md      # Vue d'ensemble
│   ├── installation.md  # Guide d'installation
│   ├── usage.md         # Guide d'utilisation
│   ├── database.md      # Schéma de la base de données
│   ├── tools.md         # Documentation des outils
│   ├── dev.md           # Guide développeur
│   ├── contributing.md  # Comment contribuer
│   └── changelog.md     # Historique des versions
│
├── envs/                # Fichiers d'environnement Conda
│   ├── environment.yaml      # Environnement de base
│   └── dev-environment.yaml  # Environnement de développement
│
├── results/             # Résultats des analyses
│   ├── itsx/           # Résultats d'extraction ITS
│   ├── alignments/     # Alignements multiples
│   ├── sativa/         # Résultats de validation taxonomique
│   └── augmented/      # Données augmentées
│
├── scripts/             # Scripts personnalisés
│   ├── utils/          # Utilitaires
│   └── analysis/       # Scripts d'analyse
│
├── tests/               # Tests automatisés
│   ├── unit/           # Tests unitaires
│   └── integration/    # Tests d'intégration
│
├── .gitignore          # Fichiers ignorés par Git
├── mkdocs.yml          # Configuration de la documentation
├── readme.md           # Ce fichier
└── setup_env.sh        # Script de configuration de l'environnement
```

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