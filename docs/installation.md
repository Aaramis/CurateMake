# Installation du pipeline

Ce guide vous explique comment installer et configurer l'environnement nécessaire pour exécuter le pipeline bioinformatique.

## Prérequis

- Système d'exploitation : Linux
- Gestionnaire de paquets : Conda ou Mamba (recommandé)
- Espace disque : Au moins 100 Go d'espace libre (pour les bases de données et résultats)
- Mémoire : 16 Go RAM minimum (32+ Go recommandé)
- Processeur : 4 cœurs minimum (8+ recommandés)

## 1. Installation de Mamba (recommandé)

Mamba est une alternative plus rapide à Conda pour la résolution des dépendances :

```bash
# Installation de Mamba dans un environnement Conda existant
conda install -n base -c conda-forge mamba

# Ou installation complète de Miniforge avec Mamba
# Télécharger depuis https://github.com/conda-forge/miniforge#mambaforge
# Puis exécuter le script d'installation
```

## 2. Configuration de l'environnement

### 2.1. Clonage du dépôt

```bash
git clone https://github.com/votre-utilisateur/bio-pipeline.git
cd bio-pipeline
```

### 2.2. Création de l'environnement Conda

Le pipeline utilise des environnements Conda isolés pour chaque règle. L'environnement principal contient Snakemake et les dépendances de base :

```bash
# Création de l'environnement principal
mamba env create -f envs/environment.yaml

# Activation de l'environnement
conda activate bio-pipeline
```

## 3. Configuration de la base de données

### 3.1. Installation de PostgreSQL

```bash
# Sur Ubuntu/Debian
sudo apt-get update
sudo apt-get install postgresql postgresql-contrib

# Sur CentOS/RHEL
sudo yum install postgresql-server postgresql-contrib
sudo postgresql-setup initdb
sudo systemctl start postgresql
sudo systemctl enable postgresql
```

### 3.2. Configuration de la base de données

```bash
# Se connecter à PostgreSQL
sudo -u postgres psql

# Créer un utilisateur et une base de données
CREATE DATABASE bio_pipeline;
CREATE USER pipeline_user WITH PASSWORD 'votre_mot_de_passe';
GRANT ALL PRIVILEGES ON DATABASE bio_pipeline TO pipeline_user;
\q
```

### 3.3. Initialisation du schéma

```bash
# Créer les tables et les relations
psql -U pipeline_user -d bio_pipeline -f db/init_db.sql
```

## 4. Installation des dépendances externes

### 4.1. Outils bioinformatiques

Le pipeline gère automatiquement l'installation des outils via Conda. Les versions recommandées sont :

- ITSx : 1.1.3
- MAFFT : 7.520
- SATIVA : 1.10
- TRICHE : 2.0.1

### 4.2. Dépendances Python

Toutes les dépendances Python sont listées dans `envs/environment.yaml` et seront installées automatiquement.

## 5. Configuration du pipeline

### 5.1. Fichier de configuration

Copiez le fichier de configuration par défaut et modifiez-le selon vos besoins :

```bash
cp config/config.default.yaml config/config.yaml
```

### 5.2. Paramètres importants à vérifier

Dans `config/config.yaml`, vérifiez notamment :

- Chemins d'accès aux données
- Paramètres de la base de données
- Nombre de cœurs à utiliser
- Paramètres spécifiques à chaque outil

## 6. Vérification de l'installation

Pour vérifier que tout est correctement installé :

```bash
# Vérifier que Snakemake est disponible
snakemake --version

# Tester la configuration (sans exécution réelle)
snakemake -n --use-conda
```

## Dépannage

### Problèmes courants

1. **Erreurs de permissions PostgreSQL**
   - Vérifiez que l'utilisateur a bien les droits sur la base de données
   - Consultez les logs : `sudo tail -f /var/log/postgresql/postgresql-*.log`

2. **Problèmes avec Conda**
   - Mettez à jour Conda : `conda update -n base -c defaults conda`
   - Essayez avec Mamba : `mamba install -c conda-forge mamba`

3. **Espace disque insuffisant**
   - Le pipeline nécessite beaucoup d'espace pour les fichiers temporaires
   - Utilisez `--directory` pour spécifier un répertoire temporaire sur une partition avec assez d'espace

## Mise à jour du pipeline

Pour mettre à jour vers une version plus récente :

```bash
git pull origin main
mamba env update -f envs/environment.yaml
```

## Support

Pour toute question ou problème, veuillez ouvrir une issue sur le [dépôt GitHub](https://github.com/votre-utilisateur/bio-pipeline/issues).
