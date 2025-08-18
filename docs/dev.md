# Guide du développeur

Ce document fournit des informations pour les développeurs souhaitant contribuer au pipeline ou l'étendre avec de nouvelles fonctionnalités.

## Table des matières

- [Architecture du code](#architecture-du-code)
- [Ajout d'une nouvelle règle](#ajout-dune-nouvelle-règle)
- [Tests](#tests)
- [Bonnes pratiques](#bonnes-pratiques)
- [Workflow de développement](#workflow-de-développement)
- [Dépannage](#dépannage)

## Architecture du code

### Structure des dossiers

```
.
├── Snakefile           # Point d'entrée principal
├── config/             # Fichiers de configuration
│   └── config.yaml     # Configuration principale
├── rules/              # Règles SnakeMake
│   ├── download.smk    # Téléchargement des données
│   ├── taxonomy.smk    # Traitement taxonomique
│   ├── itsx.smk        # Extraction ITS
│   ├── msa.smk         # Alignement multiple
│   ├── sativa.smk      # Validation taxonomique
│   ├── triche.smk      # Augmentation de données
│   ├── report.smk      # Génération de rapports
│   └── datasets.smk    # Création des jeux de données
├── scripts/            # Scripts personnalisés
│   ├── utils/          # Utilitaires partagés
│   └── analysis/       # Scripts d'analyse
├── envs/               # Environnements Conda
│   ├── base.yaml       # Environnement de base
│   └── tools/          # Environnements par outil
├── tests/              # Tests automatisés
│   ├── unit/           # Tests unitaires
│   └── integration/    # Tests d'intégration
└── docs/               # Documentation
```

### Workflow principal

1. **Téléchargement des données** : Récupération des données sources
2. **Prétraitement** : Nettoyage et préparation des données
3. **Analyse** : Traitement principal (ITSx, alignement, etc.)
4. **Validation** : Contrôle qualité et validation
5. **Post-traitement** : Génération des rapports et jeux de données finaux

## Ajout d'une nouvelle règle

### 1. Créer un nouveau fichier de règles

Créez un nouveau fichier dans le dossier `rules/`, par exemple `rules/nouvelle_regle.smk` :

```python
# rules/nouvelle_regle.smk

rule nouvelle_regle:
    input:
        "{sample}.fasta"
    output:
        "results/nouveau/{sample}.out"
    params:
        param1=config["nouvelle_regle"]["param1"],
        param2=config["nouvelle_regle"]["param2"]
    conda:
        "envs/tools/nouvel_outil.yaml"
    shell:
        "nouvel_outil --input {input} --output {output} --param1 {params.param1}"
```

### 2. Mettre à jour la configuration

Ajoutez les paramètres par défaut dans `config/config.default.yaml` :

```yaml
# Configuration pour la nouvelle règle
nouvelle_regle:
  param1: valeur1
  param2: valeur2
```

### 3. Mettre à jour l'environnement Conda

Créez un fichier `envs/tools/nouvel_outil.yaml` :

```yaml
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - nouvel_outil=1.0.0
  - python=3.8
```

### 4. Mettre à jour le Snakefile principal

Ajoutez l'import du nouveau fichier de règles dans `Snakefile` :

```python
# Inclure les règles
include: "rules/nouvelle_regle.smk"
```

### 5. Tester la nouvelle règle

```bash
# Tester la règle en mode sec
snakemake --use-conda -n nouvelle_regle

# Exécuter la règle
snakemake --use-conda --cores 1 nouvelle_regle
```

## Tests

### Tests unitaires

Les tests unitaires vérifient le bon fonctionnement des fonctions individuelles :

```bash
# Exécuter tous les tests unitaires
pytest tests/unit/

# Exécuter un test spécifique
pytest tests/unit/test_utils.py -v
```

### Tests d'intégration

Les tests d'intégration vérifient le bon fonctionnement des règles SnakeMake :

```bash
# Exécuter les tests d'intégration
pytest tests/integration/ -v

# Exécuter avec Snakemake
snakemake --use-conda -j 1 --use-singularity --configfile tests/integration/config.yaml
```

### Couverture de code

Pour générer un rapport de couverture :

```bash
pytest --cov=scripts tests/
coverage html
```

## Bonnes pratiques

### Style de code

- Suivez le [PEP 8](https://www.python.org/dev/peps/pep-0008/) pour le code Python
- Utilisez des noms de variables explicites
- Documentez les fonctions avec des docstrings
- Commentez les parties complexes du code

### Gestion des dépendances

- Utilisez toujours des versions spécifiques dans les fichiers YAML des environnements
- Mettez à jour régulièrement les dépendances pour corriger les failles de sécurité
- Testez toujours après une mise à jour de dépendances

### Documentation

- Mettez à jour la documentation quand vous modifiez le code
- Ajoutez des exemples d'utilisation
- Documentez les paramètres des fonctions et des règles

## Workflow de développement

### 1. Créer une branche

```bash
git checkout -b feature/nouvelle-fonctionnalite
```

### 2. Développer et tester

- Écrivez le code
- Ajoutez des tests
- Vérifiez que les tests passent
- Mettez à jour la documentation

### 3. Soumettre une demande de fusion

1. Poussez vos modifications :
   ```bash
   git push origin feature/nouvelle-fonctionnalite
   ```
2. Créez une Pull Request sur GitHub
3. Attendez la revue de code
4. Effectuez les modifications demandées si nécessaire
5. Une fois approuvé, fusionnez la branche

## Dépannage

### Problèmes courants

#### 1. Échec de création d'environnement Conda

```bash
# Nettoyer les caches Conda
conda clean --all

# Réessayer avec Mamba
mamba env update -f envs/base.yaml
```

#### 2. Règles qui échouent silencieusement

Activez le mode debug :

```bash
snakemake --use-conda --cores 1 --verbose --printshellcmds --debug-dag
```

#### 3. Problèmes de verrouillage

Si Snakemake se bloque, supprimez les verrous :

```bash
rm -f .snakemake/locks/*
```

### Débogage avec PDB

Pour déboguer un script Python :

```python
import pdb; pdb.set_trace()
```

## Performance

### Profilage

Pour profiler une règle :

```bash
# Activer le profilage
export SNAKEMAKE_PROFILE=profiles/performance/

# Exécuter avec profilage
snakemake --profile profiles/performance/
```

### Optimisation

- Utilisez des conteneurs pour les outils lourds
- Limitez l'utilisation de la mémoire avec `resources.mem_mb`
- Utilisez des fichiers temporaires pour les données volumineuses

## Contribution

Vos contributions sont les bienvenues ! Voici comment contribuer :

1. Forkez le dépôt
2. Créez une branche pour votre fonctionnalité
3. Committez vos modifications
4. Poussez vers votre fork
5. Créez une Pull Request

Consultez le fichier `CONTRIBUTING.md` pour plus de détails.
