# Guide de contribution

Merci de votre intérêt pour ce projet ! Nous apprécions toutes les contributions, qu'il s'agisse de signaler des bogues, de proposer des fonctionnalités ou de soumettre des corrections de code.

## Table des matières

- [Comment contribuer](#comment-contribuer)
- [Signalement de bogues](#signalement-de-bogues)
- [Propositions de fonctionnalités](#propositions-de-fonctionnalités)
- [Soumission de modifications](#soumission-de-modifications)
- [Conventions de code](#conventions-de-code)
- [Tests](#tests)
- [Revue de code](#revue-de-code)
- [Questions et support](#questions-et-support)

## Comment contribuer

1. **Signaler un problème**
   - Vérifiez d'abord si le problème n'a pas déjà été signalé
   - Ouvrez une issue avec un titre clair et descriptif
   - Suivez le modèle d'issue fourni

2. **Proposer une amélioration**
   - Décrivez clairement la fonctionnalité ou l'amélioration
   - Expliquez pourquoi elle serait utile
   - Donnez un exemple d'utilisation si possible

3. **Soumettre des modifications**
   - Suivez le [workflow de développement](dev.md#workflow-de-développement)
   - Assurez-vous que les tests passent
   - Mettez à jour la documentation si nécessaire

## Signalement de bogues

Pour signaler un bug, ouvrez une issue en suivant ces étapes :

1. Utilisez le modèle "Bug Report"
2. Donnez un titre clair et descriptif
3. Décrivez le problème de manière détaillée
4. Indiquez les étapes pour reproduire le bug
5. Précisez le comportement attendu
6. Ajoutez des captures d'écran si nécessaire
7. Indiquez votre configuration (OS, version de Python, etc.)

## Propositions de fonctionnalités

Pour proposer une nouvelle fonctionnalité :

1. Vérifiez d'abord si une demande similaire existe déjà
2. Utilisez le modèle "Feature Request"
3. Expliquez clairement la fonctionnalité souhaitée
4. Décrivez les cas d'utilisation
5. Proposez une ébauche d'implémentation si possible

## Soumission de modifications

### Prérequis

- Python 3.8+
- Git
- Un compte GitHub
- Connaissance de base de Snakemake

### Processus de contribution

1. **Fork** le dépôt
2. Clonez votre fork localement :
   ```bash
   git clone https://github.com/votre-utilisateur/bio-pipeline.git
   cd bio-pipeline
   ```
3. Créez une branche pour votre fonctionnalité :
   ```bash
   git checkout -b feature/nom-de-la-fonctionnalite
   ```
4. Installez les dépendances de développement :
   ```bash
   pip install -e ".[dev]"
   ```
5. Effectuez vos modifications
6. Vérifiez le style de code :
   ```bash
   black .
   flake8
   ```
7. Exécutez les tests :
   ```bash
   pytest
   ```
8. Committez vos modifications avec un message clair :
   ```bash
   git commit -m "Ajout: Nouvelle fonctionnalité pour résoudre #123"
   ```
9. Poussez vers votre fork :
   ```bash
   git push origin feature/nom-de-la-fonctionnalite
   ```
10. Créez une Pull Request (PR) vers la branche `main`

## Conventions de code

### Style de code

- Suivez [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Utilisez [Black](https://black.readthedocs.io/) pour le formatage
- Vérifiez avec [Flake8](https://flake8.pycqa.org/)
- Documentez les fonctions avec des docstrings (format Google)

### Structure des commits

- Utilisez des messages de commit clairs et descriptifs
- Suivez le format [Conventional Commits](https://www.conventionalcommits.org/)
- Faites référence aux numéros d'issue quand c'est pertinent

Exemple :
```
feat: ajouter la fonctionnalité X

Ajoute une nouvelle fonctionnalité qui permet de faire X.
Résout le problème #123
```

Types de commits :
- `feat`: Nouvelle fonctionnalité
- `fix`: Correction de bug
- `docs`: Modification de la documentation
- `style`: Formatage, point-virgule manquant, etc.
- `refactor`: Modification du code qui ne corrige pas un bug ni n'ajoute de fonctionnalité
- `perf`: Modification du code qui améliore les performances
- `test`: Ajout ou modification de tests
- `chore`: Mise à jour de tâches de construction, configuration, etc.

## Tests

### Exécution des tests

```bash
# Tous les tests
pytest

# Tests avec couverture
pytest --cov=.

# Un fichier de test spécifique
pytest tests/test_utils.py

# Un test spécifique
pytest tests/test_utils.py::test_ma_fonction
```

### Écriture de tests

- Un fichier de test par module : `test_<nom_du_module>.py`
- Une fonction de test par fonctionnalité
- Utilisez des fixtures pour les données de test courantes
- Testez les cas limites et les erreurs

## Revue de code

Toutes les contributions passent par un processus de revue de code. Voici comment cela fonctionne :

1. Créez une Pull Request (PR)
2. Assurez-vous que tous les tests passent
3. Un mainteneur examinera votre code
4. Des modifications peuvent être demandées
5. Une fois approuvée, votre PR sera fusionnée

### Critères d'acceptation

- Le code est bien testé
- La documentation est à jour
- Les tests passent
- Le code suit les conventions du projet
- La fonctionnalité est utile pour d'autres utilisateurs

## Questions et support

Pour toute question ou problème, vous pouvez :

1. Consulter la [documentation](index.md)
2. Ouvrir une issue sur GitHub
3. Contacter les mainteneurs via email (voir le fichier AUTHORS.md)

## Reconnaissance

Tous les contributeurs sont remerciés pour leur temps et leurs efforts. Votre nom sera ajouté à la liste des contributeurs dans le fichier CONTRIBUTORS.md.
