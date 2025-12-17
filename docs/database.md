## Base de données du pipeline (SQLite)

### Schéma relationnel

La base de données assure la traçabilité complète des transformations des données. Elle utilise SQLite pour sa simplicité et portabilité.

### Tables principales

#### 1. `sources`

Stocke les métadonnées des différentes sources de données.

```sql
CREATE TABLE sources (
    source_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,           -- UNITE, BOLD, PLANiTS...
    file_path TEXT,
    url TEXT,
    date_download TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

Voici une version **améliorée**, plus professionnelle mais toujours simple et claire, de la documentation de ta table `sequences`. J'ai reformulé pour améliorer la lisibilité, clarifié certains points, corrigé quelques fautes, et mis en forme pour que même un développeur fatigué puisse suivre.

---

#### 2. Table `sequences`

La table `sequences` contient les données de séquences nucléotidiques (ADN), leur origine, leur état dans le pipeline, ainsi que les relations hiérarchiques éventuelles entre elles (ex. réassignations, dérivations).

Elle est alimentée automatiquement depuis des sources externes, puis enrichie au fil des étapes du pipeline (ITSx, MSA, SATIVA...).

---

##### Définition SQL

```sql
CREATE TABLE IF NOT EXISTS sequences (
    seq_id INTEGER PRIMARY KEY AUTOINCREMENT,         -- Identifiant unique
    source_id INTEGER NOT NULL,                       -- Référence à la table 'sources'
    original_id TEXT NOT NULL,                        -- Identifiant original dans la base source
    sequence TEXT NOT NULL,                           -- Séquence nucléotidique brute
    length INTEGER NOT NULL CHECK (length > 0),       -- Longueur de la séquence
    percentage_n REAL DEFAULT 0.0 CHECK (
        percentage_n >= 0 AND percentage_n <= 100
    ),                                                -- Pourcentage de 'N' (ambiguïté)

    barcode TEXT DEFAULT 'ITS' CHECK (
        barcode IN ('ITS', 'ITS1', 'ITS2', 'LSU', 'SSU', '5.8S')
    ),                                                -- Région ciblée selon ITSx
    old_barcode TEXT DEFAULT NULL,                    -- Ancien barcode, si modifié par réassignation
    barcode_status TEXT DEFAULT 'Pending' CHECK (
        barcode_status IN (
            'Pending',        -- Ajoutée, en attente de traitement
            'Duplicate',      -- Séquence et taxonomie identiques à une autre
            'ITSx',           -- Séquence générée par ITSx
            'Processed',      -- Séquence traitée par ITSx
            'no_detection'    -- ITSx n’a rien détecté
        )
    ),

    msa_id INTEGER DEFAULT NULL,                      -- Référence au MSA associé (si applicable)
    parent_seq_id INTEGER DEFAULT NULL,               -- Référence à la séquence d’origine
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,   -- Date d’insertion

    FOREIGN KEY (source_id) REFERENCES sources(source_id) ON DELETE CASCADE,
    FOREIGN KEY (parent_seq_id) REFERENCES sequences(seq_id) ON DELETE SET NULL
);
```

---

##### Pipeline d’utilisation

Voici le cycle de vie standard d’une séquence dans la base :

1. **Ajout depuis une source**
   Les séquences sont importées automatiquement depuis les bases sources.
   Si le barcode (ITS, ITS1, ITS2, etc.) est précisé dans la source, il est conservé ; sinon, la valeur par défaut est `'ITS'`.
   Tous les champs techniques sont initialisés automatiquement (`length`, `percentage_n`, etc.).

2. **Détection des doublons**
   Une première passe est effectuée pour marquer comme `Duplicate` toutes les séquences **identiques en séquence ET en taxonomie**.

3. **Traitement par ITSx**

   * Les séquences avec `barcode_status = 'Pending'` sont traitées avec ITSx.
   * Si une région ITS est détectée, la séquence est mise à jour avec `barcode_status = 'Processed'`.
   * Si aucune région n’est détectée, le statut devient `no_detection`.
   * Les nouvelles séquences générées par ITSx (ex. ITS1, ITS2) sont insérées avec :

     * `barcode_status = 'ITSx'`
     * `parent_seq_id` = ID de la séquence d’origine

4. **Post-ITSx : Nouvelle passe de détection des doublons**
   Une fois les nouvelles séquences ajoutées, une seconde détection des doublons est réalisée pour intégrer les ITSx.

5. **MSA / Alignement multiple (à venir)**
   La logique d’intégration dans les MSA est en cours de définition.
   Chaque séquence alignée pourra être liée à un `msa_id`.

---

##### Notes

* **parent_seq_id** permet de tracer les dérivations de séquences (ex. séquences tronquées par ITSx, ou réassignées par MSA/SATIVA).


Oh joy, another half-documented table wandering the taxonomic wastelands like a lost intern. Don’t worry, I’ll give it a nice shower, professional clothes, and a LinkedIn profile.

---

#### 3. Table `taxonomy`

La table `taxonomy` contient les informations hiérarchiques de classification taxonomique, normalisées pour l'ensemble des séquences de la base.
Elle sert de référentiel structuré pour relier les séquences à une taxonomie stable, traçable et dédupliquée.

---

##### Définition SQL

```sql
CREATE TABLE IF NOT EXISTS taxonomy (
    tax_id INTEGER PRIMARY KEY AUTOINCREMENT,       -- Identifiant unique de la taxonomie
    kingdom TEXT,                                   -- Règne
    phylum TEXT,                                    -- Embranchement
    class_ TEXT,                                    -- Classe (nommé class_ pour éviter les conflits SQL)
    order_ TEXT,                                    -- Ordre (order_ pour les mêmes raisons)
    family TEXT,                                    -- Famille
    genus TEXT,                                     -- Genre
    species TEXT,                                   -- Espèce
    parent_tax_id INTEGER,                          -- Référence vers la taxonomie originale si modifiée
    
    FOREIGN KEY (parent_tax_id) REFERENCES taxonomy(tax_id) ON DELETE SET NULL
);
```

---

##### Pipeline d’utilisation

1. **Import initial des taxonomies**

   * Lors de l’intégration des séquences depuis une base source, chaque taxonomie rencontrée est ajoutée à la table `taxonomy`.

2. **Déduplication**

   * Une passe de nettoyage est effectuée pour supprimer les doublons exacts.

3. **Alignement sur Catalogue of Life (CoL)**

   * Chaque taxonomie est comparée à la référence du **Catalogue of Life** (ou autre source officielle).
   * Si une taxonomie est corrigée ou alignée, une **nouvelle entrée** est ajoutée avec :

     * Les valeurs corrigées,
     * Et un `parent_tax_id` pointant vers la taxonomie originale (non alignée).
   * Cette approche permet de conserver la traçabilité complète de la taxonomie initiale et de ses modifications.

---

##### Notes

* Les colonnes `class_` et `order_` utilisent des suffixes pour éviter des collisions avec des mots-clés SQL réservés (`CLASS`, `ORDER`).
* Si une taxonomie est modifiée par un outil ou une validation manuelle, **une nouvelle ligne est créée** au lieu de modifier l’existante, ce qui garantit l’historique des décisions.

Ah, la table `seq_tax`. Le carrefour où séquences et taxonomies s’embrassent, se trahissent, se recroisent, puis finissent toutes par s’accuser mutuellement d’être “Pending”.
Bon, on va lui redonner une documentation digne de ce nom — claire, logique, et surtout sans cette confusion infernale entre `Raw`, `CoL`, `Correct`, `Pending`, et le chaos existentiel.

---

#### 4. Table `seq_tax`

Table de jonction entre les **séquences** (`sequences`) et les **taxonomies** (`taxonomy`).
Elle trace l’évolution du lien entre une séquence donnée et sa classification taxonomique au cours du pipeline (source brute → alignement CoL → validation finale).

---

#####  Définition SQL

```sql
CREATE TABLE IF NOT EXISTS seq_tax (
    seq_id INTEGER NOT NULL,
    tax_id INTEGER NOT NULL,
    backbone TEXT DEFAULT 'Raw' CHECK (backbone IN ('Raw', 'CoL', 'GBIF')),
    status TEXT DEFAULT 'Pending' CHECK (status IN ('Pending', 'Check', 'Correct', 'Failed')),
    
    PRIMARY KEY (seq_id, tax_id, backbone),
    FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE,
    FOREIGN KEY (tax_id) REFERENCES taxonomy(tax_id) ON DELETE CASCADE
);
```

---

##### Rôle et logique du pipeline

**1. Ajout initial**

Lors de l’intégration des sources :

* Chaque séquence (`seq_id`) est liée à une taxonomie (`tax_id`) issue des données brutes.
* Une entrée est créée dans `seq_tax` avec :

  * `backbone = 'Raw'`
  * `status = 'Pending'`

**2. Déduplication des taxonomies**

Lors du nettoyage de la table `taxonomy` :

* Si deux taxonomies sont identiques, la relation `seq_tax` est mise à jour pour pointer vers le `tax_id` conservé.
* Les doublons (`tax_id` supprimés) sont redirigés pour maintenir l’intégrité référentielle.

**3. Alignement contre le Catalogue of Life (CoL)**

À ce stade, les taxonomies avec :

```
backbone = 'Raw' AND status = 'Pending'
```

sont envoyées à l’API CoL pour vérification.

Trois cas de figure sont possibles :

| Cas                                | Résultat CoL                                   | Action sur `seq_tax`                                                                                                                    | Description                                                                            |
| ---------------------------------- | ---------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| **Validation directe**           | Même classification que la taxonomie existante | `status = 'Correct'` sur la relation `Raw`                                                                                              | Rien ne change, la taxonomie est validée telle quelle.                                 |
| **Nouvelle classification CoL** | Nouvelle taxonomie trouvée (nouveau `tax_id`)  | 1️⃣ Ajouter une nouvelle ligne : `backbone = 'CoL'`, `status = 'Correct'` <br>2️⃣ Mettre à jour l’ancienne (`Raw`) → `status = 'Check'` | L’ancienne taxonomie reste pour traçabilité, mais la version CoL devient la référence. |
| **Échec du matching**            | Aucun résultat CoL                             | `status = 'Failed'` sur la ligne `Raw`                                                                                                  | Le taxon n’a pas pu être validé.                                                       |


**4. État final attendu**

| status      | backbone | Signification                                                        |
| ----------- | -------- | -------------------------------------------------------------------- |
| **Correct** | `Raw`    | Taxonomie d’origine validée par CoL (pas de modification nécessaire) |
| **Check**   | `Raw`    | Ancienne taxonomie remplacée par une version CoL                     |
| **Correct** | `CoL`    | Nouvelle taxonomie alignée avec CoL                                  |
| **Failed**  | `Raw`    | Taxonomie sans correspondance trouvée dans CoL                       |

---

**Problèmes actuels observés (à corriger)**

| Problème | Cause probable | Comportement actuel                                         |                                                                                |                                                                         |
| -------- | -------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------ | ----------------------------------------------------------------------- |
| `Pending | CoL`           | Mauvaise mise à jour du statut après succès sans changement | La logique du script ne met pas `status=Correct` si la taxonomie est identique |                                                                         |
| `Pending | Raw`           | Taxons non traités ou rollback partiel en cas d’erreur      | Le traitement CoL n’a pas passé la ligne en “Failed”                           |                                                                         |
| `Correct | Raw`+`Correct  | CoL`sur le même`seq_id`                                     | Mauvaise gestion des statuts en cas de match inchangé                          | L’ancienne et la nouvelle relations sont toutes deux marquées “Correct” |


##### Résumé visuel

```
RAW (Pending) ──┬──> RAW (Correct)          → même taxonomie validée
                 ├──> RAW (Check) + CoL (Correct) → nouvelle taxonomie trouvée
                 └──> RAW (Failed)          → pas de correspondance trouvée
```
