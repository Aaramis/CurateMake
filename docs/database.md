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
    version TEXT NOT NULL,
    url TEXT,
    date_download TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### 2. `sequences`

Contient les séquences d’ADN et leur état.

```sql
CREATE TABLE sequences (
    seq_id TEXT PRIMARY KEY,      -- UUID généré côté Python
    source_id INTEGER,
    original_id TEXT,
    sequence TEXT NOT NULL,
    length INTEGER,               -- calculé côté application ou via trigger
    md5 CHAR(32) UNIQUE,
    status TEXT,                  -- 'raw','itsx','aligned','augmented'
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (source_id) REFERENCES sources(source_id) ON DELETE CASCADE
);

CREATE INDEX idx_seq_md5 ON sequences(md5);
```

#### 3. `taxonomy`

Hiérarchie taxonomique normalisée.

```sql
CREATE TABLE taxonomy (
    tax_id INTEGER PRIMARY KEY AUTOINCREMENT,
    kingdom TEXT,
    phylum TEXT,
    class TEXT,
    "order" TEXT,
    family TEXT,
    genus TEXT,
    species TEXT,
    scientific_name TEXT UNIQUE
);
```

#### 4. `sequence_taxonomy`

Table de jonction entre séquences et taxonomie.

```sql
CREATE TABLE sequence_taxonomy (
    seq_id TEXT,
    tax_id INTEGER,
    method TEXT,
    confidence FLOAT,               -- entre 0 et 1, validation côté application
    PRIMARY KEY (seq_id, tax_id, method),
    FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE,
    FOREIGN KEY (tax_id) REFERENCES taxonomy(tax_id) ON DELETE CASCADE
);
```

#### 5. `transformations`

Historique des transformations appliquées aux séquences.

```sql
CREATE TABLE transformations (
    transform_id INTEGER PRIMARY KEY AUTOINCREMENT,
    seq_id TEXT,
    step TEXT NOT NULL,            -- ITSx, MSA, SATIVA, TRICHE
    parameters JSON,               -- JSON simple
    tool_version TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE
);

CREATE INDEX idx_transform_seq ON transformations(seq_id);
```

#### 6. `datasets`

Gestion des jeux de données (train/test/validation).

```sql
CREATE TABLE datasets (
    dataset_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### 7. `dataset_sequences`

Lien entre séquences et jeux de données.

```sql
CREATE TABLE dataset_sequences (
    dataset_id INTEGER,
    seq_id TEXT,
    PRIMARY KEY (dataset_id, seq_id),
    FOREIGN KEY (dataset_id) REFERENCES datasets(dataset_id) ON DELETE CASCADE,
    FOREIGN KEY (seq_id) REFERENCES sequences(seq_id) ON DELETE CASCADE
);
```

---

### Bonnes pratiques

1. **Normalisation** : pas de duplication des séquences.
2. **Traçabilité** : chaque transformation loggée avec paramètres et version de l’outil.
3. **Flexibilité** : ajout facile de nouvelles sources ou étapes.
4. **Scalabilité** : indices sur `md5`, `seq_id`, `tax_id` pour accélérer les jointures.
