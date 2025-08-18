## Base de données du pipeline

### Schéma relationnel

La base de données assure la traçabilité complète des transformations des données. Elle utilise PostgreSQL pour ses performances et fonctionnalités avancées.

### Tables principales

#### 1. `sources`

Stocke les métadonnées des différentes sources de données.

```sql
CREATE TABLE sources (
    source_id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,           -- UNITE, BOLD, PLANiTS...
    version TEXT NOT NULL,
    url TEXT,
    date_download TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);
```

#### 2. `sequences`

Contient les séquences d’ADN et leur état.

```sql
CREATE TABLE sequences (
    seq_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    source_id INT REFERENCES sources(source_id) ON DELETE CASCADE,
    original_id TEXT,
    sequence TEXT NOT NULL,
    length INT GENERATED ALWAYS AS (char_length(sequence)) STORED,
    md5 CHAR(32) UNIQUE,
    status TEXT CHECK (status IN ('raw','itsx','aligned','augmented')),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX idx_seq_md5 ON sequences(md5);
```

#### 3. `taxonomy`

Hiérarchie taxonomique normalisée.

```sql
CREATE TABLE taxonomy (
    tax_id SERIAL PRIMARY KEY,
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
    seq_id UUID REFERENCES sequences(seq_id) ON DELETE CASCADE,
    tax_id INT REFERENCES taxonomy(tax_id) ON DELETE CASCADE,
    method TEXT,
    confidence FLOAT CHECK (confidence BETWEEN 0 AND 1),
    PRIMARY KEY (seq_id, tax_id, method)
);
```

#### 5. `transformations`

Historique des transformations appliquées aux séquences.

```sql
CREATE TABLE transformations (
    transform_id SERIAL PRIMARY KEY,
    seq_id UUID REFERENCES sequences(seq_id) ON DELETE CASCADE,
    step TEXT NOT NULL,            -- ITSx, MSA, SATIVA, TRICHE
    parameters JSONB,
    tool_version TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX idx_transform_seq ON transformations(seq_id);
```

#### 6. `datasets`

Gestion des jeux de données (train/test/validation).

```sql
CREATE TABLE datasets (
    dataset_id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### 7. `dataset_sequences`

Lien entre séquences et jeux de données.

```sql
CREATE TABLE dataset_sequences (
    dataset_id INT REFERENCES datasets(dataset_id) ON DELETE CASCADE,
    seq_id UUID REFERENCES sequences(seq_id) ON DELETE CASCADE,
    PRIMARY KEY (dataset_id, seq_id)
);
```

---

### Bonnes pratiques

1. **Normalisation** : pas de duplication des séquences.
2. **Traçabilité** : chaque transformation loggée avec paramètres et version de l’outil.
3. **Flexibilité** : ajout facile de nouvelles sources ou étapes.
4. **Scalabilité** : indices sur `md5`, `seq_id`, `tax_id` pour accélérer les jointures.
