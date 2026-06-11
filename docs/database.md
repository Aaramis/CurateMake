# Database Configuration

The workflow stores every record in a normalised relational schema and supports
two interchangeable backends. Switching backend requires only a configuration
change — no code edits.

| Backend | When to use | Configuration keys |
| ------- | ----------- | ------------------ |
| SQLite | Local prototyping / CI / small datasets | `database.type: sqlite`, path in `database.sqlite.path` + `database.sqlite.name` |
| PostgreSQL | Production scale, concurrent writers | `database.type: postgresql`, connection info under `database.postgresql` |

## 1. Schema

Five tables (defined in `workflow/scripts/utils/db_models.py`):

| Table | Primary key | Role |
| ----- | ----------- | ---- |
| `sources` | `source_id` | One row per input dataset (UNITE, BOLD, PLANiTS, CALeDNA, PERSONAL) |
| `sequences` | `seq_id` | Sequence text, `length`, `percentage_n`, `barcode`, `barcode_status`, optional `parent_seq_id` |
| `taxonomy` | `tax_id` | Unique lineages `kingdom → species` (uniqueness enforced across all ranks) |
| `seq_tax` | (`seq_id`, `tax_id`, `backbone`) | Sequence↔taxonomy links with provenance and validation `status` |
| `sativa_warnings` | `warning_id` | Modifications flagged by SATIVA integration for manual review |

The composite key of `seq_tax` is what makes the **two-stage curation**
auditable. The `backbone` column records where an assignment came from:

| `backbone` value | Meaning |
| ---------------- | ------- |
| `Raw` | Original label as ingested from the source database |
| `CoL` | Label reconciled against the Catalogue of Life |
| `Sativa` | Label validated/corrected by SATIVA |
| `GBIF` | Reserved for an optional GBIF backbone |

Because the three revisions coexist per sequence, downstream reports (notably
`sativa_taxonomy_report.html`) can present a four-way comparison of the original,
CoL, and SATIVA assignments.

## 2. Environment variables

Secrets — host, port, user, password, database — are referenced via
`{{ env.VAR }}` inside `config/config.yaml` and resolved from a `.env` file in
the **project root** (loaded by the `Snakefile` at start-up). Copy the template
and edit it:

```bash
cp config/.env.example .env
```

```env
POSTGRES_HOST=db.internal
POSTGRES_PORT=5433
POSTGRES_USER=its_pipeline
POSTGRES_PASSWORD=***
POSTGRES_DB=its_curate
```

Missing variables raise a clear error at start-up.

## 3. Config keys

```yaml
database:
  type: "postgresql"        # or "sqlite"
  sqlite:
    path: "data/db"
    name: "its_mini_database.db"
  postgresql:
    host: "{{ env.POSTGRES_HOST }}"
    port: "{{ env.POSTGRES_PORT }}"
    user: "{{ env.POSTGRES_USER }}"
    password: "{{ env.POSTGRES_PASSWORD }}"
    database: "{{ env.POSTGRES_DB }}"
    schema: "its_curate"
    pool_size: 10
    max_overflow: 20
    echo: false
  writer_retry:
    max_attempts: 100        # retries for transient errors on distributed DBs
  backup_retention_days: 30
  vacuum_threshold_mb: 100
```

- `schema` lets ingestion and SATIVA integration qualify tables automatically
  (a `SET search_path` is issued on each PostgreSQL connection).
- `pool_size` / `max_overflow` control SQLAlchemy connection pooling for
  parallel jobs.
- `writer_retry.max_attempts` makes writes resilient on distributed backends
  (e.g. YugabyteDB), which return transient errors under contention.

## 4. Creating / initialising the schema

```bash
snakemake --use-conda --cores 1 create_database
```

This rule runs `workflow/scripts/create_database.py` and is idempotent
(`CREATE TABLE IF NOT EXISTS`). For PostgreSQL, the user must hold `CREATE`
privileges in the target schema.

## 5. Concurrency model

- Formatting rules acquire a `db_lock` resource. For SQLite the lock level is `1`
  (a single writer, since SQLite serialises writes); for PostgreSQL it is
  effectively unlimited (MVCC allows concurrent writers).
- All ingestion, ITSx, MSA, and SATIVA scripts go through
  `workflow/scripts/utils/db_manager.py`, which abstracts the backend.

## 6. Backups & hygiene

- SQLite snapshot:
  ```bash
  sqlite3 data/db/its_mini_database.db ".backup data/db/backup_$(date +%Y%m%d).db"
  ```
- PostgreSQL dump:
  ```bash
  pg_dump -h "$POSTGRES_HOST" -U "$POSTGRES_USER" "$POSTGRES_DB" > backup.sql
  ```

Keep `data/db/` (SQLite) or the remote schema under a retention policy if your
project requires it.
