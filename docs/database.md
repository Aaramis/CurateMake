# Database Configuration

The workflow targets two backends:

| Backend | When to use | Configuration keys |
| ------- | ----------- | ------------------ |
| SQLite  | Local prototyping / CI | `database.type: sqlite`, path defined in `database.sqlite.path/name` |
| PostgreSQL | Production scale, concurrent writers | `database.type: postgresql`, connection info under `database.postgresql` |

## 1. Environment variables

Secrets such as host, port, user, and password are referenced via `{{ env.VAR }}` inside `config/config.yaml`. Populate them in `config/.env`:

```env
POSTGRES_HOST=db.internal
POSTGRES_PORT=5433
POSTGRES_USER=its_pipeline
POSTGRES_PASSWORD=***
POSTGRES_DB=its_curate
```

Snakemake loads `.env` at startup (`workflow/Snakefile`). Missing variables will raise a clear error.

## 2. Config keys

```yaml
database:
  type: "postgresql"        # or "sqlite"
  sqlite:
    path: "data/db"
    name: "its.db"
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
  backup_retention_days: 30
  vacuum_threshold_mb: 100
```

- `schema` lets the SATIVA integration and ingestion scripts qualify tables automatically.
- `pool_size`/`max_overflow` control SQLAlchemy’s connection pooling when running multiple jobs per process.

## 3. Creating/initializing the schema

```bash
# Create tables on the selected backend
snakemake --use-conda create_database
```

This rule relies on `workflow/scripts/create_database.py`. For PostgreSQL, ensure the user has `CREATE` privileges within the target schema.

## 4. Runtime behavior

- Formatting rules (`workflow/rules/format_databases.smk`) use `ancient(get_db_path())` to serialize writes when using SQLite, and allow higher parallelism via `db_lock` resources when using PostgreSQL.
- ITSx/MSA/SATIVA scripts use `workflow/scripts/utils/db_manager.py` to abstract away the backend. Switching between SQLite and PostgreSQL requires no code changes beyond the config.
- SATIVA integration writes warning summaries and optional backups according to `sativa_integration.database.*` and `sativa_integration.reporting.*` blocks.

## 5. Backups & hygiene

- Enable automatic backups before SATIVA integration by toggling `sativa_integration.database.backup_before_integration`.
- Manual backup example (SQLite):
  ```bash
  sqlite3 data/db/its.db ".backup data/db/its.db.$(date +%Y%m%d)"
  ```
- Manual backup example (PostgreSQL):
  ```bash
  pg_dump -h $POSTGRES_HOST -U $POSTGRES_USER $POSTGRES_DB > backup.sql
  ```

Keep `data/db` (SQLite) or the remote schema under versioned retention policies if regulatory compliance is required.
