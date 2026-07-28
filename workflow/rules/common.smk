"""
Common helper functions for Snakemake rules.
"""

def get_db_path():
    """
    Get database path based on configured backend.

    Returns:
        str: Path to SQLite database file, or sentinel file for PostgreSQL
    """
    db_type = config["database"].get("type", "sqlite")

    if db_type == "sqlite":
        db_config = config["database"]["sqlite"]
        return os.path.join(db_config["path"], db_config["name"])
    else:
        # PostgreSQL doesn't create a local file
        # Return a sentinel file to track database initialization
        return "results/db/.postgresql_initialized"
