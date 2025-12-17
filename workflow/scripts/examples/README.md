# Database Profiling Examples

This directory contains examples demonstrating how to use the database profiling features.

## Quick Start

```python
from utils.db_manager import DatabaseManager

# Create database manager
db = DatabaseManager("data/db/your_database.db")

# Enable profiling
db.enable_profiling(True)
db.set_slow_query_threshold(1.0)  # seconds

# Your database operations...
db.create_schema()
db.add_source("MySource", "path/to/file.fasta")

# Generate report
report = db.generate_profiling_report("profiling_report.json")
print(f"Total operations: {report['total_operations']}")
print(f"Total time: {report['total_time']:.2f}s")
```

## Running the Examples

```bash
# Run all examples
python workflow/scripts/examples/db_profiling_example.py
```

## Features

### 1. Automatic Operation Tracking

All major database operations are automatically profiled:
- `create_schema()` - Schema creation
- `add_source()` - Source additions
- `add_taxonomy()` - Taxonomy additions
- `get_source_stats()` - Statistics queries
- `validate_integrity()` - Integrity checks
- `vacuum_database()` - Database optimization
- And many more...

### 2. Slow Query Detection

Set a threshold to automatically detect slow operations:

```python
db.set_slow_query_threshold(0.5)  # 0.5 seconds

# Operations slower than threshold are flagged
slow_queries = db.get_slow_queries(limit=10)
for query in slow_queries:
    print(f"{query['operation']}: {query['duration_seconds']:.2f}s")
```

### 3. Statistics and Reports

Get detailed statistics:

```python
# All operations
stats = db.get_profiling_stats()

# Specific operation
stats = db.get_profiling_stats("add_taxonomy")

# Full report
report = db.generate_profiling_report("report.json")
```

### 4. Database Size Metrics

Monitor database growth:

```python
metrics = db.get_database_size_metrics()
print(f"Database size: {metrics['size_mb']} MB")
print(f"Tables: {metrics['table_row_counts']}")
```

## Integration with Snakemake

Add profiling to your Snakemake scripts:

```python
# In your Snakemake script
from utils.db_manager import DatabaseManager

# Get database from config
db = DatabaseManager(snakemake.config["database"]["path"])

# Enable profiling
db.enable_profiling(True)

# Your operations...

# Generate report at the end
report_path = snakemake.log[0].replace(".log", "_profiling.json")
db.generate_profiling_report(report_path)
```

## Configuration

You can configure profiling behavior:

```python
# Enable/disable profiling
DatabaseManager.enable_profiling(True)  # or False

# Set slow query threshold
DatabaseManager.set_slow_query_threshold(1.0)  # seconds

# Reset profiling data
DatabaseManager.reset_profiling()
```

## Output Format

Profiling reports are saved as JSON with the following structure:

```json
{
  "generated_at": "2025-11-26T10:30:00",
  "total_operations": 150,
  "total_time": 45.23,
  "slow_query_threshold": 1.0,
  "statistics": {
    "add_taxonomy": {
      "count": 100,
      "total_time": 12.5,
      "avg_time": 0.125,
      "min_time": 0.05,
      "max_time": 0.8,
      "slow_count": 5,
      "slow_percentage": 5.0
    }
  },
  "slow_queries": [...],
  "top_operations_by_total_time": [...]
}
```

## Best Practices

1. **Enable profiling in development** to identify bottlenecks early
2. **Use conditional profiling in production** - enable only when debugging
3. **Set appropriate thresholds** - 1.0s for production, 0.1s for development
4. **Review reports regularly** to catch performance regressions
5. **Reset profiling data** between runs to avoid memory buildup

## Performance Impact

- **Profiling enabled**: ~1-2% overhead
- **Profiling disabled**: Near-zero overhead

The profiling system is designed to have minimal impact on performance.
