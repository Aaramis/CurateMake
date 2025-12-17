#!/usr/bin/env python3
"""
Example script demonstrating database profiling usage.

This script shows how to:
1. Enable profiling
2. Execute database operations
3. Generate profiling reports
4. Analyze slow queries
"""

import sys
from pathlib import Path

# Add the scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.db_manager import DatabaseManager
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def example_basic_profiling():
    """Example: Basic profiling usage."""

    print("=" * 60)
    print("Example 1: Basic Database Profiling")
    print("=" * 60)

    # Create database manager with a test database
    db = DatabaseManager("data/db/test_profiling.db")

    # Enable profiling
    db.enable_profiling(True)
    db.set_slow_query_threshold(0.1)  # Queries slower than 0.1s are "slow"

    # Create schema (will be profiled)
    db.create_schema()

    # Add some test data (will be profiled)
    source_id = db.add_source("TestSource", "/path/to/test.fasta")

    # Add taxonomy (will be profiled)
    tax_id = db.add_taxonomy(
        kingdom="Fungi",
        phylum="Ascomycota",
        genus="Aspergillus",
        species="niger"
    )

    # Get statistics (will be profiled)
    stats = db.get_source_stats()

    # Get profiling statistics
    print("\n--- Profiling Statistics ---")
    profiling_stats = db.get_profiling_stats()
    for operation, data in profiling_stats.items():
        print(f"\n{operation}:")
        print(f"  Count: {data['count']}")
        print(f"  Avg time: {data['avg_time']:.4f}s")
        print(f"  Total time: {data['total_time']:.4f}s")
        print(f"  Min/Max: {data['min_time']:.4f}s / {data['max_time']:.4f}s")
        if data['slow_count'] > 0:
            print(f"  ⚠️  Slow queries: {data['slow_count']} ({data['slow_percentage']:.1f}%)")


def example_profiling_report():
    """Example: Generate comprehensive profiling report."""

    print("\n" + "=" * 60)
    print("Example 2: Generate Profiling Report")
    print("=" * 60)

    db = DatabaseManager("data/db/test_profiling.db")
    db.enable_profiling(True)

    # Run some operations
    db.get_source_stats()
    db.get_backbone_stats()

    # Generate report
    report_path = "data/profiling_report.json"
    report = db.generate_profiling_report(report_path)

    print(f"\n✓ Report saved to: {report_path}")
    print(f"\nReport Summary:")
    print(f"  Total operations: {report['total_operations']}")
    print(f"  Total time: {report['total_time']:.4f}s")
    print(f"  Slow query threshold: {report['slow_query_threshold']}s")

    if report.get('top_operations_by_total_time'):
        print("\n--- Top 3 Operations by Total Time ---")
        for i, op in enumerate(report['top_operations_by_total_time'][:3], 1):
            print(f"{i}. {op['operation']}: {op['total_time']:.4f}s ({op['count']} calls)")


def example_slow_queries():
    """Example: Identify slow queries."""

    print("\n" + "=" * 60)
    print("Example 3: Identify Slow Queries")
    print("=" * 60)

    db = DatabaseManager("data/db/test_profiling.db")
    db.enable_profiling(True)
    db.set_slow_query_threshold(0.05)  # More aggressive threshold

    # Run some operations that might be slow
    db.validate_integrity()
    db.vacuum_database()

    # Get slow queries
    slow_queries = db.get_slow_queries(limit=5)

    if slow_queries:
        print(f"\n⚠️  Found {len(slow_queries)} slow queries:")
        for i, query in enumerate(slow_queries, 1):
            print(f"\n{i}. {query['operation']}")
            print(f"   Duration: {query['duration_seconds']:.4f}s")
            print(f"   Timestamp: {query['timestamp']}")
            if query.get('metadata'):
                print(f"   Metadata: {query['metadata']}")
    else:
        print("\n✓ No slow queries detected!")


def example_database_size_metrics():
    """Example: Get database size metrics."""

    print("\n" + "=" * 60)
    print("Example 4: Database Size Metrics")
    print("=" * 60)

    db = DatabaseManager("data/db/test_profiling.db")

    metrics = db.get_database_size_metrics()

    if metrics['exists']:
        print(f"\nDatabase: {metrics['path']}")
        print(f"Size: {metrics['size_mb']} MB ({metrics['size_bytes']:,} bytes)")

        if 'table_row_counts' in metrics:
            print("\nTable row counts:")
            for table, count in metrics['table_row_counts'].items():
                print(f"  {table}: {count:,} rows")
    else:
        print("\n⚠️  Database does not exist yet")


def example_profiling_in_production():
    """Example: Profiling pattern for production use."""

    print("\n" + "=" * 60)
    print("Example 5: Production Profiling Pattern")
    print("=" * 60)

    # In production, you might want to:
    # 1. Enable profiling at the start
    # 2. Run your pipeline
    # 3. Generate report at the end

    db = DatabaseManager("data/db/production.db")

    # Enable profiling with appropriate threshold
    db.enable_profiling(True)
    db.set_slow_query_threshold(1.0)  # Production threshold

    try:
        # Your pipeline operations here
        db.create_schema()
        # ... more operations ...

        print("\n✓ Pipeline completed successfully")

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")

    finally:
        # Always generate report at the end
        report_path = "logs/profiling_report.json"
        report = db.generate_profiling_report(report_path)

        print(f"\n✓ Profiling report saved to: {report_path}")
        print(f"  Operations profiled: {report['total_operations']}")
        print(f"  Total execution time: {report['total_time']:.2f}s")

        # Check for slow queries
        slow_count = sum(
            data['slow_count']
            for data in report['statistics'].values()
        )
        if slow_count > 0:
            print(f"  ⚠️  {slow_count} slow queries detected - review report for optimization")

        # Optionally reset profiling data for next run
        db.reset_profiling()


def example_conditional_profiling():
    """Example: Enable profiling only when needed."""

    print("\n" + "=" * 60)
    print("Example 6: Conditional Profiling")
    print("=" * 60)

    # You can enable/disable profiling dynamically
    # Useful for debugging specific issues

    db = DatabaseManager("data/db/test.db")

    # Disable by default (minimal overhead)
    db.enable_profiling(False)

    # ... normal operations ...

    # Enable when you need to debug
    print("\nEnabling profiling for debugging...")
    db.enable_profiling(True)

    # Operations you want to profile
    db.get_source_stats()

    # Disable again
    db.enable_profiling(False)

    # Get stats only for profiled operations
    stats = db.get_profiling_stats()
    print(f"\nProfiled {len(stats)} operations")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Database Profiling Examples")
    print("=" * 60)

    # Run all examples
    try:
        example_basic_profiling()
        example_profiling_report()
        example_slow_queries()
        example_database_size_metrics()
        example_profiling_in_production()
        example_conditional_profiling()

        print("\n" + "=" * 60)
        print("✓ All examples completed successfully!")
        print("=" * 60)

    except Exception as e:
        logger.error(f"Example failed: {e}", exc_info=True)
        sys.exit(1)
