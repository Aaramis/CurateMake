"""
Database profiling utilities for performance monitoring and optimization.

This module provides tools to:
- Profile slow SQL queries with EXPLAIN QUERY PLAN
- Measure database operation execution times
- Track database connections and size metrics
- Generate profiling reports for optimization

Usage:
    from db_profiler import profile_db_operation, DBProfiler

    # Using decorator
    @profile_db_operation("insert_sequences")
    def insert_sequences_batch(conn, data):
        # ... code ...
        pass

    # Using context manager
    with DBProfiler.profile("complex_query"):
        conn.execute(query)

    # Generate report
    DBProfiler.generate_report("profiling_report.json")
"""

import sqlite3
import time
import functools
import json
import logging
from pathlib import Path
from typing import Optional, Callable, Any, Dict, List
from contextlib import contextmanager
from datetime import datetime
from collections import defaultdict

logger = logging.getLogger(__name__)


class DBProfiler:
    """Singleton class for database profiling and performance monitoring."""

    _instance = None
    _profiles: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    _enabled: bool = False  # Disabled by default - enable via DatabaseManager.enable_profiling(True)
    _slow_query_threshold: float = 1.0  # seconds

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(DBProfiler, cls).__new__(cls)
        return cls._instance

    @classmethod
    def enable(cls, enabled: bool = True):
        """Enable or disable profiling."""
        cls._enabled = enabled
        logger.info(f"Database profiling {'enabled' if enabled else 'disabled'}")

    @classmethod
    def set_slow_query_threshold(cls, threshold: float):
        """Set threshold in seconds for slow query detection."""
        cls._slow_query_threshold = threshold
        logger.info(f"Slow query threshold set to {threshold}s")

    @classmethod
    def record_operation(cls, operation_name: str, duration: float,
                        metadata: Optional[Dict[str, Any]] = None):
        """Record a database operation with timing and metadata."""
        if not cls._enabled:
            return

        record = {
            'timestamp': datetime.now().isoformat(),
            'operation': operation_name,
            'duration_seconds': duration,
            'is_slow': duration > cls._slow_query_threshold,
            'metadata': metadata or {}
        }

        cls._profiles[operation_name].append(record)

        # Log slow queries
        if record['is_slow']:
            logger.warning(
                f"Slow query detected: {operation_name} took {duration:.2f}s "
                f"(threshold: {cls._slow_query_threshold}s)"
            )

    @classmethod
    @contextmanager
    def profile(cls, operation_name: str, metadata: Optional[Dict[str, Any]] = None):
        """Context manager for profiling database operations."""
        if not cls._enabled:
            yield
            return

        start_time = time.time()
        try:
            yield
        finally:
            duration = time.time() - start_time
            cls.record_operation(operation_name, duration, metadata)

    @classmethod
    def get_statistics(cls, operation_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Get profiling statistics for an operation or all operations.

        Args:
            operation_name: Specific operation to get stats for, or None for all

        Returns:
            Dictionary with statistics (count, total_time, avg_time, min, max, slow_count)
        """
        if operation_name:
            records = cls._profiles.get(operation_name, [])
            operations = {operation_name: records}
        else:
            operations = cls._profiles

        stats = {}
        for op_name, records in operations.items():
            if not records:
                continue

            durations = [r['duration_seconds'] for r in records]
            stats[op_name] = {
                'count': len(records),
                'total_time': sum(durations),
                'avg_time': sum(durations) / len(durations),
                'min_time': min(durations),
                'max_time': max(durations),
                'slow_count': sum(1 for r in records if r['is_slow']),
                'slow_percentage': (sum(1 for r in records if r['is_slow']) / len(records)) * 100
            }

        return stats

    @classmethod
    def get_slow_queries(cls, limit: int = 10) -> List[Dict[str, Any]]:
        """Get the slowest queries recorded."""
        all_records = []
        for op_name, records in cls._profiles.items():
            all_records.extend(records)

        # Sort by duration descending
        slow_queries = sorted(
            all_records,
            key=lambda x: x['duration_seconds'],
            reverse=True
        )[:limit]

        return slow_queries

    @classmethod
    def generate_report(cls, output_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Generate a comprehensive profiling report.

        Args:
            output_path: Optional path to save JSON report

        Returns:
            Report dictionary
        """
        report = {
            'generated_at': datetime.now().isoformat(),
            'enabled': cls._enabled,
            'slow_query_threshold': cls._slow_query_threshold,
            'statistics': cls.get_statistics(),
            'slow_queries': cls.get_slow_queries(20),
            'total_operations': sum(len(records) for records in cls._profiles.values()),
            'total_time': sum(
                r['duration_seconds']
                for records in cls._profiles.values()
                for r in records
            )
        }

        # Add top operations by time
        stats = report['statistics']
        if stats:
            top_by_time = sorted(
                stats.items(),
                key=lambda x: x[1]['total_time'],
                reverse=True
            )[:10]
            report['top_operations_by_total_time'] = [
                {'operation': op, **data} for op, data in top_by_time
            ]

        if output_path:
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Profiling report saved to {output_path}")

        return report

    @classmethod
    def reset(cls):
        """Clear all profiling data."""
        cls._profiles.clear()
        logger.info("Profiling data cleared")

    @classmethod
    def explain_query(cls, conn: sqlite3.Connection, query: str,
                     params: Optional[tuple] = None) -> List[Dict[str, Any]]:
        """
        Run EXPLAIN QUERY PLAN on a SQL query.

        Args:
            conn: SQLite connection
            query: SQL query to explain
            params: Optional query parameters

        Returns:
            List of query plan steps
        """
        try:
            explain_query = f"EXPLAIN QUERY PLAN {query}"
            cursor = conn.cursor()

            if params:
                cursor.execute(explain_query, params)
            else:
                cursor.execute(explain_query)

            plan = cursor.fetchall()

            # Format results (SQLite EXPLAIN QUERY PLAN columns)
            formatted_plan = []
            for row in plan:
                formatted_plan.append({
                    'id': row[0],
                    'parent': row[1],
                    'notused': row[2],
                    'detail': row[3]
                })

            return formatted_plan

        except Exception as e:
            logger.error(f"Failed to explain query: {e}")
            return []

    @classmethod
    def get_db_size(cls, db_path: str) -> Dict[str, Any]:
        """
        Get database size metrics.

        Args:
            db_path: Path to SQLite database file

        Returns:
            Dictionary with size metrics
        """
        db_file = Path(db_path)

        if not db_file.exists():
            return {'exists': False}

        file_size_bytes = db_file.stat().st_size
        file_size_mb = file_size_bytes / (1024 * 1024)
        file_size_gb = file_size_mb / 1024

        # Get page count and page size from database
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()

            cursor.execute("PRAGMA page_count")
            page_count = cursor.fetchone()[0]

            cursor.execute("PRAGMA page_size")
            page_size = cursor.fetchone()[0]

            # Get table sizes
            cursor.execute("""
                SELECT name,
                       (SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=m.name) as table_count
                FROM sqlite_master m
                WHERE type='table' AND name NOT LIKE 'sqlite_%'
            """)
            tables = cursor.fetchall()

            table_sizes = {}
            for table_name, _ in tables:
                cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
                row_count = cursor.fetchone()[0]
                table_sizes[table_name] = row_count

            conn.close()

            return {
                'exists': True,
                'path': str(db_file),
                'size_bytes': file_size_bytes,
                'size_mb': round(file_size_mb, 2),
                'size_gb': round(file_size_gb, 2),
                'page_count': page_count,
                'page_size': page_size,
                'table_row_counts': table_sizes
            }

        except Exception as e:
            logger.error(f"Failed to get database size metrics: {e}")
            return {
                'exists': True,
                'path': str(db_file),
                'size_bytes': file_size_bytes,
                'size_mb': round(file_size_mb, 2),
                'size_gb': round(file_size_gb, 2),
                'error': str(e)
            }


def profile_db_operation(operation_name: str, explain_query: bool = False):
    """
    Decorator to profile database operations.

    Args:
        operation_name: Name of the operation for profiling
        explain_query: If True, attempt to explain SQL queries (for debugging)

    Usage:
        @profile_db_operation("insert_sequences")
        def insert_sequences_batch(conn, data):
            # ... code ...
            pass
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if not DBProfiler._enabled:
                return func(*args, **kwargs)

            # Extract metadata
            metadata = {
                'function': func.__name__,
                'args_count': len(args),
                'kwargs_count': len(kwargs)
            }

            # Try to extract row count if in kwargs
            if 'data' in kwargs and hasattr(kwargs['data'], '__len__'):
                metadata['data_size'] = len(kwargs['data'])
            elif len(args) > 1 and hasattr(args[1], '__len__'):
                metadata['data_size'] = len(args[1])

            with DBProfiler.profile(operation_name, metadata):
                result = func(*args, **kwargs)

            return result

        return wrapper
    return decorator


# Convenience function for one-off profiling
def profile_query(conn: sqlite3.Connection, query: str,
                 params: Optional[tuple] = None,
                 operation_name: str = "custom_query") -> Any:
    """
    Execute and profile a SQL query.

    Args:
        conn: SQLite connection
        query: SQL query to execute
        params: Optional query parameters
        operation_name: Name for profiling

    Returns:
        Query result
    """
    with DBProfiler.profile(operation_name):
        cursor = conn.cursor()
        if params:
            result = cursor.execute(query, params)
        else:
            result = cursor.execute(query)
        return result.fetchall()


# Example usage and testing
if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.INFO)

    # Enable profiling
    DBProfiler.enable()
    DBProfiler.set_slow_query_threshold(0.1)

    # Example: Profile a mock database operation
    @profile_db_operation("test_insert")
    def mock_insert(conn, data):
        time.sleep(0.15)  # Simulate slow operation
        return len(data)

    @profile_db_operation("test_select")
    def mock_select(conn):
        time.sleep(0.05)  # Simulate fast operation
        return []

    # Run some operations
    fake_conn = None
    for i in range(5):
        mock_insert(fake_conn, list(range(100)))
        mock_select(fake_conn)

    # Get statistics
    stats = DBProfiler.get_statistics()
    print("\n=== Profiling Statistics ===")
    print(json.dumps(stats, indent=2))

    # Get slow queries
    slow = DBProfiler.get_slow_queries(limit=3)
    print("\n=== Slowest Queries ===")
    print(json.dumps(slow, indent=2))

    # Generate report
    report = DBProfiler.generate_report("test_profiling_report.json")
    print(f"\n=== Report generated with {report['total_operations']} operations ===")
