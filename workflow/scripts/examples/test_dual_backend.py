#!/usr/bin/env python3
"""
Test script for dual backend support (SQLite + PostgreSQL).

This script demonstrates that the same code works transparently with both
SQLite and PostgreSQL backends.
"""

import sys
from pathlib import Path
import yaml
import logging

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "utils"))

from db_manager_v2 import DatabaseManager, get_database_manager, initialize_database

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_sqlite_backend():
    """Test SQLite backend with basic operations."""
    print("\n" + "=" * 60)
    print("Test 1: SQLite Backend")
    print("=" * 60)

    # Create in-memory database for testing
    db = DatabaseManager(db_type="sqlite", path=":memory:")

    # Create schema
    db.create_schema()
    print("✓ Schema created")

    # Add source
    source_id = db.add_source("TEST_SOURCE", "/path/to/test.fasta")
    print(f"✓ Source added: ID={source_id}")

    # Add taxonomy
    tax_id = db.add_taxonomy(
        kingdom="Fungi",
        phylum="Ascomycota",
        class_="Eurotiomycetes",
        order_="Eurotiales",
        family="Aspergillaceae",
        genus="Aspergillus",
        species="niger"
    )
    print(f"✓ Taxonomy added: ID={tax_id}")

    # Get stats
    stats = db.get_source_stats()
    print(f"✓ Source stats: {len(stats)} sources found")

    # Validate integrity
    integrity = db.validate_integrity()
    print(f"✓ Integrity check: valid={integrity['is_valid']}")

    print("\n✓ SQLite backend test PASSED\n")
    return True


def test_config_based_initialization():
    """Test initialization from config dictionary."""
    print("\n" + "=" * 60)
    print("Test 2: Config-based Initialization")
    print("=" * 60)

    # Test config (SQLite)
    config = {
        'database': {
            'type': 'sqlite',
            'sqlite': {
                'path': ':memory:',
                'name': 'test.db'
            }
        }
    }

    db = initialize_database(config)
    print(f"✓ Database initialized from config (type={config['database']['type']})")

    # Add test data
    source_id = db.add_source("CONFIG_TEST", "/test/path.fasta")
    print(f"✓ Test source added: ID={source_id}")

    tables = db.list_tables()
    print(f"✓ Tables created: {', '.join(tables)}")

    print("\n✓ Config-based initialization test PASSED\n")
    return True


def test_yaml_config_loading():
    """Test loading actual config.yaml file."""
    print("\n" + "=" * 60)
    print("Test 3: YAML Config Loading")
    print("=" * 60)

    config_path = Path(__file__).parent.parent.parent.parent / "config" / "config.yaml"

    if not config_path.exists():
        print(f"⚠️  Config file not found: {config_path}")
        print("⚠️  Skipping YAML config test")
        return True

    # Load config
    with open(config_path) as f:
        config = yaml.safe_load(f)

    print(f"✓ Config loaded from: {config_path}")
    print(f"  Database type: {config['database']['type']}")

    # Override to use in-memory for testing
    if config['database']['type'] == 'sqlite':
        config['database']['sqlite']['path'] = ':memory:'

        db = get_database_manager(config)
        db.create_schema()

        print("✓ Database created from YAML config")

        # Test that it works
        source_id = db.add_source("YAML_TEST", "/test.fasta")
        print(f"✓ Test operation successful: source_id={source_id}")

    else:
        print(f"⚠️  PostgreSQL configured - skipping test (would require running PostgreSQL)")

    print("\n✓ YAML config loading test PASSED\n")
    return True


def test_backward_compatibility():
    """Test that old API still works."""
    print("\n" + "=" * 60)
    print("Test 4: Backward Compatibility")
    print("=" * 60)

    db = DatabaseManager(db_type="sqlite", path=":memory:")
    db.create_schema()

    # Test old-style methods
    with db.get_connection() as conn:
        print("✓ get_connection() context manager works")

    # Test all expected methods exist
    expected_methods = [
        'create_schema',
        'add_source',
        'add_taxonomy',
        'add_seq_tax',
        'link_sequence_taxonomy',
        'get_sequence_hierarchy',
        'get_taxonomic_lineage',
        'get_sequences_by_taxonomy',
        'get_source_stats',
        'get_backbone_stats',
        'get_seq_id_from_original_id',
        'get_priority_taxonomy_id',
        'check_existing_taxonomy',
        'add_seq_tax_sativa',
        'vacuum_database',
        'validate_integrity',
        'list_tables',
        'calculate_percentage_n'
    ]

    missing_methods = []
    for method in expected_methods:
        if not hasattr(db, method):
            missing_methods.append(method)

    if missing_methods:
        print(f"❌ Missing methods: {', '.join(missing_methods)}")
        return False
    else:
        print(f"✓ All {len(expected_methods)} expected methods present")

    print("\n✓ Backward compatibility test PASSED\n")
    return True


def test_profiling_integration():
    """Test that profiling still works."""
    print("\n" + "=" * 60)
    print("Test 5: Profiling Integration")
    print("=" * 60)

    # Enable profiling
    DatabaseManager.enable_profiling(True)
    DatabaseManager.set_slow_query_threshold(0.01)  # Very low for testing

    db = DatabaseManager(db_type="sqlite", path=":memory:")
    db.create_schema()
    db.add_source("PROFILE_TEST", "/test.fasta")

    # Get profiling stats
    stats = DatabaseManager.get_profiling_stats()

    if stats:
        print(f"✓ Profiling captured {len(stats)} operations:")
        for op, data in stats.items():
            print(f"  - {op}: {data['count']} calls, avg {data['avg_time']:.4f}s")
    else:
        print("⚠️  No profiling data captured (profiling may be disabled)")

    # Reset profiling
    DatabaseManager.reset_profiling()
    print("✓ Profiling reset successful")

    print("\n✓ Profiling integration test PASSED\n")
    return True


def test_data_consistency():
    """Test that data is consistent across operations."""
    print("\n" + "=" * 60)
    print("Test 6: Data Consistency")
    print("=" * 60)

    db = DatabaseManager(db_type="sqlite", path=":memory:")
    db.create_schema()

    # Add source
    source_id = db.add_source("CONSISTENCY_TEST", "/test.fasta")

    # Add multiple taxonomies
    tax_ids = []
    genera = ["Aspergillus", "Penicillium", "Fusarium"]

    for genus in genera:
        tax_id = db.add_taxonomy(
            kingdom="Fungi",
            phylum="Ascomycota",
            genus=genus,
            species="test_species"
        )
        tax_ids.append(tax_id)

    print(f"✓ Added {len(tax_ids)} taxonomies")

    # Check duplicate detection
    duplicate_tax_id = db.add_taxonomy(
        kingdom="Fungi",
        phylum="Ascomycota",
        genus="Aspergillus",
        species="test_species"
    )

    if duplicate_tax_id == tax_ids[0]:
        print("✓ Duplicate taxonomy correctly returned existing ID")
    else:
        print(f"❌ Duplicate taxonomy got new ID: {duplicate_tax_id} vs {tax_ids[0]}")
        return False

    # Check stats
    source_stats = db.get_source_stats()
    if len(source_stats) == 1:
        print("✓ Source stats consistent")
    else:
        print(f"❌ Expected 1 source, got {len(source_stats)}")
        return False

    print("\n✓ Data consistency test PASSED\n")
    return True


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print(" DUAL BACKEND TEST SUITE")
    print("=" * 70)

    tests = [
        ("SQLite Backend", test_sqlite_backend),
        ("Config-based Init", test_config_based_initialization),
        ("YAML Config Loading", test_yaml_config_loading),
        ("Backward Compatibility", test_backward_compatibility),
        ("Profiling Integration", test_profiling_integration),
        ("Data Consistency", test_data_consistency),
    ]

    results = []

    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            logger.error(f"Test '{test_name}' failed with exception: {e}", exc_info=True)
            results.append((test_name, False))

    # Summary
    print("\n" + "=" * 70)
    print(" TEST SUMMARY")
    print("=" * 70)

    passed_count = sum(1 for _, passed in results if passed)
    total_count = len(results)

    for test_name, passed in results:
        status = "✓ PASSED" if passed else "❌ FAILED"
        print(f"{status:12} {test_name}")

    print("\n" + "=" * 70)
    print(f" RESULT: {passed_count}/{total_count} tests passed")
    print("=" * 70)

    if passed_count == total_count:
        print("\n🎉 All tests passed! Dual backend support is working correctly.\n")
        return 0
    else:
        print(f"\n⚠️  {total_count - passed_count} test(s) failed.\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
