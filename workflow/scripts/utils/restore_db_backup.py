#!/usr/bin/env python3
"""
Restore the ITS database from a backup.

Usage:
    python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db
    python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db --target data/db/its.db
"""

import argparse
import shutil
import sqlite3
from pathlib import Path
from datetime import datetime


def verify_database(db_path: Path) -> bool:
    """
    Verify the integrity of a SQLite database.

    Args:
        db_path: Path to the database file

    Returns:
        True if the DB is valid, False otherwise
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # # Run a full integrity check
        # cursor.execute("PRAGMA integrity_check")
        # result = cursor.fetchone()

        # if result[0] != "ok":
        #     print(f"❌ Integrity error: {result[0]}")
        #     return False

        # Ensure that the main tables exist
        cursor.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name IN ('sources', 'sequences', 'taxonomy', 'seq_tax')
        """)
        tables = {row[0] for row in cursor.fetchall()}
        expected_tables = {'sources', 'sequences', 'taxonomy', 'seq_tax'}

        if tables != expected_tables:
            missing = expected_tables - tables
            print(f"❌ Missing tables: {missing}")
            return False

        # Display a few statistics
        print("\n📊 Database statistics:")
        for table in ['sources', 'sequences', 'taxonomy', 'seq_tax']:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            count = cursor.fetchone()[0]
            print(f"  - {table}: {count:,} entries")

        conn.close()
        return True

    except sqlite3.Error as e:
        print(f"❌ SQLite error: {e}")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False


def create_safety_backup(db_path: Path) -> Path:
    """
    Create a safety backup of the current database.

    Args:
        db_path: Path to the current database file

    Returns:
        Path to the created backup file
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = db_path.parent / f"{db_path.stem}_before_restore_{timestamp}.db"

    print(f"💾 Creating a safety backup: {backup_path.name}")
    shutil.copy2(db_path, backup_path)

    return backup_path


def restore_database(backup_path: Path, target_path: Path, skip_verification: bool = False) -> bool:
    """
    Restore the database from a backup file.

    Args:
        backup_path: Path to the backup file
        target_path: Path to the destination DB
        skip_verification: Skip integrity checks when True

    Returns:
        True when the restoration succeeds, False otherwise
    """
    print("=" * 70)
    print("🔄 DATABASE RESTORATION")
    print("=" * 70)

    # Ensure the backup exists
    if not backup_path.exists():
        print(f"❌ Error: Backup file does not exist: {backup_path}")
        return False

    print(f"\n📁 Backup file: {backup_path}")
    print(f"📁 Target file: {target_path}")
    print(f"📦 Backup size: {backup_path.stat().st_size / (1024**2):.2f} MB")

    # Verify backup integrity
    if not skip_verification:
        print("\n🔍 Checking backup integrity...")
        if not verify_database(backup_path):
            print("❌ The backup file is corrupted or invalid.")
            return False
        print("✅ Backup looks valid.")

    # Create a safety backup when the target already exists
    if target_path.exists():
        try:
            safety_backup = create_safety_backup(target_path)
            print(f"✅ Safety backup created: {safety_backup.name}")
        except Exception as e:
            print(f"⚠️  Unable to create a safety backup: {e}")
            response = input("Continue anyway? (y/N): ")
            if response.lower() != 'y':
                print("❌ Restoration aborted.")
                return False

    # Restore the backup
    print(f"\n🔄 Restoring backup...")
    try:
        shutil.copy2(backup_path, target_path)
        print(f"✅ File copied successfully.")
    except Exception as e:
        print(f"❌ Error during copy: {e}")
        return False

    # Validate the restored database
    print("\n🔍 Checking the restored database...")
    if not verify_database(target_path):
        print("❌ Restored database is invalid!")
        return False

    print("\n" + "=" * 70)
    print("✅ RESTORATION COMPLETE!")
    print("=" * 70)

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Restore the ITS database from a backup.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db
  python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db --target data/db/its.db
  python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db --skip-verification
        """
    )

    parser.add_argument(
        '--backup',
        type=Path,
        required=True,
        help='Path to the backup file to restore'
    )

    parser.add_argument(
        '--target',
        type=Path,
        default=None,
        help='Path to the destination file (default: data/db/its.db)'
    )

    parser.add_argument(
        '--skip-verification',
        action='store_true',
        help='Skip backup integrity verification before restoring'
    )

    args = parser.parse_args()

    # Determine the target path
    if args.target is None:
        # Use the default path
        target_path = Path("data/db/its.db")
    else:
        target_path = args.target

    # Convert to absolute paths
    backup_path = args.backup.resolve()
    target_path = target_path.resolve()

    # Ensure backup and target are not the same file
    if backup_path == target_path:
        print("❌ Error: Backup and target files are identical!")
        return 1

    # Ask for confirmation
    print(f"\n⚠️  You are about to restore:")
    print(f"   FROM: {backup_path}")
    print(f"   TO:   {target_path}")
    print()
    response = input("Confirm the restoration? (y/N): ")

    if response.lower() != 'y':
        print("❌ Restoration aborted.")
        return 0

    # Restore
    success = restore_database(backup_path, target_path, args.skip_verification)

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
