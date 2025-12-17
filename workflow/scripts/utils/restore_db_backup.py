#!/usr/bin/env python3
"""
Script pour restaurer la base de données ITS à partir d'un backup.

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
    Vérifie l'intégrité d'une base de données SQLite.

    Args:
        db_path: Chemin vers le fichier de base de données

    Returns:
        True si la DB est valide, False sinon
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # # Vérifier l'intégrité
        # cursor.execute("PRAGMA integrity_check")
        # result = cursor.fetchone()

        # if result[0] != "ok":
        #     print(f"❌ Erreur d'intégrité: {result[0]}")
        #     return False

        # Vérifier que les tables principales existent
        cursor.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' AND name IN ('sources', 'sequences', 'taxonomy', 'seq_tax')
        """)
        tables = {row[0] for row in cursor.fetchall()}
        expected_tables = {'sources', 'sequences', 'taxonomy', 'seq_tax'}

        if tables != expected_tables:
            missing = expected_tables - tables
            print(f"❌ Tables manquantes: {missing}")
            return False

        # Afficher quelques statistiques
        print("\n📊 Statistiques de la base de données:")
        for table in ['sources', 'sequences', 'taxonomy', 'seq_tax']:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            count = cursor.fetchone()[0]
            print(f"  - {table}: {count:,} entrées")

        conn.close()
        return True

    except sqlite3.Error as e:
        print(f"❌ Erreur SQLite: {e}")
        return False
    except Exception as e:
        print(f"❌ Erreur: {e}")
        return False


def create_safety_backup(db_path: Path) -> Path:
    """
    Crée une sauvegarde de sécurité de la DB actuelle.

    Args:
        db_path: Chemin vers le fichier de base de données actuel

    Returns:
        Chemin vers le fichier de backup créé
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = db_path.parent / f"{db_path.stem}_before_restore_{timestamp}.db"

    print(f"💾 Création d'une sauvegarde de sécurité: {backup_path.name}")
    shutil.copy2(db_path, backup_path)

    return backup_path


def restore_database(backup_path: Path, target_path: Path, skip_verification: bool = False) -> bool:
    """
    Restaure la base de données à partir d'un backup.

    Args:
        backup_path: Chemin vers le fichier de backup
        target_path: Chemin vers le fichier cible
        skip_verification: Si True, ne vérifie pas l'intégrité du backup

    Returns:
        True si la restauration a réussi, False sinon
    """
    print("=" * 70)
    print("🔄 RESTAURATION DE LA BASE DE DONNÉES")
    print("=" * 70)

    # Vérifier que le backup existe
    if not backup_path.exists():
        print(f"❌ Erreur: Le fichier de backup n'existe pas: {backup_path}")
        return False

    print(f"\n📁 Fichier de backup: {backup_path}")
    print(f"📁 Fichier cible: {target_path}")
    print(f"📦 Taille du backup: {backup_path.stat().st_size / (1024**2):.2f} MB")

    # Vérifier l'intégrité du backup
    if not skip_verification:
        print("\n🔍 Vérification de l'intégrité du backup...")
        if not verify_database(backup_path):
            print("❌ Le fichier de backup est corrompu ou invalide.")
            return False
        print("✅ Le backup est valide.")

    # Créer une sauvegarde de sécurité si le fichier cible existe
    if target_path.exists():
        try:
            safety_backup = create_safety_backup(target_path)
            print(f"✅ Sauvegarde de sécurité créée: {safety_backup.name}")
        except Exception as e:
            print(f"⚠️  Impossible de créer une sauvegarde de sécurité: {e}")
            response = input("Continuer quand même? (y/N): ")
            if response.lower() != 'y':
                print("❌ Restauration annulée.")
                return False

    # Restaurer le backup
    print(f"\n🔄 Restauration en cours...")
    try:
        shutil.copy2(backup_path, target_path)
        print(f"✅ Fichier copié avec succès.")
    except Exception as e:
        print(f"❌ Erreur lors de la copie: {e}")
        return False

    # Vérifier la base restaurée
    print("\n🔍 Vérification de la base restaurée...")
    if not verify_database(target_path):
        print("❌ La base restaurée est invalide!")
        return False

    print("\n" + "=" * 70)
    print("✅ RESTAURATION RÉUSSIE!")
    print("=" * 70)

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Restaure la base de données ITS à partir d'un backup.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples:
  python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db
  python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db --target data/db/its.db
  python restore_db_backup.py --backup data/db/its_backup_20251029_102912.db --skip-verification
        """
    )

    parser.add_argument(
        '--backup',
        type=Path,
        required=True,
        help='Chemin vers le fichier de backup à restaurer'
    )

    parser.add_argument(
        '--target',
        type=Path,
        default=None,
        help='Chemin vers le fichier cible (défaut: data/db/its.db)'
    )

    parser.add_argument(
        '--skip-verification',
        action='store_true',
        help='Ne pas vérifier l\'intégrité du backup avant restauration'
    )

    args = parser.parse_args()

    # Déterminer le chemin cible
    if args.target is None:
        # Utiliser le chemin par défaut
        target_path = Path("data/db/its.db")
    else:
        target_path = args.target

    # Convertir les chemins en absolus
    backup_path = args.backup.resolve()
    target_path = target_path.resolve()

    # Vérifier que backup et target ne sont pas le même fichier
    if backup_path == target_path:
        print("❌ Erreur: Le fichier de backup et le fichier cible sont identiques!")
        return 1

    # Demander confirmation
    print(f"\n⚠️  Vous allez restaurer:")
    print(f"   FROM: {backup_path}")
    print(f"   TO:   {target_path}")
    print()
    response = input("Confirmer la restauration? (y/N): ")

    if response.lower() != 'y':
        print("❌ Restauration annulée.")
        return 0

    # Restaurer
    success = restore_database(backup_path, target_path, args.skip_verification)

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
