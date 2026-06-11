#!/usr/bin/env python3
"""
Generate task manifest for the MSA pipeline.

Extracts the grouping / splitting logic from HierarchicalAligner so that
Snakemake can orchestrate each (barcode, taxon) task individually.
"""

import argparse
import json
import logging
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger("msa_task_manifest")


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )


@dataclass
class TaskConfig:
    barcodes: List[str]
    hierarchy_levels: List[str]
    tax_columns: List[str]
    domain: str
    min_domain: str
    min_seq_threshold: int
    max_seq_threshold: int
    large_task_threshold: int
    id_col: str
    seq_col: str
    barcode_col: str


def get_next_tax_level(current_level: str, tax_columns: List[str]) -> Optional[str]:
    try:
        idx = tax_columns.index(current_level)
        if idx < len(tax_columns) - 1:
            return tax_columns[idx + 1]
    except (ValueError, IndexError):
        return None
    return None


def split_group_into_tasks(
    group: pd.DataFrame,
    level: str,
    clade_name: str,
    max_sequences: int,
    id_col: str
) -> List[Dict]:
    tasks: List[Dict] = []
    if id_col in group.columns:
        group = group.sort_values(id_col, kind="mergesort")
    num_splits = (len(group) + max_sequences - 1) // max_sequences

    for idx in range(num_splits):
        start_idx = idx * max_sequences
        end_idx = min((idx + 1) * max_sequences, len(group))
        sub_group = group.iloc[start_idx:end_idx].copy()
        sub_clade_name = f"{clade_name}_part{idx + 1}"

        tasks.append({
            "level": level,
            "clade_name": sub_clade_name,
            "original_clade_name": clade_name,
            "part_number": idx + 1,
            "total_parts": num_splits,
            "sequence_count": int(len(sub_group)),
            "sequence_ids": sub_group[id_col].tolist(),
            "is_split": True
        })

    logger.info("  - %s split into %d parts (~%d seq each)", clade_name, num_splits, max_sequences)
    return tasks


def create_tasks_recursive(
    df: pd.DataFrame,
    level: str,
    barcode: str,
    cfg: TaskConfig,
    processed_indices: set
) -> List[Dict]:
    tasks: List[Dict] = []
    if level not in df.columns:
        return tasks

    for clade_name, group in df.groupby(level, sort=True):
        name = str(clade_name).replace(" ", "_")
        if cfg.id_col in group.columns:
            group = group.sort_values(cfg.id_col, kind="mergesort")

        if len(group) < cfg.min_seq_threshold:
            continue
        if isinstance(clade_name, str) and "incertae sedis" in clade_name.lower():
            continue

        if len(group) <= cfg.max_seq_threshold:
            tasks.append({
                "level": level,
                "clade_name": name,
                "sequence_count": int(len(group)),
                "sequence_ids": group[cfg.id_col].tolist(),
                "is_split": False
            })
            processed_indices.update(group.index)
        else:
            if level != cfg.min_domain:
                next_level = get_next_tax_level(level, cfg.tax_columns)
                if next_level and next_level in group.columns and not group[next_level].dropna().empty:
                    tasks.extend(create_tasks_recursive(group, next_level, barcode, cfg, processed_indices))
                else:
                    logger.warning("%s - Cannot descend from %s, splitting artificially", barcode, name)
                    tasks.extend(split_group_into_tasks(group, level, name, cfg.max_seq_threshold, cfg.id_col))
                    processed_indices.update(group.index)
            else:
                tasks.extend(split_group_into_tasks(group, level, name, cfg.max_seq_threshold, cfg.id_col))
                processed_indices.update(group.index)
    return tasks


def create_tasks_for_barcode(barcode_df: pd.DataFrame, barcode: str, cfg: TaskConfig) -> Tuple[List[Dict], List[Dict], pd.DataFrame]:
    logger.info("%s - total sequences %d", barcode, len(barcode_df))
    domain_counts = barcode_df[cfg.domain].dropna().value_counts()
    valid_domains = domain_counts[domain_counts >= cfg.min_seq_threshold].index
    valid_df = barcode_df[barcode_df[cfg.domain].isin(valid_domains)]

    processed_indices: set = set()
    all_tasks = create_tasks_recursive(valid_df, cfg.domain, barcode, cfg, processed_indices)

    normal_tasks: List[Dict] = []
    large_tasks: List[Dict] = []
    total_seq = 0

    for task in all_tasks:
        total_seq += task["sequence_count"]
        is_large = task["sequence_count"] >= cfg.large_task_threshold
        task["is_large"] = bool(is_large)
        if is_large:
            large_tasks.append(task)
        else:
            normal_tasks.append(task)

    logger.info("%s - %d sequences assigned to %d tasks (%d normal / %d large)",
                barcode, total_seq, len(all_tasks), len(normal_tasks), len(large_tasks))

    unprocessed_df = barcode_df[~barcode_df.index.isin(processed_indices)].copy()
    logger.info("%s - %d sequences unprocessed due to thresholds", barcode, len(unprocessed_df))

    return normal_tasks, large_tasks, unprocessed_df


def build_manifest(df: pd.DataFrame, cfg: TaskConfig) -> Dict:
    manifest = {
        "parameters": asdict(cfg),
        "barcodes": {},
        "total_sequences": int(len(df))
    }

    available_barcodes = cfg.barcodes or sorted(df[cfg.barcode_col].dropna().unique())

    for barcode in available_barcodes:
        barcode_df = df[df[cfg.barcode_col] == barcode]
        if barcode_df.empty:
            logger.warning("%s - no data found", barcode)
            manifest["barcodes"][barcode] = {
                "tasks": [],
                "large_tasks": [],
                "unprocessed_sequence_ids": []
            }
            continue

        normal_tasks, large_tasks, unprocessed_df = create_tasks_for_barcode(barcode_df, barcode, cfg)
        manifest["barcodes"][barcode] = {
            "tasks": normal_tasks,
            "large_tasks": large_tasks,
            "unprocessed_sequence_ids": sorted(unprocessed_df[cfg.id_col].tolist(), key=str)
        }

    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate MSA task manifest.")
    parser.add_argument("--input-csv", required=True, help="Enriched CSV with taxonomy + sequences")
    parser.add_argument("--output-json", required=True, help="Path to manifest JSON")
    parser.add_argument("--params-json", required=True, help="JSON file with task configuration")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser.parse_args()


def load_config(path: Path) -> TaskConfig:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    missing = [key for key in [
        "barcodes", "hierarchy_levels", "tax_columns", "domain", "min_domain",
        "min_seq_threshold", "max_seq_threshold", "large_task_threshold",
        "id_col", "seq_col", "barcode_col"
    ] if key not in data]
    if missing:
        raise ValueError(f"Missing keys in params-json: {missing}")
    return TaskConfig(
        barcodes=data["barcodes"],
        hierarchy_levels=data["hierarchy_levels"],
        tax_columns=data["tax_columns"],
        domain=data["domain"],
        min_domain=data["min_domain"],
        min_seq_threshold=int(data["min_seq_threshold"]),
        max_seq_threshold=int(data["max_seq_threshold"]),
        large_task_threshold=int(data["large_task_threshold"]),
        id_col=data["id_col"],
        seq_col=data["seq_col"],
        barcode_col=data["barcode_col"]
    )


def write_manifest(input_csv: str, output_json: str, cfg: TaskConfig) -> Dict:
    """Generate and persist the manifest to disk."""
    df = pd.read_csv(input_csv)
    manifest = build_manifest(df, cfg)

    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    logger.info("Manifest written to %s", output_path)
    return manifest


def run_cli():
    """Entry point when executed as a standalone CLI script."""
    args = parse_args()
    setup_logging(args.verbose)

    cfg = load_config(Path(args.params_json))
    write_manifest(args.input_csv, args.output_json, cfg)


def run_from_snakemake(snakemake):
    """Entry point used when Snakemake executes this script."""
    params = snakemake.params
    setup_logging(bool(params.get("verbose", False)))

    cfg = TaskConfig(
        barcodes=params.get("barcodes", []),
        hierarchy_levels=params["hierarchy_levels"],
        tax_columns=params["tax_columns"],
        domain=params["domain"],
        min_domain=params["min_domain"],
        min_seq_threshold=int(params["min_seq_threshold"]),
        max_seq_threshold=int(params["max_seq_threshold"]),
        large_task_threshold=int(params["large_task_threshold"]),
        id_col=params["id_col"],
        seq_col=params["seq_col"],
        barcode_col=params["barcode_col"],
    )

    input_csv = getattr(snakemake.input, "csv", None)
    if input_csv is None:
        input_csv = snakemake.input[0]

    output_json = getattr(snakemake.output, "manifest", None)
    if output_json is None:
        output_json = snakemake.output[0]

    if not input_csv:
        raise ValueError("Snakemake input `csv` is required for task manifest generation")

    if not output_json:
        raise ValueError("Snakemake output `manifest` is required for task manifest generation")

    write_manifest(input_csv, output_json, cfg)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake(snakemake)
    else:
        run_cli()
