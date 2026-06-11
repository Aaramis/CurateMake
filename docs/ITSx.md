# ITSx Stage

The ITSx stage extracts ITS, ITS1, and ITS2 regions from curated sequences and prepares downstream manifests for the MSA phase.

## Data flow

1. `deduplicate_sequences_pre_itsx` queries the DB (via `workflow/scripts/deduplicate_sequences.py`) to produce `data/processed/itsx_deduplicate_pre.json`.
2. `prepare_itsx` writes barcode FASTA files into `data/ITSx/<barcode>.fasta`.
3. `run_itsx_barcode` (one job per barcode) launches ITSx with the CPU count defined in `config["itsx_cpu"]`.
4. `process_barcode_itsx` parses the ITSx outputs, updates `barcode_status` flags in the DB, and generates barcode-level JSON summaries.
5. `aggregate_itsx_results` collapses per-barcode stats into `data/processed/itsx_process.json`.
6. `deduplicate_sequences_post_itsx` performs a second DB sweep to drop duplicates introduced during extraction.

## Configuration knobs

Defined in `config/config.yaml`:

```yaml
itsx_cpu: 30                     # threads per ITSx invocation
itsx_deterministic: false        # set true for reproducible (single-threaded) runs
itsx_output_dir: "data"          # where raw FASTA inputs live (defaults to "data")
msa:
  barcodes: ["ITS","ITS1","ITS2"]  # reused by ITSx rules
sequence_validation:
  min_length: 50
  max_n_fraction: 0.5
```

Key points:

- `itsx_cpu` is applied to both Snakemake’s `threads` directive and the `ITSx --cpu` flag.
- `itsx_deterministic: true` forces a single thread (`--cpu 1 --multi_thread F`)
  for reproducible output, at the cost of speed.
- The rule auto-discovers barcodes by inspecting `data/ITSx/`, so additional barcodes can be injected without editing the rule file.
- Sequence validation thresholds come from the shared `sequence_validation` block; any sequence failing length or `N`-content checks will not be reinserted.

## Outputs

- `data/processed/itsx_run_<barcode>.done` – sentinel per barcode.
- `logs/ITSX/...` – three log families: deduplication, prepare, run/process.
- `data/processed/itsx_<barcode>_result.json` – stats for each barcode (counts, failures, DB updates).
- `data/processed/itsx_process.json` – merged metrics consumed by `generate_sequence_flow_report.py`.

## Tips

- Rerun ITSx for a specific barcode:
  ```bash
  snakemake --use-conda --cores 16 data/processed/itsx_run_ITS2.done
  ```
- The `barcode_status` column in the `sequences` table tracks each sequence
  through extraction (`Pending` → `ITSx`/`no_detection` → `Processed`), so you
  can audit ITSx outcomes directly in the database.
