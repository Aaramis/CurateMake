# Adding Your Own Data (PERSONAL Database)

This guide explains how to add your own ITS sequences to the pipeline alongside the standard databases (UNITE, BOLD, PLANiTS, CALeDNA).

## Overview

The PERSONAL database feature allows you to integrate custom sequences with their taxonomic annotations into the pipeline. Your data will go through the same curation process as the other databases: ITSx extraction, MSA alignment, SATIVA validation, etc.

## FASTA Header Format

Your FASTA file must follow a specific header format where taxonomic information is encoded as key-value pairs separated by semicolons.

### Basic Format

```
>SEQUENCE_ID;barcode=ITS;kingdom=Fungi;phylum=Basidiomycota;class=Agaricomycetes;order=Russulales;family=Russulaceae;genus=Russula;species=Russula_abietina
ATCGATCG...
```

### Required Fields

| Field | Description | Values |
|-------|-------------|--------|
| `barcode` | The ITS region | `ITS`, `ITS1`, or `ITS2` |

### Optional Taxonomic Fields

All taxonomic fields are optional. Missing fields will be stored as `NULL` in the database.

| Field | Description | Example |
|-------|-------------|---------|
| `kingdom` | Kingdom | `Fungi`, `Plantae` |
| `phylum` | Phylum | `Basidiomycota`, `Ascomycota` |
| `class` | Class | `Agaricomycetes` |
| `order` | Order | `Russulales` |
| `family` | Family | `Russulaceae` |
| `genus` | Genus | `Russula` |
| `species` | Species (can include genus) | `Russula_abietina` or `Russula abietina` |

### Handling Special Cases

The parser is designed to handle various edge cases:

#### 1. Spaces in Headers

Spaces after the sequence ID (before the first semicolon) are tolerated:

```
>4b6ed8c0-5a11-4103-405c-b686dc74c49b ;barcode=ITS;kingdom=Fungi;...
```

#### 2. Additional Metadata

Extra information after the taxonomy fields is automatically ignored:

```
>SEQ001;barcode=ITS;kingdom=Fungi;species=Russula_abietina,+strand,0-702 length=718 error-free_length=728 read_identity=94.681%
```

The parser will extract only the taxonomic fields and ignore `,+strand,0-702 length=718...`.

#### 3. Missing Taxonomic Ranks

You can omit any taxonomic level:

```
>SEQ001;barcode=ITS;kingdom=Fungi;genus=Russula;species=Russula_sp
```

#### 4. Unknown/NA Values

The following values are treated as missing (NULL):
- `NA`, `N/A`
- `none`, `null`
- `unknown`, `unidentified`
- Empty values (`field=`)

```
>SEQ001;barcode=ITS;kingdom=Fungi;phylum=NA;class=unknown;genus=Russula
```

#### 5. Underscores in Names

Underscores are automatically converted to spaces in taxonomic names:

```
species=Russula_abietina  ->  "Russula abietina"
```

## Configuration

### Single File

For a single FASTA file, use the `path` option:

```yaml
databases:
  personal:
    path: "/path/to/your/custom_sequences.fasta"
```

### Multiple Files

For multiple FASTA files, use the `files` option:

```yaml
databases:
  personal:
    files:
      - "/path/to/custom_sequences_1.fasta"
      - "/path/to/custom_sequences_2.fasta"
      - "/path/to/custom_sequences_3.fasta"
```

### Disabling PERSONAL Database

To disable the PERSONAL database, either:
- Remove the `personal` section from the config
- Comment it out
- Set `path` to an empty string or `null`

```yaml
databases:
  # personal:
  #   path: "/path/to/your/sequences.fasta"
```

## File Location

Place your FASTA file(s) in the `data/raw/PERSONAL/` directory (recommended) or any accessible path.

```
data/
  raw/
    PERSONAL/
      my_sequences.fasta
    UNITE/
    BOLD/
    ...
```

## Running the Pipeline

Once configured, run the pipeline normally. The PERSONAL database will be processed automatically:

```bash
# Dry run to verify configuration
snakemake --use-conda -np

# Run the pipeline
snakemake --use-conda --cores 32
```

### Verify PERSONAL Processing

Check that the PERSONAL database is included in the formatting step:

```bash
# Look for the personal formatting report
ls data/processed/personal_*_formatted_report.json
```

The report contains processing statistics:
- Total sequences
- Filtered sequences (by quality or barcode)
- Barcode distribution (ITS, ITS1, ITS2)
- Processing rate

## Example FASTA File

```fasta
>SEQ001;barcode=ITS;kingdom=Fungi;phylum=Basidiomycota;class=Agaricomycetes;order=Russulales;family=Russulaceae;genus=Russula;species=Russula_abietina
TAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTATCGTACCACCGAGGTGCAAGGG
CTGTCGCTGACCTTTGAAGGTCGTGCACGCCCGAGCGCTCTCGCCACAATCCATCTCACCCCTTTGTGCATCACCGCGTG
GGTCCCCCTTTGCGGGGAGGGCTCGCGTTTTCACACAAAACTTGATACAGTCTAGAATGTTTATTTTTGCGGTAACACGC

>SEQ002;barcode=ITS1;kingdom=Fungi;phylum=Ascomycota;genus=Penicillium;species=Penicillium_chrysogenum
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG

>SEQ003;barcode=ITS2;kingdom=Fungi;genus=Unknown_genus
GCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
GCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
```

## Sequence Validation

Sequences are validated against the same criteria as other databases:

| Parameter | Config Key | Default |
|-----------|------------|---------|
| Minimum length | `sequence_validation.min_length` | 50 bp |
| Maximum N fraction | `sequence_validation.max_n_fraction` | 0.5 (50%) |
| Valid characters | `sequence_validation.valid_chars` | `ACGTURYSWKMBDHVN-` |

Sequences failing validation are filtered out and counted in the report.

## Troubleshooting

### "No sequences processed"

1. Check that the file path in config is correct
2. Verify the FASTA format (headers start with `>`)
3. Ensure `barcode` field is present and valid (`ITS`, `ITS1`, or `ITS2`)
4. Check that sequences meet minimum length requirements

### "Parsing errors"

Check the log file at `logs/format_databases/personal_<filename>_formatting.log` for details on which headers failed to parse.

### "Barcode filtered"

Sequences with barcodes other than ITS, ITS1, or ITS2 are filtered out. Check your `barcode=` field values.

## Technical Notes

- The PERSONAL parser is implemented in `workflow/scripts/format_personal.py`
- Data is stored in the same normalized schema as other databases
- Source tracking allows querying sequences by origin (`PERSONAL` vs `UNITE`, etc.)
- Duplicate taxonomies are deduplicated across all databases
