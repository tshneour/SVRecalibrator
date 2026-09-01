# SVRecalibrator

Tools to collect reads around structural-variant (SV) breakpoints and refine predictions using split-read evidence and optional de-novo scaffolds.

## Prerequisites

* Linux / macOS (bash)
* **git**
* **conda** / **mamba**

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/tshneour/SVRecalibrator.git
cd SVRecalibrator
```

### 2. Create conda environment and install dependencies

Use `environment.yml` to create the `sv-analysis` environment:

```bash
conda env create -f environment.yml
conda activate sv-analysis
```

### 3. Install SVRecalibrator package CLI

Install in editable mode so the unified `sv-recalibrator` command is available on your PATH:

```bash
pip install -e .
```

---

## Quickstart

### Single Sample Execution

To refine breakpoints for a single sample BAM and SV summary directory:

```bash
sv-recalibrator run \
  --bam /path/to/sample.bam \
  --sv-dir /path/to/AA_summaries \
  -o /path/to/output_dir \
  -f /path/to/GRCh38.fa \
  --strict
```

> **Note on Disk Space**: By default, `sv-recalibrator run` automatically cleans up temporary FASTQ files and SPAdes assembly build directories upon completion to save disk space on large BAM runs. To keep these intermediate files, pass `--keep-intermediates`.

---

## Unified CLI Reference (`sv-recalibrator`)

### 1. `sv-recalibrator run` (Single Sample)

Runs the complete collection, alignment, and breakpoint refinement workflow for one sample.

```bash
sv-recalibrator run \
  --bam path/to/sample.bam \
  --sv-dir path/to/AA_summaries \
  --outdir path/to/output_dir \
  --fasta path/to/ref.fa \
  [--sample SAMPLE_NAME] \
  [--mode {split,scaffold,both}] \
  [--radius 350] \
  [--threads 16] \
  [--spades-timeout 2.0] \
  [--strict] \
  [--keep-intermediates]
```

**Options:**

| Option | Description |
|---|---|
| `--bam BAM` | **(Required)** Path to coordinate-sorted, indexed BAM file. |
| `--sv-dir DIR` | **(Required)** Directory containing AmpliconArchitect SV summary TSVs (or single file). |
| `-o, --outdir DIR` | **(Required)** Output directory for output TSVs and logs. |
| `-f, --fasta FASTA` | **(Required for scaffold/both)** Indexed reference genome FASTA. |
| `-s, --sample NAME` | Sample prefix for filtering TSVs. If omitted, derived from BAM filename. |
| `-m, --mode MODE` | Refinement mode: `split`, `scaffold`, or `both` (default: `both`). |
| `-r, --radius INT` | Collection window radius around breakpoints in bp (default: 350). |
| `-t, --threads INT` | SPAdes assembly threads (default: 16). |
| `--spades-timeout HR` | Per-breakpoint SPAdes assembly timeout in hours (default: 2.0). |
| `--strict` | Filter for read pairs fully aligning within region of interest. |
| `--keep-intermediates` | Preserve intermediate FASTQs and SPAdes folders (default: auto-clean). |

### 2. `sv-recalibrator batch` (Multi-Sample Batch)

Runs batch execution across multiple samples listed in a sample file.

```bash
sv-recalibrator batch \
  --samples phase2_samples.txt \
  --sv-dir /path/to/SV_summaries \
  --outdir /path/to/batch_outputs \
  --fasta /path/to/ref.fa \
  --strict
```

### 3. `sv-recalibrator summarize` (Summary & Plots)

Generates summary CSV tables and Euler/histogram figures comparing refined breakpoints against baseline AA calls.

```bash
sv-recalibrator summarize /path/to/batch_outputs -o /path/to/summary_prefix
```

---

## Input format

### SV summary (`sum`) files

`collect.py` expects one or more **tab-separated (`.tsv`) files** describing structural-variant breakpoints with a fixed set of required columns.

Each TSV in the `sum/` directory represents a collection of SV breakpoints which are required to be **coordinate-sorted**. When `--sample` is provided, `collect.py` filters the SV summary directory to TSV files matching `<SAMPLE>_*.tsv` or `<SAMPLE>.tsv`.

#### Required columns

Each input TSV **must** contain the following columns (case-sensitive):

| Column name         | Type   | Description                                                                                                                  |
| ------------------- | ------ | -----------------------------------------------------------------------------------------------------------------------------|
| `chrom1`            | string | Chromosome of the first breakpoint end (e.g. `chr8`)                                                                         |
| `pos1`              | int    | 0-based genomic coordinate of the first breakpoint                                                                           |
| `chrom2`            | string | Chromosome of the second breakpoint end                                                                                      |
| `pos2`              | int    | 0-based genomic coordinate of the second breakpoint                                                                          |
| `sv_type`           | string | Structural variant type (i.e. `deletion`, `duplication`, `interchromosomal`, `inversion`, `foldback`)                        |
| `orientation`       | string | Breakpoint orientation as a 2-character string (`++`, `--`, `+-`, `-+`)                                                      |

Coordinates in the output tables are reported as **1-based**, but the input `pos1` / `pos2` values are treated as **0-based** internally. Note that `sv_type` is a case-sensitive field.

---

#### Optional columns

If provided, these columns (case-sensitive) will be included in the final output for comparison purposes.

| Column name         | Type   | Description                                                                               |
| ------------------- | ------ | ----------------------------------------------------------------------------------------- |
| `features`          | string | Arbitrary annotation or metadata for the breakpoint                                       |
| `read_support`      | int    | Number of reads supporting the breakpoint (used for reporting only)                       |
| `homology_length`   | int    | Length of homology reported for the breakpoint (may be 0)                                 |
| `homology_sequence` | string | Homology or inserted sequence (may be empty)                                              |

---

## Underlying Engine Tools

### `collect.py`

Collect reads around SV breakpoints for refinement and write per-breakpoint paired FASTQs plus a combined TSV of read evidence.

```
usage: collect.py [-h] [-s SAMPLE] [-v] [--strict] [-f FILE] refine sum bam
```

### `refine.py`

Refine SV breakpoints using split-read evidence and/or local de-novo scaffold reconstruction.

```
usage: refine.py FILE [--mode {split,scaffold,both}]
                     [--out-table OUT] [--split-log PATH]
                     [--scaffold-log PATH] [--outdir DIR]
                     [--clean] [--keep-intermediates]
                     [-l | --list] [-b IDX [IDX ...]] [-v ...]
                     [--fasta FASTA]
```

---

## Working assumptions & tips

* **BAM**: Coordinate-sorted, indexed (`.bai` present).
* **SV summaries**: Input TSVs must conform to the column specification above.
* **Mapping quality**: Reads with MAPQ > 15 (or mapped status) are retained.
* **FASTA**: Required for scaffold mode; must be indexed first:

  ```bash
  samtools faidx /path/to/genome.fa
  ```
