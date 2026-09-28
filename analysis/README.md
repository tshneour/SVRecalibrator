# Delly comparison analysis

Standalone scripts for consolidating SVRecalibrator and Delly outputs and
performing informal concordance analysis.

## Consolidating Delly calls

Build a spreadsheet-readable CSV with one row per existing SVRecalibrator row:

```bash
python analysis/consolidate_delly.py /path/to/batch_outputs /path/to/delly_vcfs \
  -o /path/to/comparison --window 100
```

Input VCFs must be single-sample Delly files named `<sample>.vcf` or
`<sample>.vcf.gz`; sample names match SVRecalibrator directory names exactly.
The original `sample` column is preserved, while `sample_id` comes from the
parent directory (older outputs sometimes use amplicon labels in `sample`).
No additional dependencies are required. These standalone analysis scripts are
not installed as part of the SVRecalibrator package. Run the commands below
from the repository root.

Outputs:

- `consolidated.csv`: all original columns plus the nearest compatible Delly
  call, distances, quality flags, and ambiguity/reuse counts. Unmatched SV rows
  remain in the table.
- `candidates.csv`: all nearby candidate pairs, including orientation conflicts.
- `delly_calls.csv`: all Delly calls for the SVRecalibrator sample directories,
  including unselected calls and complete INFO/FORMAT fields.
- `sample_inventory.csv`: nonempty, empty, or missing SV tables and missing VCFs.
- `manifest.json`: matching settings, input checksums, script checksum, and VCF
  filenames without corresponding SVRecalibrator sample directories.

Matching uses the **original** `break_pos1/2`, with both ends within the window.
The script compares the reported coordinates directly, strips an optional `chr`
for matching, and checks both endpoint orders. Delly CT maps `3` to `+` and `5`
to `-`; swapping endpoints swaps their signs. BND/TRA partners use `CHR2/POS2`,
other calls use `END`. Coordinates, sequences, and chromosome names in the
source fields remain unchanged. Candidates with compatible orientation are
ranked by maximum endpoint distance, then summed distance, then record key.
Ties and Delly calls selected by multiple SV rows are flagged; this is not a
one-to-one assignment. All filters and genotypes are retained, including
LowQual calls; variant type is preserved but is not an additional matching gate.

Run the informal analysis separately after reviewing the consolidation:

```bash
python analysis/analyze_delly_concordance.py /path/to/comparison/consolidated.csv \
  -o /path/to/comparison/analysis
```

This writes per-row details and a JSON summary with coordinate-tolerance counts
and homology-length comparisons for all matches and the PASS/PRECISE subset.
It compares split/scaffold sequences directly and as reverse complements.
Delly `CONSENSUS` is a full junction sequence, **not** a homology sequence, so it
is retained for downstream work rather than treated as directly equivalent to
SVRecalibrator `hom` or `sc_hom`. Delly `HOMLEN` permits edits; differences from
SVRecalibrator lengths need further interpretation. No coordinate-offset or
microhomology-shift normalization is applied. Missing lengths remain unknown.
These are SVRecalibrator-ascertained comparisons, not sensitivity or precision
estimates; genome-wide Delly-only calls are not automatically false positives.

## Comparing junctions across callers

`compare_callers.py` matches every SVRecalibrator row to Delly2, SvABA and
GRIDSS2 calls and compares junction homology/insertion lengths. See
[junction_comparison_methods.md](junction_comparison_methods.md) for the
matching and junction-comparison methods.

```bash
python compare_callers.py SVR_BATCH_DIR -f GENOME.fa -o OUTDIR \
    --delly DELLY_DIR --svaba SVABA_DIR --gridss GRIDSS_DIR
```

The Turner 2017 inputs are packaged per tool in `turner2017/` (see its README); run them
with `python compare_callers.py --bundle -f GENOME.fa -o OUTDIR`.

Expected inputs per sample `S`: `SVR_BATCH_DIR/S/final_augmented.tsv`,
`DELLY_DIR/S.vcf[.gz]`, `SVABA_DIR/S/S.svaba.unfiltered.sv.vcf`, and
`GRIDSS_DIR/S.vcf.gz`. The reference FASTA (indexed) is used to rebuild Delly2
junctions from `CONSENSUS`. By default only samples with every caller's VCF are
included (`--all-samples` to override); `--pass-only` keeps only PASS calls and
`--window` sets the matching distance (default 100 bp). Outputs:
`callers_consolidated.csv`, `sample_inventory.csv`, `summary.json`, and
junction-length histograms and per-tool category bar charts for SVs identified
by all callers and for all rows.

## Tests

Run the small fixture tests from this directory:

```bash
cd analysis
python -m unittest discover -s tests -v
```
