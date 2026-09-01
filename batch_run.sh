#!/usr/bin/env bash
# batch_run.sh — Run SVRecalibrator (collect + refine) across multiple samples.
#
# Each sample requires a BAM file (coordinate-sorted, indexed) and one or more
# SV summary TSV files named <SAMPLE>_amplicon*_SV_summary.tsv (or just
# <SAMPLE>_SV_summary.tsv) in a shared SV summary directory.
#
# The sample-to-BAM mapping is provided as a tab-separated file with columns:
#   SAMPLE_NAME  /path/to/sample.bam  [any extra columns are ignored]
#
# Outputs per sample:
#   <outdir>/<SAMPLE>/alignments.tsv           collect.py output
#   <outdir>/<SAMPLE>/fastq/                   per-breakpoint FASTQ pairs
#   <outdir>/<SAMPLE>/final_augmented.tsv      refine.py output
#   <outdir>/<SAMPLE>/split_read_alignments    split-read log
#   <outdir>/<SAMPLE>/scaffold_alignments      scaffold log
#   <outdir>/<SAMPLE>/collect.{stdout,stderr}.log
#   <outdir>/<SAMPLE>/refine.{stdout,stderr}.log
#   <outdir>/<SAMPLE>/done.flag                written on success; skip if present
#
# Usage:
#   bash batch_run.sh [OPTIONS]
#
# Required options:
#   -s, --samples FILE    Tab-separated sample list: SAMPLE  BAM_PATH  [...]
#   -v, --sv-dir DIR      Directory containing all SV summary TSV files
#   -o, --outdir DIR      Output root directory (created if absent)
#   -f, --fasta FILE      Indexed reference FASTA (required for scaffold/both mode)
#
# Optional options:
#   -m, --mode MODE       Refinement mode: split | scaffold | both  (default: both)
#   -r, --radius INT      Collect radius in bp around each breakpoint (default: 350)
#   -t, --threads INT     SPAdes threads per sample (default: 16)
#       --spades-timeout HOURS  Per-breakpoint SPAdes time limit in hours; 0 = no limit (default: 2)
#       --strict          Pass --strict to collect.py (keep only reads fully within region)
#       --sample NAME     Run a single named sample instead of all
#       --clean           Clean intermediate FASTQ and SPAdes directories after run
#       --keep-intermediates Preserve intermediate directories (overrides --clean)
#   -h, --help            Show this help and exit
#
# Example:
#   bash batch_run.sh \
#     -s phase2_samples.txt \
#     -v /data/SV_summaries \
#     -o /data/SVRecalibrator_outputs \
#     -f /data/ref/hg38.fa \
#     --strict
#
# To run a single sample:
#   bash batch_run.sh ... --sample ACHN

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── defaults ───────────────────────────────────────────────────────────────────
SAMPLES_FILE=""
SV_SUMMARY_DIR=""
OUTPUT_BASE=""
FASTA=""
MODE="both"
RADIUS=350
THREADS=16
SPADES_TIMEOUT=2
STRICT=""
SINGLE_SAMPLE=""
CLEAN=""
KEEP_INTERMEDIATES=""

# ── argument parsing ───────────────────────────────────────────────────────────
usage() {
    sed -n '/^# Usage:/,/^[^#]/{ /^#/{ s/^# \?//; p } }' "$0"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samples)   SAMPLES_FILE="$2"; shift 2 ;;
        -v|--sv-dir)    SV_SUMMARY_DIR="$2"; shift 2 ;;
        -o|--outdir)    OUTPUT_BASE="$2"; shift 2 ;;
        -f|--fasta)     FASTA="$2"; shift 2 ;;
        -m|--mode)      MODE="$2"; shift 2 ;;
        -r|--radius)    RADIUS="$2"; shift 2 ;;
        -t|--threads)   THREADS="$2"; shift 2 ;;
        --spades-timeout) SPADES_TIMEOUT="$2"; shift 2 ;;
        --strict)       STRICT="--strict"; shift ;;
        --sample)       SINGLE_SAMPLE="$2"; shift 2 ;;
        --clean)        CLEAN="true"; shift ;;
        --keep-intermediates) KEEP_INTERMEDIATES="true"; shift ;;
        -h|--help)      usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

# ── validate required args ─────────────────────────────────────────────────────
missing=()
[[ -z "$SAMPLES_FILE" ]]   && missing+=("-s/--samples")
[[ -z "$SV_SUMMARY_DIR" ]] && missing+=("-v/--sv-dir")
[[ -z "$OUTPUT_BASE" ]]    && missing+=("-o/--outdir")
if [[ "$MODE" != "split" && -z "$FASTA" ]]; then
    missing+=("-f/--fasta (required for mode=$MODE)")
fi
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "Missing required options: ${missing[*]}" >&2
    usage 1
fi

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "Samples file not found: $SAMPLES_FILE" >&2; exit 1
fi
if [[ ! -d "$SV_SUMMARY_DIR" ]]; then
    echo "SV summary directory not found: $SV_SUMMARY_DIR" >&2; exit 1
fi
if [[ -n "$FASTA" && ! -f "$FASTA" ]]; then
    echo "FASTA not found: $FASTA" >&2; exit 1
fi

# ── helpers ────────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%H:%M:%S')] $*"; }
err() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; }

run_sample() {
    local SAMPLE="$1"

    # Resolve BAM from sample list
    local BAM
    BAM=$(awk -v s="$SAMPLE" '$1 == s { print $2; exit }' "$SAMPLES_FILE")
    if [[ -z "$BAM" ]]; then
        err "$SAMPLE: not found in $SAMPLES_FILE — skipping"
        return 1
    fi
    if [[ ! -f "$BAM" ]]; then
        err "$SAMPLE: BAM not found: $BAM — skipping"
        return 1
    fi
    if [[ ! -f "${BAM}.bai" ]]; then
        err "$SAMPLE: BAM index not found: ${BAM}.bai — skipping"
        return 1
    fi

    # Find SV summary files for this sample (supports both naming conventions)
    local SV_FILES
    mapfile -t SV_FILES < <(find "$SV_SUMMARY_DIR" -maxdepth 1 \
        \( -name "${SAMPLE}_amplicon*_SV_summary.tsv" -o -name "${SAMPLE}_SV_summary.tsv" \) \
        | sort)
    if [[ ${#SV_FILES[@]} -eq 0 ]]; then
        err "$SAMPLE: no SV summary files in $SV_SUMMARY_DIR — skipping"
        return 1
    fi

    local OUT_DIR="${OUTPUT_BASE}/${SAMPLE}"
    local DONE_FLAG="${OUT_DIR}/done.flag"
    if [[ -f "$DONE_FLAG" ]]; then
        log "$SAMPLE: already complete (remove ${DONE_FLAG} to rerun)"
        return 0
    elif [[ -d "$OUT_DIR" ]]; then
        log "$SAMPLE: incomplete output directory found — clearing and rerunning"
        rm -rf "$OUT_DIR"
    fi

    log "$SAMPLE: starting (${#SV_FILES[@]} SV file(s), BAM=$(basename "$BAM"))"
    mkdir -p "$OUT_DIR"

    # Create a per-sample SV summary subdir with symlinks so collect.py sees only
    # this sample's files
    local SV_SUBDIR="${OUT_DIR}/sv_summaries"
    mkdir -p "$SV_SUBDIR"
    for f in "${SV_FILES[@]}"; do
        ln -sf "$f" "$SV_SUBDIR/$(basename "$f")"
    done

    local ALIGNMENTS_TSV="${OUT_DIR}/alignments.tsv"
    local OUT_TABLE_STEM="${OUT_DIR}/final_augmented"
    local SPADES_DIR="${OUT_DIR}/spades"

    # ── Step 1: collect ────────────────────────────────────────────────────────
    log "$SAMPLE: collect.py (radius=${RADIUS}bp${STRICT:+, strict})"

    # collect.py writes fastq/ relative to cwd; run from output dir
    pushd "$OUT_DIR" > /dev/null
    if ! python "${SCRIPT_DIR}/collect.py" \
            "$RADIUS" \
            "$SV_SUBDIR" \
            "$BAM" \
            $STRICT \
            -f "$ALIGNMENTS_TSV" \
            > "${OUT_DIR}/collect.stdout.log" \
            2> "${OUT_DIR}/collect.stderr.log"; then
        popd > /dev/null
        err "$SAMPLE: collect.py failed — see ${OUT_DIR}/collect.stderr.log"
        return 1
    fi
    popd > /dev/null

    if [[ ! -f "$ALIGNMENTS_TSV" ]]; then
        err "$SAMPLE: collect.py produced no alignments.tsv"
        return 1
    fi

    # ── Step 2: refine ─────────────────────────────────────────────────────────
    log "$SAMPLE: refine.py (mode=${MODE})"

    local REFINE_ARGS=(
        "$ALIGNMENTS_TSV"
        --mode "$MODE"
        --out-table "$OUT_TABLE_STEM"
        --split-log "${OUT_DIR}/split_read_alignments"
        --scaffold-log "${OUT_DIR}/scaffold_alignments"
        --outdir "$SPADES_DIR"
        -t "$THREADS"
        --spades-timeout "$SPADES_TIMEOUT"
    )
    [[ -n "$FASTA" ]] && REFINE_ARGS+=(--fasta "$FASTA")
    [[ "$CLEAN" == "true" ]] && REFINE_ARGS+=(--clean)
    [[ "$KEEP_INTERMEDIATES" == "true" ]] && REFINE_ARGS+=(--keep-intermediates)

    # refine.py resolves FASTQ paths as ./fastq/ relative to CWD; must match
    # where collect.py wrote them (OUT_DIR)
    pushd "$OUT_DIR" > /dev/null
    if ! python "${SCRIPT_DIR}/refine.py" \
            "${REFINE_ARGS[@]}" \
            > "${OUT_DIR}/refine.stdout.log" \
            2> "${OUT_DIR}/refine.stderr.log"; then
        popd > /dev/null
        err "$SAMPLE: refine.py failed — see ${OUT_DIR}/refine.stderr.log"
        return 1
    fi
    popd > /dev/null

    if [[ "$CLEAN" == "true" && "$KEEP_INTERMEDIATES" != "true" ]]; then
        rm -rf "${OUT_DIR}/fastq" "${SPADES_DIR}" "${OUT_DIR}/sv_summaries"
    fi

    touch "$DONE_FLAG"
    log "$SAMPLE: done → ${OUT_TABLE_STEM}.tsv"
}

# ── main ───────────────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_BASE"

if [[ -n "$SINGLE_SAMPLE" ]]; then
    run_sample "$SINGLE_SAMPLE"
    exit $?
fi

# Discover all unique sample names from SV summary files in the directory
mapfile -t ALL_SAMPLES < <(
    find "$SV_SUMMARY_DIR" -maxdepth 1 -name "*_SV_summary.tsv" -printf '%f\n' \
    | sed 's/_amplicon[0-9]*_SV_summary\.tsv//; s/_SV_summary\.tsv//' \
    | sort -u
)

TOTAL=${#ALL_SAMPLES[@]}
COUNT=0
FAILED=()

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    COUNT=$((COUNT + 1))
    log "[$COUNT/$TOTAL] $SAMPLE"
    if ! run_sample "$SAMPLE"; then
        FAILED+=("$SAMPLE")
    fi
done

log "Finished: $((TOTAL - ${#FAILED[@]}))/$TOTAL succeeded."
if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "Failed samples (${#FAILED[@]}):"
    printf '  %s\n' "${FAILED[@]}"
    exit 1
fi
