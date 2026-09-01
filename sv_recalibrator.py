#!/usr/bin/env python3
"""
sv_recalibrator.py — Unified CLI entry point for SVRecalibrator pipeline.

Subcommands:
    run       Run end-to-end breakpoint refinement for a single sample.
    batch     Orchestrate batch execution across a sample list.
    summarize Summarize microhomology gain and produce reporting artifacts.
"""

__version__ = "0.3.0"

import argparse
import os
import sys
import shutil
import subprocess

# Ensure active python environment bin directory is on PATH for samtools/spades.py
env_bin = os.path.dirname(sys.executable)
if env_bin not in os.environ.get("PATH", "").split(os.pathsep):
    os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")

import collect
import refine
import summarize_homology


def run_sample_cli(args):
    """Execute collect -> refine workflow for a single sample."""
    bam = os.path.abspath(args.bam)
    sv_dir = os.path.abspath(args.sv_dir)
    outdir = os.path.abspath(args.outdir)
    fasta = os.path.abspath(args.fasta) if args.fasta else None

    if not os.path.isfile(bam):
        sys.exit(f"Error: BAM file not found: {bam}")
    if not os.path.exists(sv_dir):
        sys.exit(f"Error: SV summary directory or file not found: {sv_dir}")
    if args.mode in ("scaffold", "both") and not fasta:
        sys.exit("Error: --fasta is required for scaffold or both mode")
    if fasta and not os.path.isfile(fasta):
        sys.exit(f"Error: FASTA reference not found: {fasta}")

    sample = args.sample if args.sample else os.path.basename(bam).split(".")[0]
    os.makedirs(outdir, exist_ok=True)

    alignments_tsv = os.path.join(outdir, "alignments.tsv")
    out_table_stem = os.path.join(outdir, "final_augmented")
    spades_dir = os.path.join(outdir, "spades")
    done_flag = os.path.join(outdir, "done.flag")

    print(f"[sv-recalibrator] Starting sample: {sample}")
    print(f"[sv-recalibrator] Output directory: {outdir}")

    cwd_orig = os.getcwd()
    try:
        # Step 1: collect reads around breakpoints
        os.chdir(outdir)

        collect_args = [
            str(args.radius),
            sv_dir,
            bam,
            "-s",
            sample,
            "-f",
            alignments_tsv,
        ]
        if args.strict:
            collect_args.append("--strict")
        if args.verbose:
            collect_args.append("-v")

        print(f"[sv-recalibrator] Running collect step...")
        collect.main(collect_args)

        if not os.path.isfile(alignments_tsv):
            sys.exit(f"Error: collection failed, {alignments_tsv} not generated.")

        # Step 2: refine breakpoints
        refine_args = [
            alignments_tsv,
            "--mode",
            args.mode,
            "--out-table",
            out_table_stem,
            "--split-log",
            os.path.join(outdir, "split_read_alignments"),
            "--scaffold-log",
            os.path.join(outdir, "scaffold_alignments"),
            "--outdir",
            spades_dir,
            "-t",
            str(args.threads),
            "--spades-timeout",
            str(args.spades_timeout),
        ]
        if fasta:
            refine_args.extend(["--fasta", fasta])
        if not args.keep_intermediates:
            refine_args.append("--clean")
        else:
            refine_args.append("--keep-intermediates")
        if args.verbose:
            refine_args.append("-v")

        print(f"[sv-recalibrator] Running refine step (mode={args.mode})...")
        refine.main(refine_args)

        # Cleanup if required
        if not args.keep_intermediates:
            print("[sv-recalibrator] Cleaning up intermediate FASTQ and SPAdes files...")
            shutil.rmtree(os.path.join(outdir, "fastq"), ignore_errors=True)
            shutil.rmtree(spades_dir, ignore_errors=True)

        with open(done_flag, "w") as f:
            f.write("done\n")

        print(f"[sv-recalibrator] Done! Output table: {out_table_stem}.tsv")

    finally:
        os.chdir(cwd_orig)


def run_batch_cli(args):
    """Delegate batch execution to batch_run.sh script."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    batch_script = os.path.join(script_dir, "batch_run.sh")

    if not os.path.isfile(batch_script):
        sys.exit(f"Error: batch_run.sh not found at {batch_script}")

    cmd = [
        "bash",
        batch_script,
        "-s",
        args.samples,
        "-v",
        args.sv_dir,
        "-o",
        args.outdir,
        "-m",
        args.mode,
        "-r",
        str(args.radius),
        "-t",
        str(args.threads),
        "--spades-timeout",
        str(args.spades_timeout),
    ]
    if args.fasta:
        cmd.extend(["-f", args.fasta])
    if args.strict:
        cmd.append("--strict")
    if args.sample:
        cmd.extend(["--sample", args.sample])
    if not args.keep_intermediates:
        cmd.append("--clean")
    else:
        cmd.append("--keep-intermediates")

    print(f"[sv-recalibrator] Launching batch pipeline: {' '.join(cmd)}")
    res = subprocess.run(cmd)
    sys.exit(res.returncode)


def run_summarize_cli(args):
    """Delegate to summarize_homology."""
    summ_args = [args.outdir]
    if args.sample:
        summ_args.extend(["--sample"] + args.sample)
    if args.min_hom_len:
        summ_args.extend(["--min-hom-len", str(args.min_hom_len)])
    if args.output_prefix:
        summ_args.extend(["-o", args.output_prefix])

    summarize_homology.main(summ_args)


def main():
    parser = argparse.ArgumentParser(
        description="SVRecalibrator — Structural Variant Breakpoint Refinement Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # ── Subcommand: run ────────────────────────────────────────────────────────
    p_run = subparsers.add_parser(
        "run", help="Run end-to-end refinement for a single sample"
    )
    p_run.add_argument(
        "--bam", required=True, help="Path to sample coordinate-sorted BAM file"
    )
    p_run.add_argument(
        "--sv-dir",
        required=True,
        help="Path to SV summaries directory (or summary file)",
    )
    p_run.add_argument(
        "-s",
        "--sample",
        help="Sample name/prefix (derived from BAM filename if omitted)",
    )
    p_run.add_argument(
        "-o", "--outdir", required=True, help="Output directory path"
    )
    p_run.add_argument(
        "-f",
        "--fasta",
        help="Indexed reference FASTA (required for scaffold or both mode)",
    )
    p_run.add_argument(
        "-m",
        "--mode",
        choices=["split", "scaffold", "both"],
        default="both",
        help="Refinement mode (default: both)",
    )
    p_run.add_argument(
        "-r",
        "--radius",
        type=int,
        default=350,
        help="Refinement radius around breakpoints in bp (default: 350)",
    )
    p_run.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="SPAdes assembly threads (default: 16)",
    )
    p_run.add_argument(
        "--spades-timeout",
        type=float,
        default=2.0,
        help="SPAdes timeout per breakpoint in hours; 0=no limit (default: 2.0)",
    )
    p_run.add_argument(
        "--strict",
        action="store_true",
        help="Keep only reads fully aligning within region of interest",
    )
    p_run.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Preserve temporary FASTQ and SPAdes assembly files (default: auto-clean)",
    )
    p_run.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose log output"
    )

    # ── Subcommand: batch ──────────────────────────────────────────────────────
    p_batch = subparsers.add_parser(
        "batch", help="Run batch refinement across multiple samples"
    )
    p_batch.add_argument(
        "-s",
        "--samples",
        required=True,
        help="Tab-separated sample list file (SAMPLE BAM_PATH)",
    )
    p_batch.add_argument(
        "-v",
        "--sv-dir",
        required=True,
        help="Shared SV summary directory containing sample TSVs",
    )
    p_batch.add_argument(
        "-o", "--outdir", required=True, help="Batch output root directory"
    )
    p_batch.add_argument(
        "-f",
        "--fasta",
        help="Indexed reference FASTA (required for scaffold/both mode)",
    )
    p_batch.add_argument(
        "-m",
        "--mode",
        choices=["split", "scaffold", "both"],
        default="both",
        help="Refinement mode (default: both)",
    )
    p_batch.add_argument(
        "-r",
        "--radius",
        type=int,
        default=350,
        help="Refinement radius in bp (default: 350)",
    )
    p_batch.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="SPAdes threads per sample (default: 16)",
    )
    p_batch.add_argument(
        "--spades-timeout",
        type=float,
        default=2.0,
        help="Per-breakpoint SPAdes timeout in hours (default: 2.0)",
    )
    p_batch.add_argument(
        "--strict", action="store_true", help="Filter for strict region overlap"
    )
    p_batch.add_argument(
        "--sample", help="Run a single named sample from the sample list"
    )
    p_batch.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Preserve intermediate FASTQ and SPAdes files (default: auto-clean)",
    )

    # ── Subcommand: summarize ──────────────────────────────────────────────────
    p_summ = subparsers.add_parser(
        "summarize", help="Summarize homology improvement over baseline AA calls"
    )
    p_summ.add_argument(
        "outdir", help="Batch output directory containing sample results"
    )
    p_summ.add_argument(
        "--sample", nargs="+", help="Restrict summary to specified sample(s)"
    )
    p_summ.add_argument(
        "--min-hom-len",
        type=int,
        default=1,
        help="Minimum homology length to count (default: 1)",
    )
    p_summ.add_argument(
        "-o",
        "--output-prefix",
        help="Prefix for summary CSV and visualization plots",
    )

    args = parser.parse_args()

    if not args.subcommand:
        parser.print_help()
        sys.exit(1)

    if args.subcommand == "run":
        run_sample_cli(args)
    elif args.subcommand == "batch":
        run_batch_cli(args)
    elif args.subcommand == "summarize":
        run_summarize_cli(args)


if __name__ == "__main__":
    main()
