#!/usr/bin/env python3
"""
summarize_homology.py — Summarize homology-detection improvement over AA baseline.

Loads final_augmented.tsv files produced by SVRecalibrator (refine.py --mode both)
from a batch output directory and reports how many SVs gained microhomology by
split-read and/or scaffold approaches that were not present in the AA input calls.

Usage:
    python summarize_homology.py <outdir> [--sample SAMPLE ...] [--min-hom-len N]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── helpers ────────────────────────────────────────────────────────────────────

def parse_sp_hom_len(val):
    """Convert sp_hom_len (stored as string) to int or NaN."""
    if pd.isna(val) or val == "" or val == "nan":
        return np.nan
    if str(val).strip() == "N/A":
        return 0
    try:
        return int(val)
    except (ValueError, TypeError):
        return np.nan


def load_augmented(outdir, samples=None):
    """Load all final_augmented.tsv files under outdir/*/final_augmented.tsv."""
    frames = []
    missing = []

    candidate_dirs = (
        [os.path.join(outdir, s) for s in samples]
        if samples
        else [
            os.path.join(outdir, d)
            for d in sorted(os.listdir(outdir))
            if os.path.isdir(os.path.join(outdir, d))
        ]
    )

    for sample_dir in candidate_dirs:
        sample = os.path.basename(sample_dir)
        tsv = os.path.join(sample_dir, "final_augmented.tsv")
        if not os.path.isfile(tsv):
            missing.append(sample)
            continue
        try:
            df = pd.read_csv(tsv, sep="\t", low_memory=False)
        except EmptyDataError:
            continue  # empty file — sample produced no SVs
        except Exception as e:
            print(f"Warning: could not load {tsv}: {e}", file=sys.stderr)
            missing.append(sample)
            continue
        if df.empty:
            continue
        df["_sample"] = sample
        frames.append(df)

    if missing:
        print(f"Warning: no final_augmented.tsv for {len(missing)} sample(s): {', '.join(missing)}",
              file=sys.stderr)
    if not frames:
        sys.exit("No data loaded — check outdir and sample names.")
    return pd.concat(frames, ignore_index=True)


def classify(df, min_hom_len):
    """Add boolean columns for AA / split-read / scaffold homology presence."""
    # AA homology (optional column)
    if "homology_len" in df.columns:
        df["aa_hom"] = pd.to_numeric(df["homology_len"], errors="coerce").fillna(0) >= min_hom_len
    else:
        df["aa_hom"] = False

    # Split-read homology: sp_hom_len is a string; positive = homology
    df["_sp_hom_len_int"] = df["sp_hom_len"].apply(parse_sp_hom_len)
    df["sp_hom"] = df["_sp_hom_len_int"] >= min_hom_len

    # Scaffold homology: sc_hom_len is numeric; positive = homology
    df["_sc_hom_len_num"] = pd.to_numeric(df["sc_hom_len"], errors="coerce")
    df["sc_hom"] = df["_sc_hom_len_num"] >= min_hom_len

    df["either_hom"] = df["sp_hom"] | df["sc_hom"]

    # Insertion flags (negative lengths)
    if "homology_len" in df.columns:
        aa_num = pd.to_numeric(df["homology_len"], errors="coerce")
        df["aa_ins"]   = aa_num < 0
        df["aa_blunt"] = aa_num == 0
    else:
        df["aa_ins"]   = False
        df["aa_blunt"] = False

    df["sv_ins"]   = (df["_sp_hom_len_int"] < 0) | (df["_sc_hom_len_num"] < 0)
    df["sv_blunt"] = (df["_sp_hom_len_int"] == 0) | (df["_sc_hom_len_num"] == 0)

    df["aa_any"] = df["aa_hom"] | df["aa_ins"] | df["aa_blunt"]
    df["sv_any"] = df["either_hom"] | df["sv_ins"] | df["sv_blunt"]

    return df


def pct(n, total):
    return f"{n:>5d}  ({100 * n / total:.1f}%)" if total else f"{n:>5d}  (  n/a)"


def print_section(title):
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print(f"{'─' * 60}")


def summarize(df, min_hom_len):
    total = len(df)
    has_aa = "homology_len" in df.columns

    print_section("Dataset overview")
    samples = df["_sample"].nunique()
    print(f"  Samples loaded           : {samples}")
    print(f"  Total SVs                : {total}")
    if has_aa:
        n_aa = df["aa_hom"].sum()
        print(f"  SVs with AA homology     : {pct(n_aa, total)}")
        print(f"  SVs without AA homology  : {pct(total - n_aa, total)}")
    else:
        print("  (homology_len column absent — AA homology comparison unavailable)")

    # ── SVRecalibrator detection across ALL SVs ────────────────────────────────
    print_section("SVRecalibrator homology detection (all SVs)")
    n_sp   = df["sp_hom"].sum()
    n_sc   = df["sc_hom"].sum()
    n_both = (df["sp_hom"] & df["sc_hom"]).sum()
    n_any  = df["either_hom"].sum()
    n_sp_only = (df["sp_hom"] & ~df["sc_hom"]).sum()
    n_sc_only = (~df["sp_hom"] & df["sc_hom"]).sum()

    print(f"  Detected by split reads  : {pct(n_sp,   total)}")
    print(f"  Detected by scaffold     : {pct(n_sc,   total)}")
    print(f"  Detected by both         : {pct(n_both, total)}")
    print(f"  Detected by either       : {pct(n_any,  total)}")
    print(f"    split only             : {pct(n_sp_only, total)}")
    print(f"    scaffold only          : {pct(n_sc_only, total)}")

    # ── Gain over AA ──────────────────────────────────────────────────────────
    if has_aa:
        print_section(f"Homology gain over AA (SVs without AA homology, min_hom_len={min_hom_len})")
        no_aa = df[~df["aa_hom"]]
        N = len(no_aa)
        print(f"  Universe (no AA hom)     : {N}")
        if N:
            g_sp      = no_aa["sp_hom"].sum()
            g_sc      = no_aa["sc_hom"].sum()
            g_both    = (no_aa["sp_hom"] & no_aa["sc_hom"]).sum()
            g_any     = no_aa["either_hom"].sum()
            g_sp_only = (no_aa["sp_hom"] & ~no_aa["sc_hom"]).sum()
            g_sc_only = (~no_aa["sp_hom"] & no_aa["sc_hom"]).sum()
            g_none    = (~no_aa["either_hom"]).sum()

            print(f"  Gained by either         : {pct(g_any,     N)}")
            print(f"    split reads only       : {pct(g_sp_only, N)}")
            print(f"    scaffold only          : {pct(g_sc_only, N)}")
            print(f"    both approaches        : {pct(g_both,    N)}")
            print(f"  Still no homology        : {pct(g_none,    N)}")

        # Concordance on SVs that DO have AA homology
        print_section("Concordance on SVs with AA homology")
        with_aa = df[df["aa_hom"]]
        M = len(with_aa)
        print(f"  Universe (AA hom)        : {M}")
        if M:
            c_any  = with_aa["either_hom"].sum()
            c_none = M - c_any
            print(f"  Confirmed by SVRecalib   : {pct(c_any,  M)}")
            print(f"  Not confirmed            : {pct(c_none, M)}")

    # ── Homology length distributions ─────────────────────────────────────────
    print_section("Homology length distributions (detected SVs only)")

    def len_stats(series, label):
        s = pd.to_numeric(series, errors="coerce").dropna()
        s = s[s > 0]
        if s.empty:
            print(f"  {label}: no data")
            return
        print(f"  {label} (n={len(s)}):")
        print(f"    median={s.median():.0f}  mean={s.mean():.1f}  "
              f"min={s.min():.0f}  max={s.max():.0f}  "
              f"p25={s.quantile(0.25):.0f}  p75={s.quantile(0.75):.0f}")

    len_stats(df.loc[df["sp_hom"], "_sp_hom_len_int"], "split-read hom len")
    len_stats(df.loc[df["sc_hom"], "_sc_hom_len_num"], "scaffold   hom len")
    if has_aa:
        len_stats(df.loc[df["aa_hom"], "homology_len"], "AA         hom len")

    # ── Per-sample breakdown ───────────────────────────────────────────────────
    print_section("Per-sample breakdown")
    grp = df.groupby("_sample")
    rows = []
    for sample, g in grp:
        n = len(g)
        aa  = g["aa_hom"].sum() if has_aa else np.nan
        sp  = g["sp_hom"].sum()
        sc  = g["sc_hom"].sum()
        any_ = g["either_hom"].sum()
        no_aa_g = g[~g["aa_hom"]] if has_aa else g
        gain = no_aa_g["either_hom"].sum() if has_aa else np.nan
        rows.append({"sample": sample, "SVs": n, "AA_hom": aa,
                     "sp_hom": sp, "sc_hom": sc, "either": any_, "new_gain": gain})
    tbl = pd.DataFrame(rows).set_index("sample")
    if not has_aa:
        tbl = tbl.drop(columns=["AA_hom", "new_gain"])
    print(tbl.to_string())
    print()
    return tbl


# ── plotting ──────────────────────────────────────────────────────────────────

def _circle_intersection_area(r1, r2, d):
    """Intersection area of two circles with radii r1, r2 and centre distance d."""
    if d >= r1 + r2:
        return 0.0
    if d + min(r1, r2) <= max(r1, r2):
        return np.pi * min(r1, r2) ** 2
    a = np.arccos(np.clip((d**2 + r1**2 - r2**2) / (2 * d * r1), -1, 1))
    b = np.arccos(np.clip((d**2 + r2**2 - r1**2) / (2 * d * r2), -1, 1))
    c = 0.5 * np.sqrt(max((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2), 0))
    return r1**2 * a + r2**2 * b - c


def _find_centre_distance(r1, r2, target_area):
    """Binary-search for the centre distance that yields target_area overlap."""
    if target_area <= 0:
        return r1 + r2
    max_area = np.pi * min(r1, r2) ** 2
    if target_area >= max_area:
        return abs(r1 - r2)
    lo, hi = float(abs(r1 - r2)), float(r1 + r2)
    for _ in range(64):
        mid = (lo + hi) / 2
        if _circle_intersection_area(r1, r2, mid) > target_area:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _draw_venn(ax, sets, labels, colors, neither_label=None):
    """
    Proportional Euler diagram: rectangle + two circles with geometrically
    correct areas and overlap.
    sets  : (total, n_left, n_right, n_both)
    labels: (title, left_label, right_label)
    colors: (left_color, right_color)
    """
    total, n_left, n_right, n_both = sets
    title, left_label, right_label = labels
    lc, rc = colors

    W, H   = 10.0, 7.0
    pad    = 0.45
    rx, ry = pad, pad
    rw, rh = W - 2 * pad, H - 2 * pad

    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")

    # Rectangle — all SVs
    ax.add_patch(mpatches.FancyBboxPatch(
        (rx, ry), rw, rh,
        boxstyle="round,pad=0.05",
        linewidth=1.5, edgecolor="#444", facecolor="whitesmoke", zorder=0,
    ))

    # Scale: map SV count → display area so areas are truly proportional
    k = (rw * rh) / total
    r_l = np.sqrt(n_left  * k / np.pi) if n_left  else 0.0
    r_r = np.sqrt(n_right * k / np.pi) if n_right else 0.0

    # Clamp so neither circle overflows the rectangle
    max_r = min(rw, rh) / 2 * 0.97
    sf = min(1.0, max_r / max(r_l, r_r, 1e-9))
    r_l *= sf
    r_r *= sf

    # Solve for the centre distance that gives the correct overlap area
    target_overlap = n_both * k * sf ** 2
    d = _find_centre_distance(r_l, r_r, target_overlap)

    # Centre both circles symmetrically about the rectangle centre
    cx = rx + rw / 2
    cy = ry + rh / 2
    cx_l = cx - d / 2
    cx_r = cx + d / 2

    ax.add_patch(plt.Circle((cx_l, cy), r_l, color=lc, alpha=0.40, zorder=2))
    ax.add_patch(plt.Circle((cx_r, cy), r_r, color=rc, alpha=0.40, zorder=2))

    # ── region counts ──────────────────────────────────────────────────────────
    n_l_only  = n_left  - n_both
    n_r_only  = n_right - n_both
    n_neither = total   - n_left - n_right + n_both

    kw = dict(ha="center", va="center", fontsize=11, fontweight="bold", zorder=4)

    # "both" — centre of the overlap zone
    ax.text((cx_l + cx_r) / 2, cy, str(n_both), **kw)

    # right-only — well inside right crescent
    ax.text(cx_r + r_r * 0.58, cy, str(n_r_only), **kw)

    # left-only — annotate outside with arrow if region is too small to label in place
    if n_l_only / total < 0.02:
        ax.annotate(
            str(n_l_only),
            xy=(cx_l - r_l * 0.92, cy),
            xytext=(rx + 0.3, ry + rh * 0.82),
            fontsize=10, fontweight="bold", color=lc, zorder=4,
            arrowprops=dict(arrowstyle="-|>", color=lc, lw=0.9),
        )
    else:
        ax.text(cx_l - r_l * 0.58, cy, str(n_l_only), **kw)

    # neither — bottom-left; emphasised when a label is provided
    if neither_label:
        ax.text(rx + 0.3, ry + 0.3, str(n_neither),
                fontsize=14, fontweight="bold", color="#444", zorder=4, va="bottom")
        ax.text(rx + 0.3, ry + 1.05, neither_label,
                fontsize=9, color="#444", zorder=4, va="bottom", style="italic")
    else:
        ax.text(rx + 0.3, ry + 0.35, str(n_neither),
                fontsize=9, color="dimgray", zorder=4, va="bottom")

    # total label — top-left
    ax.text(rx + 0.3, ry + rh - 0.25, f"n={total} total SVs",
            fontsize=9, va="top", color="#444", zorder=4)

    # circle labels as a legend box — avoids placement issues with nested circles
    legend_patches = [
        mpatches.Patch(facecolor=lc, alpha=0.55, label=left_label),
        mpatches.Patch(facecolor=rc, alpha=0.55, label=right_label),
    ]
    ax.legend(handles=legend_patches, loc="lower right",
              bbox_to_anchor=(0.97, 0.03), framealpha=0.88,
              fontsize=9, edgecolor="#aaa")

    ax.set_title(title, fontsize=11, pad=6)


def _save_fig(fig, prefix, suffix):
    for ext in (".png", ".pdf"):
        path = prefix + suffix + ext
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Written: {path}")
    plt.close(fig)


def _hist_series(df):
    """Return ordered (values, label, color) tuples for the three histogram tracks."""
    has_aa = "homology_len" in df.columns
    aa_has_data = has_aa and pd.to_numeric(df["homology_len"], errors="coerce").notna().any()
    sp_vals = df["_sp_hom_len_int"].dropna()
    sc_vals = df["_sc_hom_len_num"].dropna()
    series = []
    if aa_has_data:
        aa_vals = pd.to_numeric(df["homology_len"], errors="coerce").dropna()
        series.append((aa_vals, "AA", "#D65F5F"))
    series.append((sp_vals, "Split reads", "#4878CF"))
    series.append((sc_vals, "Scaffold",    "#6ACC65"))
    return series


def plot_venn(df, prefix):
    has_aa = "homology_len" in df.columns
    total  = len(df)

    if has_aa:
        # 1×3: homology / insertion / blunt
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
        fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.02, wspace=0.06)
        panels = [
            ("Homology (> 0)",  "aa_hom",   "either_hom", "#D65F5F"),
            ("Insertion (< 0)", "aa_ins",   "sv_ins",     "#E8A838"),
            ("Blunt (= 0)",     "aa_blunt", "sv_blunt",   "#6ACC65"),
        ]
        for ax, (title, aa_col, sv_col, rc) in zip(axes, panels):
            n_aa   = int(df[aa_col].sum())
            n_sv   = int(df[sv_col].sum())
            n_both = int((df[aa_col] & df[sv_col]).sum())
            _draw_venn(ax, (total, n_aa, n_sv, n_both),
                       (title, "AA", "SVRecalibrator"),
                       ("#4878CF", rc))
    else:
        # Fallback: single panel, split reads vs scaffold (homology only)
        fig, ax = plt.subplots(figsize=(6, 5.5))
        fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.02)
        n_sp   = int(df["sp_hom"].sum())
        n_sc   = int(df["sc_hom"].sum())
        n_both = int((df["sp_hom"] & df["sc_hom"]).sum())
        _draw_venn(ax, (total, n_sp, n_sc, n_both),
                   ("Homology (> 0)", "Split reads", "Scaffold"),
                   ("#4878CF", "#6ACC65"))

    fig.suptitle("SVRecalibrator — junction type: AA vs SVRecalibrator",
                 fontsize=13, fontweight="bold")
    _save_fig(fig, prefix, "_venn")


def plot_venn_combined(df, prefix):
    total = len(df)
    fig, ax = plt.subplots(figsize=(6, 5.5))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.02)

    n_aa   = int(df["aa_any"].sum())
    n_sv   = int(df["sv_any"].sum())
    n_both = int((df["aa_any"] & df["sv_any"]).sum())
    _draw_venn(ax, (total, n_aa, n_sv, n_both),
               ("Any junction call", "AA", "SVRecalibrator"),
               ("#4878CF", "#9B59B6"), neither_label="unrefined")

    fig.suptitle("SVRecalibrator — junction refinement overview",
                 fontsize=12, fontweight="bold")
    _save_fig(fig, prefix, "_venn_combined")


def plot_hist(df, prefix):
    series   = _hist_series(df)
    n_tracks = len(series)

    # x range: hard-capped at ±50; values outside pile into the endpoint bins
    lo, hi = -50, 50
    bins = np.arange(lo - 0.5, hi + 1.5, 1.0)

    fig, axes_h = plt.subplots(n_tracks, 1, figsize=(7, 2.5 * n_tracks + 0.8),
                                sharex=True)
    if n_tracks == 1:
        axes_h = [axes_h]
    fig.subplots_adjust(left=0.11, right=0.97, top=0.91, bottom=0.10, hspace=0.08)

    for i, (vals, label, color) in enumerate(series):
        ax = axes_h[i]
        n_lo = int((vals < lo).sum())
        n_hi = int((vals > hi).sum())
        ax.hist(vals.clip(lo, hi), bins=bins, color=color, alpha=0.78,
                edgecolor="white", linewidth=0.3)
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--", zorder=3)
        ax.set_ylabel("Count", fontsize=8)
        ax.text(0.01, 0.93, label, transform=ax.transAxes,
                ha="left", va="top", fontsize=9, fontweight="bold", color=color,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.6))
        notes = []
        if n_lo:
            notes.append(f"+{n_lo} in ≤−50 bin")
        if n_hi:
            notes.append(f"+{n_hi} in ≥50 bin")
        if notes:
            ax.text(0.99, 0.93, "  ".join(notes),
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=7, color="gray")
        if i == 0:
            ax.axvspan(lo - 0.5, 0, alpha=0.06, color="orange", zorder=0)
            ax.axvspan(0, hi + 0.5, alpha=0.06, color="steelblue", zorder=0)
            ax.legend(handles=[
                mpatches.Patch(facecolor="orange",    alpha=0.5, label="insertion (< 0)"),
                mpatches.Patch(facecolor="steelblue", alpha=0.5, label="homology (> 0)"),
            ], fontsize=7, loc="upper right", framealpha=0.8)
        if i < n_tracks - 1:
            ax.tick_params(labelbottom=False)
        ax.tick_params(axis="both", labelsize=8)

    # Uniform y-limits
    max_y = max(ax.get_ylim()[1] for ax in axes_h)
    for ax in axes_h:
        ax.set_ylim(0, max_y)

    # Label the ≤−50 and ≥50 endpoint ticks on the bottom axis
    bot = axes_h[-1]
    from matplotlib.ticker import FixedLocator
    ticks = sorted(set(list(bot.get_xticks()) + [lo, hi]))
    bot.xaxis.set_major_locator(FixedLocator(ticks))
    bot.set_xticklabels(
        ["≤−50" if t == lo else "≥50" if t == hi else str(int(t)) for t in ticks],
        fontsize=8,
    )

    axes_h[0].set_title("Junction length distributions", fontsize=10)
    axes_h[-1].set_xlabel("Junction length (bp)     ← insertion  |  homology →", fontsize=9)

    fig.suptitle("SVRecalibrator — microhomology summary", fontsize=12, fontweight="bold")
    _save_fig(fig, prefix, "_hist")


# ── CLI ────────────────────────────────────────────────────────────────────────

def main(args_list=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir", help="batch output directory (contains one subdirectory per sample)")
    p.add_argument("--sample", nargs="+", metavar="NAME", help="restrict to these sample(s)")
    p.add_argument(
        "--min-hom-len",
        type=int,
        default=1,
        metavar="N",
        help="minimum homology length to count as detected (default: 1)",
    )
    p.add_argument(
        "-o", "--output-prefix",
        metavar="PREFIX",
        help="output prefix; produces <PREFIX>_summary.csv, <PREFIX>_venn.png, <PREFIX>_hist.png",
    )
    args = p.parse_args(args_list)

    if not os.path.isdir(args.outdir):
        sys.exit(f"outdir not found: {args.outdir}")

    df = load_augmented(args.outdir, args.sample)
    df = classify(df, args.min_hom_len)
    tbl = summarize(df, args.min_hom_len)
    if args.output_prefix:
        tbl.to_csv(args.output_prefix + "_summary.csv")
        print(f"Written: {args.output_prefix}_summary.csv")
        plot_venn(df, args.output_prefix)
        plot_venn_combined(df, args.output_prefix)
        plot_hist(df, args.output_prefix)


if __name__ == "__main__":
    main()

