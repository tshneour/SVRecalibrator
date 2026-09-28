#!/usr/bin/env python3
"""Compare junction homology/insertion calls: AA, SVRecalibrator, Delly2, SvABA, GRIDSS2.

Every SVRecalibrator row is matched to each caller's nearest orientation-compatible
breakpoint pair (both original AA ends within --window, ends may swap). Among
compatible calls, a call with a resolved junction is preferred over a nearer
unresolved one, giving each caller its best chance at the junction.

Junction length is signed: homology > 0, blunt 0, insertion < 0. Per tool:
  AA              homology_len; AA writes None when no split reads were found
                  (empty after SVRecalibrator's pandas round trip) -> no junction call
  SVRecalibrator  split-read (sp_hom_len) and scaffold (sc_hom_len) lengths. Empty
                  means that path made no call (split: no supported candidate;
                  scaffold: SPAdes failure/timeout or no valid alignment). Split
                  'N/A' (insertion chosen without a sequence) is also no call. One
                  path alone is used; rows where both exist but disagree get their
                  own category (they were dropped from earlier figures).
  Delly2          reconstructed from CONSENSUS (delly_junctions.py); INSLEN is always 0
  SvABA           INSERTION / HOMSEQ; discordant-only (DSCRD) calls have no junction
  GRIDSS2         inserted bases in the BND ALT / HOMSEQ; IMPRECISE calls have no junction
Caller outcomes per SV: junction class, 'sv_called_no_junction', or 'sv_not_called'.
The main comparison is restricted to SVs identified (called at all) by every caller.
All FILTER values are retained by default; --pass-only restricts every caller to PASS.
Inputs are either directories or per-tool tarballs built by bundle_caller_outputs.py
(--bundle DIR; the newest tarball per tool is used, default DIR turner2017/ next to
this script).
"""
import argparse
import csv
import json
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

from caller_vcfs import read_delly, read_gridss, read_svaba
from consolidate_delly import write_csv
from delly_junctions import make_aligner, reconstruct

CALLERS = ['delly', 'svaba', 'gridss']
BUNDLE_DIR = Path(__file__).resolve().parent / 'turner2017'
LABELS = {'aa': 'AA', 'svr': 'SVRecalibrator', 'delly': 'Delly2 (reconstructed)',
          'svaba': 'SvABA', 'gridss': 'GRIDSS2'}
COLORS = {'aa': '#D65F5F', 'svr': '#4878CF', 'delly': '#6ACC65', 'svaba': '#B47CC7', 'gridss': '#C4AD66'}


def number(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def svr_length(row):
    """Return (length, status) for the combined SVRecalibrator call."""
    sp = number(row.get('sp_hom_len'))
    sc = number(row.get('sc_hom_len'))
    if sp is not None and sc is not None and sp != sc:
        return None, 'split_scaffold_disagree'
    value = sp if sp is not None else sc
    return value, 'resolved' if value is not None else 'no_junction_call'


def chrom(value):
    return value.removeprefix('chr')


class CallIndex:
    def __init__(self, calls):
        self.by_pair = defaultdict(list)
        for call in calls:
            self.by_pair[(chrom(call['chrom1']), chrom(call['chrom2']))].append(call)

    def near(self, c1, p1, c2, p2, orientation, window):
        """Yield (call, swapped, d1, d2) for orientation-compatible calls near both ends."""
        for swapped, pair in ((False, (c1, c2)), (True, (c2, c1))):
            for call in self.by_pair.get(pair, []):
                a, b = (call['pos2'], call['pos1']) if swapped else (call['pos1'], call['pos2'])
                d1, d2 = a - p1, b - p2
                if max(abs(d1), abs(d2)) > window:
                    continue
                ori = call['orientation'][::-1] if swapped else call['orientation']
                if ori == orientation:
                    yield call, swapped, d1, d2


def resolve_delly(call, fasta, aligner, cache):
    if call['key'] not in cache:
        if not call['precise']:
            out = dict(dl_status='imprecise')
        else:
            out = reconstruct(call, fasta, aligner)
        cache[call['key']] = out
        call.update(out)
        call['resolved'] = out['dl_status'] == 'resolved'
        call['junction_len'] = out.get('dl_junction_len') if call['resolved'] else None
        call['junction_seq'] = out.get('dl_junction_seq', '') if call['resolved'] else ''
    return call


def load_calls(sample, args):
    calls = {}
    if args.delly:
        path = args.delly / f'{sample}.vcf'
        path = path if path.exists() else args.delly / f'{sample}.vcf.gz'
        calls['delly'] = read_delly(path) if path.exists() else None
    if args.svaba:
        path = args.svaba / sample / f'{sample}.svaba.{args.svaba_vcf}.vcf'
        calls['svaba'] = read_svaba(path) if path.exists() else None
    if args.gridss:
        path = args.gridss / f'{sample}.vcf.gz'
        path = path if path.exists() else args.gridss / f'{sample}.vcf'
        calls['gridss'] = read_gridss(path) if path.exists() else None
    if args.pass_only:
        calls = {k: v and [c for c in v if c['filter'] == 'PASS'] for k, v in calls.items()}
    return calls


def compare(args):
    fasta = pysam.FastaFile(str(args.fasta))
    aligner, cache = make_aligner(), {}
    rows, inventory = [], []
    callers = [c for c in CALLERS if getattr(args, c)]
    for directory in sorted(p for p in args.sv_dir.iterdir() if p.is_dir()):
        sample = directory.name
        table = directory / 'final_augmented.tsv'
        records = []
        if table.exists():
            with table.open(newline='') as handle:
                records = list(csv.DictReader(handle, delimiter='\t'))
        calls = load_calls(sample, args)
        present = {c: calls[c] is not None for c in callers}
        inventory.append(dict(sample_id=sample, svr_rows=len(records),
                              **{f'{c}_vcf': present[c] for c in callers},
                              **{f'{c}_calls': len(calls[c] or []) for c in callers}))
        indexes = {c: CallIndex(calls[c]) for c in callers if present[c]}
        for n, record in enumerate(records, 1):
            out = dict(sample_id=sample, svr_row_id=f'{sample}:{n}', **record)
            out['aa_len'] = number(record.get('homology_len'))
            out['aa_status'] = 'resolved' if out['aa_len'] is not None else 'no_junction_call'
            out['svr_len'], out['svr_status'] = svr_length(record)
            c1, c2 = chrom(record['break_chrom1']), chrom(record['break_chrom2'])
            p1, p2 = int(record['break_pos1']), int(record['break_pos2'])
            for caller in callers:
                if not present[caller]:
                    out[f'{caller}_status'] = 'missing_vcf'
                    continue
                found = []
                for call, swapped, d1, d2 in indexes[caller].near(c1, p1, c2, p2, record['break_orientation'], args.window):
                    if caller == 'delly':
                        resolve_delly(call, fasta, aligner, cache)
                    found.append((not call['resolved'], max(abs(d1), abs(d2)), abs(d1) + abs(d2), call['key'],
                                  call, swapped, d1, d2))
                out[f'{caller}_candidates'] = len(found)
                if not found:
                    out[f'{caller}_status'] = 'sv_not_called'
                    continue
                *_, call, swapped, d1, d2 = min(found, key=lambda t: t[:4])
                out[f'{caller}_status'] = 'resolved' if call['resolved'] else 'sv_called_no_junction'
                out.update({f'{caller}_{k}': call.get(k, '') for k in
                            ['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'orientation', 'filter', 'qual',
                             'junction_len', 'junction_seq', 'evidence']})
                out.update({f'{caller}_ends_swapped': swapped, f'{caller}_delta_pos1': d1,
                            f'{caller}_delta_pos2': d2})
                if caller == 'delly':
                    out.update({f'delly_{k}': call.get(k, '') for k in
                                ['raw_homlen', 'raw_inslen', 'svtype', 'dl_status', 'dl_shift', 'dl_min_identity']})
            out['identified_by_all_callers'] = all(
                out[f'{c}_status'] in ('resolved', 'sv_called_no_junction') for c in callers)
            rows.append(out)
    return rows, inventory, callers


def junction_class(value):
    if value is None or value == '':
        return None
    return 'insertion' if value < 0 else 'blunt' if value == 0 else 'homology'


def track_length(row, track):
    if row.get(f'{track}_status') != 'resolved':
        return None
    return row[f'{track}_len'] if track in ('aa', 'svr') else row[f'{track}_junction_len']


def track_category(row, track):
    return junction_class(track_length(row, track)) or row.get(f'{track}_status')


def summarize(rows, callers):
    tracks = ['aa', 'svr'] + callers
    summary = {'rows': len(rows), 'samples': len({r['sample_id'] for r in rows})}
    for t in tracks:
        summary[t] = dict(Counter(track_category(r, t) for r in rows))
    agreement = {}
    for caller in callers:
        both = [(track_length(r, 'svr'), track_length(r, caller)) for r in rows]
        both = [(a, b) for a, b in both if a is not None and b is not None]
        agreement[caller] = dict(
            comparable=len(both), exact_length=sum(a == b for a, b in both),
            same_class=sum(junction_class(a) == junction_class(b) for a, b in both),
            svr_insertions=sum(a < 0 for a, _ in both),
            svr_insertions_called_insertion=sum(a < 0 and b < 0 for a, b in both))
    summary['agreement_with_svr'] = agreement
    resolved_all = [r for r in rows if all(track_length(r, t) is not None for t in ['svr'] + callers)]
    summary['junction_resolved_by_svr_and_all_callers'] = dict(
        rows=len(resolved_all),
        all_same_length=sum(len({track_length(r, t) for t in ['svr'] + callers}) == 1 for r in resolved_all),
        all_same_class=sum(len({junction_class(track_length(r, t)) for t in ['svr'] + callers}) == 1
                           for r in resolved_all))
    return summary


def track_values(rows, track):
    return np.array([v for v in (track_length(r, track) for r in rows) if v is not None])


def plot_hist(rows, callers, prefix, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator
    tracks = ['aa', 'svr'] + callers
    lo, hi = -50, 50
    bins = np.arange(lo - 0.5, hi + 1.5, 1.0)
    fig, axes = plt.subplots(len(tracks), 1, figsize=(10, 2.6 * len(tracks)), sharex=True)
    for ax, track in zip(axes, tracks):
        vals = track_values(rows, track)
        ax.hist(np.clip(vals, lo, hi), bins=bins, color=COLORS[track], alpha=0.78, edgecolor='white', linewidth=0.3)
        ax.axvline(0, color='black', linewidth=1, linestyle='--', zorder=3)
        ax.set_ylabel('SV count')
        ax.text(0.01, 0.93, f'{LABELS[track]} (n={len(vals)})', transform=ax.transAxes, ha='left', va='top',
                fontweight='bold', color=COLORS[track],
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#333', alpha=0.9))
        notes = [f'+{n} in {lab} bin' for n, lab in [((vals <= lo).sum(), '≤-50'), ((vals >= hi).sum(), '≥50')] if n]
        if notes:
            ax.text(0.99, 0.93, '  '.join(notes), transform=ax.transAxes, ha='right', va='top', fontsize=9, color='gray')
    top = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, top)
    ticks = [lo, -40, -20, 0, 20, 40, hi]
    axes[-1].xaxis.set_major_locator(FixedLocator(ticks))
    axes[-1].set_xticklabels(['≤−50' if t == lo else '≥50' if t == hi else str(t) for t in ticks])
    axes[-1].set_xlabel('Signed junction length (bp):  insertion (−)  |  blunt (0)  |  homology (+)')
    axes[0].set_title(title)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'{prefix}_hist.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_classes(summary, callers, prefix, title):
    """Stacked bars over the same SV rows: how each tool characterized every junction."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tracks = ['aa', 'svr'] + callers
    parts = [('homology', '#4C72B0', 'homology'), ('blunt', '#8C8C8C', 'blunt'),
             ('insertion', '#DD8452', 'insertion'),
             ('split_scaffold_disagree', '#F2D16B', 'split/scaffold disagree'),
             ('no_junction_call', '#E5E5E5', 'no junction call'),
             ('sv_called_no_junction', '#E5E5E5', None),
             ('sv_not_called', 'white', 'SV not called by caller')]
    fig, ax = plt.subplots(figsize=(9, 0.7 * len(tracks) + 1.5))
    for i, track in enumerate(reversed(tracks)):
        left = 0
        for part, color, label in parts:
            n = summary[track].get(part, 0)
            ax.barh(i, n, left=left, color=color, edgecolor='#555', linewidth=0.5,
                    label=label if i == 0 else None)
            if n >= summary['rows'] * 0.04:
                ax.text(left + n / 2, i, str(n), ha='center', va='center', fontsize=8,
                        color='white' if part in ('homology', 'insertion') else 'black')
            left += n
    ax.set_yticks(range(len(tracks)))
    ax.set_yticklabels([LABELS[t] for t in reversed(tracks)])
    ax.set_xlim(0, summary['rows'])
    ax.set_xlabel(f"SV junctions (n={summary['rows']}, {summary['samples']} samples)")
    ax.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.25), frameon=False, fontsize=8)
    ax.set_title(title)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'{prefix}_classes.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


# bundle tool name -> compare_callers input argument
BUNDLE_TOOLS = {'svr': 'sv_dir', 'delly2': 'delly', 'svaba': 'svaba', 'gridss2': 'gridss'}


def use_bundle(args, workdir):
    """Extract the newest tarball per tool from args.bundle and point the input arguments at them.

    Tools without a tarball are skipped (as if their directory was not given);
    returns the tarball names used.
    """
    used = []
    for tool, attr in BUNDLE_TOOLS.items():
        tarballs = sorted(Path(args.bundle).glob(f'*_{tool}_*.tar.gz'))
        if not tarballs:
            setattr(args, attr, None)
            continue
        with tarfile.open(tarballs[-1]) as tar:
            tar.extractall(workdir, filter='data')
        setattr(args, attr, Path(workdir) / tool)
        used.append(tarballs[-1].name)
    if args.sv_dir is None:
        raise SystemExit(f'No SVRecalibrator (svr) tarball in {args.bundle}')
    return used


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sv_dir', type=Path, nargs='?', help='SVRecalibrator batch output directory')
    parser.add_argument('--bundle', type=Path, nargs='?', const=BUNDLE_DIR,
                        help='Read inputs from the per-tool tarballs in this directory (default: turner2017/)')
    parser.add_argument('-f', '--fasta', type=Path, required=True, help='Reference FASTA (for Delly reconstruction)')
    parser.add_argument('-o', '--outdir', type=Path, required=True)
    parser.add_argument('--delly', type=Path, help='Directory of <sample>.vcf[.gz] Delly files')
    parser.add_argument('--svaba', type=Path, help='Directory of <sample>/<sample>.svaba.*.vcf SvABA runs')
    parser.add_argument('--svaba-vcf', default='unfiltered.sv', help='SvABA VCF flavour (default unfiltered.sv)')
    parser.add_argument('--gridss', type=Path, help='Directory of <sample>.vcf.gz GRIDSS2 files')
    parser.add_argument('--window', type=int, default=100, help='Maximum distance at each AA breakpoint (bp)')
    parser.add_argument('--pass-only', action='store_true', help='Keep only FILTER=PASS calls for every caller')
    parser.add_argument('--all-samples', action='store_true',
                        help='Plot all samples; by default only samples with every caller VCF present')
    args = parser.parse_args()
    if args.bundle:
        with tempfile.TemporaryDirectory() as workdir:
            inputs = use_bundle(args, workdir)
            rows, inventory, callers = compare(args)
    elif args.sv_dir:
        inputs = str(args.sv_dir)
        rows, inventory, callers = compare(args)
    else:
        parser.error('give SV_DIR (with caller directories) or --bundle')
    args.outdir.mkdir(parents=True, exist_ok=True)
    complete = {i['sample_id'] for i in inventory if all(i[f'{c}_vcf'] for c in callers)}
    in_samples = rows if args.all_samples else [r for r in rows if r['sample_id'] in complete]
    scopes = {'identified_by_all_callers': [r for r in in_samples if r['identified_by_all_callers']],
              'all_rows': in_samples}
    summary = dict(inputs=inputs,
                   window_bp=args.window, pass_only=args.pass_only, svaba_vcf=args.svaba_vcf,
                   samples_with_all_callers=len(complete), samples_total=len(inventory))
    suffix = ', PASS only' if args.pass_only else ''
    for name, subset in scopes.items():
        summary[name] = summarize(subset, callers)
        n, samples = summary[name]['rows'], summary[name]['samples']
        scope = (f'SVs identified by all callers: n={n}, {samples} samples' if name != 'all_rows'
                 else f'all SVRecalibrator rows: n={n}, {samples} samples') + suffix
        prefix = str(args.outdir / f'junctions_{name}')
        plot_hist(subset, callers, prefix, f'Junction length distributions\n({scope})')
        plot_classes(summary[name], callers, prefix, f'Junction characterization per tool\n({scope})')
    write_csv(args.outdir / 'callers_consolidated.csv', rows)
    write_csv(args.outdir / 'sample_inventory.csv', inventory)
    (args.outdir / 'summary.json').write_text(json.dumps(summary, indent=2, default=int) + '\n')
    print(json.dumps(summary, indent=2, default=int))


if __name__ == '__main__':
    main()
