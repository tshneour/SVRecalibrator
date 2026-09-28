#!/usr/bin/env python3
"""Reconstruct junction homology/insertion from Delly CONSENSUS sequences.

Delly v1.7.3 reports INSLEN=0 for every precise non-INS call: its consensus-to-
reference alignment scores gaps and mismatches equally, so inserted bases are
absorbed into the flanks. The split-read CONSENSUS still spans the junction, so
here it is realigned to the two reference flanks implied by the call.

Each flank is oriented so that reading the junction goes A-side -> B-side:
  A (pos1, sign +): ref[..pos1]      A (sign -): revcomp(ref[pos1..])
  B (pos2, sign -): ref[pos2..]      B (sign +): revcomp(ref[..pos2])
Flanks extend SLOP bp past each breakpoint so microhomology can align to both.
Consensus bases aligned to both flanks are homology (positive length); bases
aligned to neither are inserted sequence (negative length); 0 is blunt.
"""
from Bio.Align import PairwiseAligner

FLANK = 200
SLOP = 40
MIN_ALIGNED = 20
MIN_IDENTITY = 0.9
MAX_SHIFT = 20

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rev_comp(seq):
    return seq.translate(COMP)[::-1]


def make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2
    return aligner


def fetch(fasta, chrom, start, end):
    """0-based half-open fetch, clipped to the contig; returns (seq, clipped start)."""
    name = chrom if chrom in fasta.references else ('chr' + chrom if 'chr' + chrom in fasta.references
                                                     else chrom.removeprefix('chr'))
    start = max(0, start)
    end = min(fasta.get_reference_length(name), end)
    return fasta.fetch(name, start, end).upper(), start


def oriented_flank(fasta, chrom, pos, sign, leads_in):
    """Flank oriented along the junction read, plus a map back to 1-based ref coordinates.

    leads_in: True for the A side (flank ends at the junction), False for B.
    pos is the 1-based breakpoint base as reported in the VCF.
    """
    if sign == '+':  # retained sequence lies at and to the left of pos
        seq, start = fetch(fasta, chrom, pos - FLANK, pos + SLOP)
    else:
        seq, start = fetch(fasta, chrom, pos - 1 - SLOP, pos - 1 + FLANK)
    forward = (sign == '+') if leads_in else (sign == '-')
    if forward:
        return seq, lambda i: start + i + 1
    n = len(seq)
    return rev_comp(seq), lambda i: start + (n - 1 - i) + 1


def best_alignment(aligner, query, target):
    aln = next(iter(aligner.align(target, query)), None)
    if aln is None or aln.score <= 0:
        return None
    t0, t1 = aln.aligned[0][0][0], aln.aligned[0][-1][1]
    q0, q1 = aln.aligned[1][0][0], aln.aligned[1][-1][1]
    matches = sum(target[a] == query[b] for (ts, te), (qs, qe) in zip(*aln.aligned)
                  for a, b in zip(range(ts, te), range(qs, qe)))
    span = max(t1 - t0, q1 - q0)
    return dict(score=aln.score, t0=t0, t1=t1, q0=q0, q1=q1, identity=matches / span, length=q1 - q0)


def resolve(consensus, a_flank, a_map, b_flank, b_map, aligner):
    """Resolve one strand of the consensus against oriented flanks; None if unusable."""
    a = best_alignment(aligner, consensus, a_flank)
    b = best_alignment(aligner, consensus, b_flank)
    if a is None or b is None:
        return None
    for aln in (a, b):
        if aln['length'] < MIN_ALIGNED or aln['identity'] < MIN_IDENTITY:
            return None
    # A must precede B along the consensus, without one alignment containing the other.
    if not (a['q0'] < b['q0'] and a['q1'] < b['q1']):
        return None
    overlap = a['q1'] - b['q0']
    if overlap > 0:
        length, seq = overlap, consensus[b['q0']:a['q1']]
    else:
        length, seq = overlap, consensus[a['q1']:b['q0']]
    return dict(length=length, seq=seq, score=a['score'] + b['score'],
                a_end=a_map(a['t1'] - 1), b_start=b_map(b['t0']),
                a_len=a['length'], b_len=b['length'],
                identity=min(a['identity'], b['identity']))


def reconstruct(call, fasta, aligner=None):
    """Return dict with dl_status, dl_junction_len, dl_junction_seq, and diagnostics."""
    aligner = aligner or make_aligner()
    consensus = (call.get('consensus') or '').upper()
    ori = call.get('orientation') or ''
    if not consensus or len(ori) != 2:
        return dict(dl_status='no_consensus')
    a_flank, a_map = oriented_flank(fasta, call['chrom1'], int(call['pos1']), ori[0], True)
    b_flank, b_map = oriented_flank(fasta, call['chrom2'], int(call['pos2']), ori[1], False)
    options = [r for r in (resolve(s, a_flank, a_map, b_flank, b_map, aligner)
                           for s in (consensus, rev_comp(consensus))) if r]
    if not options:
        return dict(dl_status='unaligned')
    best = max(options, key=lambda r: r['score'])
    # The resolved junction must lie near Delly's own breakpoints (allowing for homology).
    shift = max(abs(best['a_end'] - int(call['pos1'])), abs(best['b_start'] - int(call['pos2'])))
    status = 'resolved' if shift <= MAX_SHIFT + max(best['length'], 0) else 'shifted'
    return dict(dl_status=status, dl_junction_len=best['length'], dl_junction_seq=best['seq'],
                dl_pos1=best['a_end'], dl_pos2=best['b_start'], dl_shift=shift,
                dl_flank_len1=best['a_len'], dl_flank_len2=best['b_len'],
                dl_min_identity=round(best['identity'], 3))
