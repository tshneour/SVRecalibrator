# Junction comparison methods: AA, SVRecalibrator, Delly2, SvABA, GRIDSS2

Implemented in `compare_callers.py`, with VCF parsing in `caller_vcfs.py` and the
Delly2 junction reconstruction in `delly_junctions.py`.

Tools are compared on each junction's **signed length** and the **category**
derived from it. Junction sequences are recorded (`*_junction_seq`) but are not
compared base by base.

## 1. Reference set and junction values for AA and SVRecalibrator

The reference set is the SVRecalibrator (SVR) output rows (`final_augmented.tsv`),
one per AA breakpoint pair. Each row keeps the original AA coordinates
(`break_chrom1/2`, `break_pos1/2`, `break_orientation`). AA and SVR junction
values come from that same row, so they need no matching.

Every junction is described by one **signed length**: homology > 0, blunt = 0,
insertion < 0. The **category** follows from the sign.

- **AA:** `homology_len`. AA writes `None` when it found no split reads; that
  counts as *no junction call*.
- **SVR:** combines the split-read length (`sp_hom_len`) and the assembly
  (scaffold) length (`sc_hom_len`):
  - if only one is present, that value is used;
  - if both are present and equal, that value is used;
  - if both are present and differ, the row gets its own category,
    *split/scaffold disagree*;
  - if neither is present, it is *no junction call*. A split value of `N/A`
    (insertion chosen without a sequence) also counts as no call.

## 2. Loading each caller's calls

All callers were run on the same BAMs. Every FILTER value is kept unless
`--pass-only` is set. Each call becomes one breakpoint pair (chrom1, pos1,
chrom2, pos2, orientation), using AA's orientation convention: `+` means the
retained sequence lies to the left of the breakpoint. Chromosome names are
compared without the `chr` prefix.

- **GRIDSS2:** breakend records are paired using `MATEID`/`PARID`, and each pair
  is kept once. Single breakends are dropped.
  - Orientation comes from the bracket notation in the ALT field
    (`t[p[`→`+-`, `t]p]`→`++`, `]p]t`→`-+`, `[p[t`→`--`).
  - Inserted bases are the ALT bases beyond the REF base.
  - `IMPRECISE` calls have no junction.
- **SvABA:** `*.svaba.unfiltered.sv.vcf`, keeping only the `:1` record of each
  `:1`/`:2` pair.
  - Discordant-read-only calls (`EVDNC=DSCRD`) have no junction.
  - Otherwise the length is taken from the `INSERTION` or `HOMSEQ` field.
- **Delly2:** standalone insertion calls (`SVTYPE=INS`) are excluded because they
  are not breakpoint pairs.
  - Delly's own `INSLEN` is always 0 for breakpoint calls, so junctions are
    rebuilt from its `CONSENSUS` sequence (section 4).
  - `IMPRECISE` calls have no junction.

The length for GRIDSS2 and SvABA is set as follows (inserted bases take priority):

```
junction_len(inserted, homseq, homlen):
    if inserted:      return -len(inserted)
    if homseq:        return +len(homseq)
    if homlen > 0:    return +homlen
    return 0          # blunt
```

## 3. Matching SVs between tools

Each SVR row is matched separately against each caller. A caller call is a
candidate if:

- both of its ends are within **W = 100 bp** (`--window`) of the corresponding
  AA ends, and
- its orientation is the same as AA's.

The two ends may be swapped (caller pos1 ↔ AA pos2), in which case the caller's
orientation string is reversed before comparing. The 1-bp difference between
AA's coordinate convention and the 1-based VCF convention is immaterial at this
window.

If there are several candidates, the choice favours a call with a resolved
junction. This gives every caller its best chance of scoring the junction.

```
for each SVR row r (AA ends: c1,p1 ; c2,p2 ; orientation o):
    for each caller C:
        cands = []
        for each call x in C, for swap in {False, True}:
            (xc1,xp1,xc2,xp2,xo) = swap ? (x.c2,x.p2,x.c1,x.p1,reverse(x.o)) : x
            if (xc1,xc2) == (c1,c2) and xo == o
               and |xp1-p1| <= W and |xp2-p2| <= W:
                cands.append(x)          # Delly calls get reconstructed here (section 4)
        if cands is empty:
            status[r,C] = "SV not called"
        else:
            best = argmin over cands of
                   ( not resolved(x),              # resolved junction first
                     max(|xp1-p1|, |xp2-p2|),      # then nearest on the farther end
                     |xp1-p1| + |xp2-p2|,          # then smallest summed distance
                     x.id )                        # deterministic tie-break
            status[r,C] = resolved(best) ? category(best.len) : "SV called, no junction"

    identified_by_all_callers[r] = every caller has ≥1 candidate
```

Matching is not one-to-one: one caller call can match more than one SVR row.

The **main comparison** uses only rows where every caller detected the SV
(`identified_by_all_callers`). A caller counts as detecting the SV even if it
gave no junction for it.

## 4. Rebuilding Delly2 junctions

Delly aligns its consensus to the reference with a scoring scheme that penalises
a gap no more than a mismatch. As a result, inserted bases are absorbed into the
flanks as mismatches and `INSLEN` stays 0. For each PRECISE call, the split-read
`CONSENSUS` is realigned to the two reference flanks implied by the call.

```
A = flank at pos1: 200 bp on the retained side + 40 bp past the breakpoint,
    reverse-complemented if sign1 is '-'
B = flank at pos2: same construction, reverse-complemented if sign2 is '+'
    # both are oriented so a junction read runs A -> B

for s in {consensus, revcomp(consensus)}:
    a = local_align(s, A);  b = local_align(s, B)
        # match +1, mismatch -3, gap open -5, gap extend -2
    require each alignment >= 20 bp at >= 90% identity,
            and a starts and ends before b on s
    overlap = a.query_end - b.query_start
    len = overlap          # > 0: bases aligned to both flanks = homology
                           # < 0: bases aligned to neither = insertion
                           #   0: blunt
keep the strand with the higher summed score

resolved if the rebuilt breakpoints lie within 20 bp (+ homology length) of
Delly's reported positions; otherwise the call counts as "SV called, no junction"
```

The 40 bp that each flank extends past the breakpoint lets microhomology align to
both flanks. Validation: rebuilt lengths matched SVR exactly for 351 of 370
junctions (54 of 57 insertions). The unit tests (`tests/test_compare_callers.py`)
use synthetic insertion, blunt and inversion junctions.

## 5. Comparing junctions

For each caller, the comparison uses only rows where both SVR and that caller
gave a junction.

- **Exact length:** SVR's signed length equals the caller's.
- **Same category:** both are homology, both blunt, or both insertion.
- **Insertion recovery:** the fraction of SVR insertions (length < 0) that the
  caller also calls as insertions.
- **Four-way agreement:** among rows where SVR and all three callers gave a
  junction, how many share one length or one category.

The histograms clip values beyond ±50 bp into the end bins for display only.

## Limitation

Because only lengths are compared, two tools can agree on "−7 bp" while
reporting different inserted bases. Also, homology cannot always be placed at a
unique position, so the same homology can be written as different sequences at
shifted positions. A sequence check (strand-aware, and allowing a shift for
homology) would turn "same length" into "same sequence".
