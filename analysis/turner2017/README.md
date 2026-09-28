# Turner 2017 tool outputs

Inputs for the junction comparison in `../compare_callers.py`, for the published
Turner et al. 2017 cell lines. One dated tarball per tool, so each can be updated
separately:

| Tarball | Contents |
|---|---|
| `turner2017_svr_<date>.tar.gz` | SVRecalibrator `final_augmented.tsv` (105 samples; includes AA junction fields) |
| `turner2017_delly2_<date>.tar.gz` | complete Delly2 VCFs (162 samples) |
| `turner2017_svaba_<date>.tar.gz` | complete SvABA SV VCFs, filtered and unfiltered (105 samples), and the run script |
| `turner2017_gridss2_<date>.tar.gz` | complete GRIDSS2 VCFs, bgzipped (105 samples), and the run script |

The caller VCFs are not subset to AA breakpoints. Each tarball unpacks to its own
directory (`svr/`, `delly2/`, `svaba/`, `gridss2/`) with a `MANIFEST.md` giving
the source, tool version and samples. GRIDSS2 indexes are not included
(`tabix -p vcf <S>.vcf.gz` rebuilds them).

Run the comparison from the newest tarball of each tool. The reference is the
GRCh38 FASTA from the AA data repo, used to rebuild Delly2 junctions:

```bash
cd analysis
python compare_callers.py --bundle -f $AA_DATA_REPO/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -o OUTDIR
```

After rerunning a tool, rebuild only its tarball; the new dated tarball replaces
the old one for that tool:

```bash
T=/kelp/research/projects/turner2017
python bundle_caller_outputs.py svr $T/SVRecalibrator_analysis/SVRecalibrator_outputs_nonstrict -o turner2017 --note '...'
python bundle_caller_outputs.py delly2 $T/delly2 -o turner2017 --note '...'
python bundle_caller_outputs.py svaba $T/svaba -o turner2017 --include $T/svaba/run_one.sh --note '...'
python bundle_caller_outputs.py gridss2 $T/gridss2 -o turner2017 --include $T/gridss2/run_one.sh --note '...'
```

Each tarball must stay under GitHub's 100 MB file limit (GRIDSS2 is the largest,
75 MB on 2026-09-28).
