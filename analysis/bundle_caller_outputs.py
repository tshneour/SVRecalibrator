#!/usr/bin/env python3
"""Package one tool's outputs into a dated tarball for compare_callers.py --bundle.

One tarball per tool, so each can be updated separately:

  svr      <prefix>_svr_<date>.tar.gz      svr/<S>/final_augmented.tsv
  delly2   <prefix>_delly2_<date>.tar.gz   delly2/<S>.vcf[.gz]
  svaba    <prefix>_svaba_<date>.tar.gz    svaba/<S>/<S>.svaba.{sv,unfiltered.sv}.vcf
  gridss2  <prefix>_gridss2_<date>.tar.gz  gridss2/<S>.vcf.gz (bgzipped; rebuild indexes with tabix)

Caller VCFs are complete (not subset to AA breakpoints). Each tarball also holds
<tool>/MANIFEST.md and any --include files under <tool>/commands/. Older
tarballs for the same tool in the output directory are removed, so a rerun
replaces that tool's bundle.
"""
import argparse
import datetime
import io
import re
import tarfile
from pathlib import Path

from caller_vcfs import open_vcf

SVABA_FLAVOURS = ['sv', 'unfiltered.sv']


def svr_files(source):
    return [(f'{d.name}/final_augmented.tsv', d / 'final_augmented.tsv', d.name)
            for d in sorted(p for p in source.iterdir() if p.is_dir())
            if (d / 'final_augmented.tsv').exists()]


def delly_files(source):
    return [(p.name, p, p.name.split('.vcf')[0])
            for p in sorted([*source.glob('*.vcf'), *source.glob('*.vcf.gz')])]


def svaba_files(source):
    return [(f'{d.name}/{d.name}.svaba.{f}.vcf', d / f'{d.name}.svaba.{f}.vcf', d.name)
            for d in sorted(p for p in source.iterdir() if p.is_dir() and not p.name.endswith('.tmp'))
            for f in SVABA_FLAVOURS if (d / f'{d.name}.svaba.{f}.vcf').exists()]


def gridss_files(source):
    return [(p.name, p, p.name.split('.vcf')[0]) for p in sorted(source.glob('*.vcf.gz'))]


# tool -> (file lister, VCF header line holding the version)
TOOLS = {'svr': (svr_files, None),
         'delly2': (delly_files, r'##fileDate=(.*)'),
         'svaba': (svaba_files, r'##source=(svaba\(\S+\))'),
         'gridss2': (gridss_files, r'##gridssVersion=(.*)')}


def header_value(path, pattern):
    with open_vcf(path) as handle:
        for line in handle:
            if not line.startswith('##'):
                break
            m = re.match(pattern, line)
            if m:
                return m.group(1).strip()
    return ''


def manifest(name, tool, source, files, note):
    samples = sorted({sample for *_, sample in files})
    pattern = TOOLS[tool][1]
    vcfs = [path for _, path, _ in files if '.vcf' in path.name]
    version = header_value(vcfs[0], pattern) if pattern and vcfs else ''
    lines = [f'# {name}', '',
             f'Built {datetime.date.today().isoformat()} by `analysis/bundle_caller_outputs.py` from `{source}`.',
             '']
    if version:
        lines += [f'Version (from VCF header): {version}', '']
    if note:
        lines += [note, '']
    lines += [f'{len(files)} files, {len(samples)} samples:', '', ', '.join(samples)]
    return '\n'.join(lines) + '\n'


def add_bytes(tar, arcname, data):
    info = tarfile.TarInfo(arcname)
    info.size = len(data)
    info.mtime = int(datetime.datetime.now().timestamp())
    info.mode = 0o644
    tar.addfile(info, io.BytesIO(data))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tool', choices=TOOLS)
    parser.add_argument('source', type=Path, help='Tool output directory')
    parser.add_argument('-o', '--outdir', type=Path, required=True)
    parser.add_argument('--prefix', default='turner2017')
    parser.add_argument('--include', type=Path, nargs='*', default=[], help='Extra files stored under commands/')
    parser.add_argument('--note', default='', help='Free text added to MANIFEST.md')
    args = parser.parse_args()

    stem = f'{args.prefix}_{args.tool}'
    name = f'{stem}_{datetime.date.today().isoformat()}'
    files = TOOLS[args.tool][0](args.source)
    if not files:
        raise SystemExit(f'No {args.tool} outputs found in {args.source}')
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f'{name}.tar.gz'
    tmp = out.with_name(out.name + '.partial')
    with tarfile.open(tmp, 'w:gz') as tar:
        add_bytes(tar, f'{args.tool}/MANIFEST.md',
                  manifest(name, args.tool, args.source, files, args.note).encode())
        for arcname, path, _ in files:
            tar.add(path.resolve(), arcname=f'{args.tool}/{arcname}', recursive=False)
        for extra in args.include:
            tar.add(extra.resolve(), arcname=f'{args.tool}/commands/{extra.name}', recursive=False)
    for old in args.outdir.glob(f'{stem}_*.tar.gz'):
        if old != out:
            old.unlink()
    tmp.rename(out)
    print(f'Wrote {out} ({out.stat().st_size / 1e6:.1f} MB, {len(files)} files)')


if __name__ == '__main__':
    main()
