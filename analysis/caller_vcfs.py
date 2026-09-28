#!/usr/bin/env python3
"""Parse Delly, SvABA and GRIDSS2 VCFs into one junction record per breakpoint pair.

Each record has: caller, key, id, chrom1, pos1, chrom2, pos2, orientation (AA
convention: + = retained sequence left of the breakpoint), filter, qual,
resolved (junction sequence available), junction_len (homology > 0, blunt 0,
insertion < 0; None when unresolved), junction_seq, plus caller-specific extras.
No filtering is applied; FILTER is carried through.
"""
import gzip
import re

from consolidate_delly import read_vcf as read_delly_vcf

BRACKET = re.compile(r'([\[\]])([^\[\]:]+):(\d+)([\[\]])')


def open_vcf(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def vcf_records(path):
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            info = dict(i.split('=', 1) if '=' in i else (i, True) for i in f[7].split(';'))
            yield f, info


def parse_bnd_alt(ref, alt):
    """Return (mate chrom, mate pos, orientation, inserted bases) for a VCF BND ALT, or None.

    t[p[ -> '+-'   t]p] -> '++'   ]p]t -> '-+'   [p[t -> '--'
    Inserted bases are those beyond the REF base on the local side of the bracket.
    """
    m = BRACKET.search(alt)
    if not m or m.group(1) != m.group(4):
        return None
    before, after = alt[:m.start()], alt[m.end():]
    s1 = '+' if before else '-'
    s2 = '-' if m.group(1) == '[' else '+'
    if before:
        inserted = before[len(ref):] if before.startswith(ref) else before[1:]
    else:
        inserted = after[:-len(ref)] if after.endswith(ref) else after[:-1]
    return m.group(2), int(m.group(3)), s1 + s2, inserted.upper()


def signed_length(inserted, homseq, homlen):
    if inserted:
        return -len(inserted), inserted
    if homseq:
        return len(homseq), homseq
    if homlen and int(homlen) > 0:
        return int(homlen), ''
    return 0, ''


def read_svaba(path):
    """SvABA SV VCF (filtered or unfiltered). DSCRD (discordant-only) calls are unresolved."""
    calls = []
    for f, info in vcf_records(path):
        if not f[2].endswith(':1'):  # each pair appears as :1 and :2
            continue
        parsed = parse_bnd_alt(f[3], f[4])
        if parsed is None:
            continue
        chrom2, pos2, ori, _ = parsed
        evidence = info.get('EVDNC', '')
        resolved = evidence != 'DSCRD'
        length, seq = (signed_length(info.get('INSERTION', '').upper(), info.get('HOMSEQ', '').upper(),
                                     info.get('HOMLEN')) if resolved else (None, ''))
        calls.append(dict(caller='svaba', key=f'{path.name}:{f[2]}', id=f[2], chrom1=f[0], pos1=int(f[1]),
                          chrom2=chrom2, pos2=pos2, orientation=ori, filter=f[6], qual=f[5],
                          resolved=resolved, junction_len=length, junction_seq=seq, evidence=evidence,
                          span=info.get('SPAN', '')))
    return calls


def read_gridss(path):
    """GRIDSS2 VCF. Single breakends are skipped; IMPRECISE calls are unresolved."""
    calls, seen = [], set()
    for f, info in vcf_records(path):
        mate = info.get('MATEID') or info.get('PARID')
        if not mate or mate in seen:
            continue
        parsed = parse_bnd_alt(f[3], f[4])
        if parsed is None:
            continue
        seen.add(f[2])
        chrom2, pos2, ori, inserted = parsed
        resolved = 'IMPRECISE' not in info
        length, seq = (signed_length(inserted, info.get('HOMSEQ', '').upper(), info.get('HOMLEN'))
                       if resolved else (None, ''))
        calls.append(dict(caller='gridss', key=f'{path.name}:{f[2]}', id=f[2], chrom1=f[0], pos1=int(f[1]),
                          chrom2=chrom2, pos2=pos2, orientation=ori, filter=f[6], qual=f[5],
                          resolved=resolved, junction_len=length, junction_seq=seq,
                          evidence=f"AS={info.get('AS', '')};RAS={info.get('RAS', '')};SR={info.get('SR', '')};RP={info.get('RP', '')}"))
    return calls


def read_delly(path):
    """Delly VCF; junctions are reconstructed later from CONSENSUS (see delly_junctions.py)."""
    calls = []
    for c in read_delly_vcf(path):
        if c['delly_svtype'] == 'INS':  # standalone insertions, not breakpoint pairs
            continue
        calls.append(dict(caller='delly', key=c['delly_key'], id=c['delly_id'], chrom1=c['delly_chrom1'],
                          pos1=c['delly_pos1'], chrom2=c['delly_chrom2'], pos2=c['delly_pos2'],
                          orientation=c['delly_orientation'], filter=c['delly_filter'], qual=c['delly_qual'],
                          resolved=None, junction_len=None, junction_seq='',
                          evidence=f"{'PRECISE' if c['delly_precise'] else 'IMPRECISE'};PE={c['delly_pe']};SR={c['delly_sr']}",
                          precise=c['delly_precise'], consensus=c['delly_consensus'],
                          raw_homlen=c['delly_homlen'], raw_inslen=c['delly_inslen'], svtype=c['delly_svtype']))
    return calls
