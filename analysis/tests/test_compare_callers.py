import argparse
import gzip
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pysam

import bundle_caller_outputs
from caller_vcfs import parse_bnd_alt, read_gridss, read_svaba, signed_length
from compare_callers import svr_length, use_bundle
from delly_junctions import reconstruct, rev_comp

HEADER = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'


class ParserTests(unittest.TestCase):
    def test_bnd_alt_orientations_and_insertions(self):
        self.assertEqual(parse_bnd_alt('A', 'A[chr2:200['), ('chr2', 200, '+-', ''))
        self.assertEqual(parse_bnd_alt('A', 'AGT]chr2:200]'), ('chr2', 200, '++', 'GT'))
        self.assertEqual(parse_bnd_alt('A', ']chr2:200]CCA'), ('chr2', 200, '-+', 'CC'))
        self.assertEqual(parse_bnd_alt('A', '[chr2:200[A'), ('chr2', 200, '--', ''))
        self.assertIsNone(parse_bnd_alt('A', 'A.'))

    def test_signed_length(self):
        self.assertEqual(signed_length('GTC', '', None), (-3, 'GTC'))
        self.assertEqual(signed_length('', 'AC', '2'), (2, 'AC'))
        self.assertEqual(signed_length('', '', '0'), (0, ''))

    def test_svr_combined_length(self):
        self.assertEqual(svr_length(dict(sp_hom_len='-3', sc_hom_len='-3')), (-3, 'resolved'))
        # split 'N/A' is an insertion without a sequence, not a blunt junction
        self.assertEqual(svr_length(dict(sp_hom_len='N/A', sc_hom_len='')), (None, 'no_junction_call'))
        self.assertEqual(svr_length(dict(sp_hom_len='0', sc_hom_len='0')), (0, 'resolved'))
        self.assertEqual(svr_length(dict(sp_hom_len='', sc_hom_len='4')), (4, 'resolved'))
        self.assertEqual(svr_length(dict(sp_hom_len='2', sc_hom_len='3')), (None, 'split_scaffold_disagree'))

    def test_svaba_pairs_and_discordant(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 's.svaba.sv.vcf'
            path.write_text(HEADER +
                'chr1\t100\t7:1\tA\tA[chr1:900[\t30\tPASS\tEVDNC=ASSMB;INSERTION=GGT;SPAN=800\tGT\t0/1\n'
                'chr1\t900\t7:2\tC\t]chr1:100]C\t30\tPASS\tEVDNC=ASSMB;INSERTION=GGT;SPAN=800\tGT\t0/1\n'
                'chr1\t300\t8:1\tA\t]chr2:50]A\t9\tLOWMAPQ\tEVDNC=DSCRD;SPAN=-1\tGT\t0/1\n'
                'chr1\t300\t8:2\tA\tA[chr1:300[\t9\tLOWMAPQ\tEVDNC=DSCRD;SPAN=-1\tGT\t0/1\n')
            calls = read_svaba(path)
            self.assertEqual(len(calls), 2)
            self.assertEqual((calls[0]['orientation'], calls[0]['junction_len']), ('+-', -3))
            self.assertFalse(calls[1]['resolved'])
            self.assertIsNone(calls[1]['junction_len'])

    def test_gridss_pairs_once_and_skips_single_breakends(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 's.vcf'
            path.write_text(HEADER +
                'chr1\t100\tg1o\tA\tAT[chr1:900[\t50\tPASS\tMATEID=g1h;HOMLEN=0;SVTYPE=BND\tGT\t.\n'
                'chr1\t900\tg1h\tC\t]chr1:100]AC\t50\tPASS\tMATEID=g1o;HOMLEN=0;SVTYPE=BND\tGT\t.\n'
                'chr1\t500\tg2o\tA\tA.\t5\tLOW_QUAL\tSVTYPE=BND\tGT\t.\n'
                'chr2\t10\tg3o\tG\tG]chr3:40]\t50\tLOW_QUAL\tMATEID=g3h;HOMLEN=2;HOMSEQ=CA;IMPRECISE;SVTYPE=BND\tGT\t.\n')
            calls = read_gridss(path)
            self.assertEqual([c['id'] for c in calls], ['g1o', 'g3o'])
            self.assertEqual((calls[0]['junction_len'], calls[0]['junction_seq']), (-1, 'T'))
            self.assertIsNone(calls[1]['junction_len'])


class BundleTests(unittest.TestCase):
    def test_bundle_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'svr/S1').mkdir(parents=True)
            (tmp / 'svr/S1/final_augmented.tsv').write_text('break_chrom1\n')
            (tmp / 'delly').mkdir()
            (tmp / 'delly/S1.vcf').write_text(HEADER)
            (tmp / 'svaba/S1').mkdir(parents=True)
            (tmp / 'svaba/S1/S1.svaba.unfiltered.sv.vcf').write_text(HEADER)
            (tmp / 'gridss').mkdir()
            with gzip.open(tmp / 'gridss/S1.vcf.gz', 'wt') as handle:
                handle.write(HEADER)
            out = tmp / 'out'
            out.mkdir()
            (out / 'turner2017_svaba_2000-01-01.tar.gz').write_text('old')
            sources = dict(svr='svr', delly2='delly', svaba='svaba', gridss2='gridss')
            for tool, source in sources.items():
                argv = ['bundle', tool, str(tmp / source), '-o', str(out)]
                with mock.patch('sys.argv', argv), mock.patch('builtins.print'):
                    bundle_caller_outputs.main()
            bundles = sorted(p.name.split('_')[1] for p in out.glob('*.tar.gz'))
            self.assertEqual(bundles, sorted(sources))  # older svaba bundle replaced

            args = argparse.Namespace(bundle=out, sv_dir=None)
            used = use_bundle(args, tmp / 'extract')
            self.assertEqual(len(used), 4)
            self.assertTrue((args.sv_dir / 'S1/final_augmented.tsv').exists())
            self.assertTrue((args.delly / 'S1.vcf').exists())
            self.assertTrue((args.svaba / 'S1/S1.svaba.unfiltered.sv.vcf').exists())
            with gzip.open(args.gridss / 'S1.vcf.gz', 'rt') as handle:
                self.assertEqual(handle.read(), HEADER)

class DellyReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(1)
        cls.tmp = tempfile.TemporaryDirectory()
        cls.ref = ''.join(random.choice('ACGT') for _ in range(4000))
        path = Path(cls.tmp.name) / 'ref.fa'
        path.write_text('>chrT\n' + cls.ref + '\n')
        pysam.faidx(str(path))
        cls.fasta = pysam.FastaFile(str(path))

    @classmethod
    def tearDownClass(cls):
        cls.fasta.close()
        cls.tmp.cleanup()

    def call(self, consensus, pos1, pos2, orientation):
        return reconstruct(dict(consensus=consensus, orientation=orientation, chrom1='chrT', pos1=pos1,
                                chrom2='chrT', pos2=pos2), self.fasta)

    def test_deletion_with_insertion(self):
        # keep ref[..1000] (1-based pos1=1000), insert, resume at 1-based 3001. Inserted
        # bases matching the adjacent reference are indistinguishable from it, so the
        # insert avoids ref[1000] (C) at its start and ref[2999] (A) at its end.
        consensus = self.ref[940:1000] + 'GATTACG' + self.ref[3000:3060]
        out = self.call(consensus, 1000, 3001, '+-')
        self.assertEqual((out['dl_status'], out['dl_junction_len'], out['dl_junction_seq']),
                         ('resolved', -7, 'GATTACG'))

    def test_blunt_deletion_either_strand(self):
        consensus = self.ref[940:1000] + self.ref[3000:3060]
        for seq in (consensus, rev_comp(consensus)):
            out = self.call(seq, 1000, 3001, '+-')
            self.assertEqual(out['dl_junction_len'], self.expected_blunt_len(1000, 3001))

    def test_inversion_with_insertion(self):
        # ++ : ref[..1000] joined to revcomp(ref[..3000])
        consensus = self.ref[940:1000] + 'TCGTG' + rev_comp(self.ref[2940:3000])
        out = self.call(consensus, 1000, 3000, '++')
        self.assertEqual((out['dl_junction_len'], out['dl_junction_seq']), (-5, 'TCGTG'))

    def expected_blunt_len(self, pos1, pos2):
        """Chance identity between bases just past each breakpoint appears as homology."""
        n = 0
        while self.ref[pos1 + n] == self.ref[pos2 - 1 + n]:
            n += 1
        m = 0
        while self.ref[pos1 - 1 - m] == self.ref[pos2 - 2 - m]:
            m += 1
        return n + m


if __name__ == '__main__':
    unittest.main()
