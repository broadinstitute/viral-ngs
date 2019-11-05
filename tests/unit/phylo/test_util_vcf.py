# Unit tests for util/vcf.py

__author__ = "dpark@broadinstitute.org"

import phylo.vcf
import util.file
import unittest
'''
TODO
make_intervals
sliding_windows
calc_maf
TabixReader - edge cases in a simple file
VcfReader - proper one-based genomic coordinates
VcfReader - haploid cases, diploid cases
VcfReader.getFullSequences and replaceAlleles - this needs a whole battery of tests,
    since it does so much and is subject to so many edge cases and weird scenarios and
    keeps breaking each time we refactor.
'''


class StubGenome:
    ''' This is a mock genome that should satisfy very simple functions in phylo.vcf
        like get_chrlens, GenomePosition, make_intervals, sliding_windows, etc.
        It simply contains a list of chromosome names and lengths.
    '''

    def __init__(self, chromlist):
        self.chrs = [c for c, clen in chromlist]
        self.chrlen = dict(chromlist)
        self.totlen = sum(clen for c, clen in chromlist)

    def chrlens(self):
        return [(c, self.chrlen[c]) for c in self.chrs]


class TestGenomePosition(unittest.TestCase):
    ''' Test the GenomePosition class which maps chr,pos pairs to a single gpos int and vice versa '''

    def test_fail_OOB_get_gpos(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        invalids = [('SDF', 0), ('SDF', 124), ('ASDF', -1), ('lala', 48), ('lala', 200), ('sdf', 80), ('la', 2),
                    (None, 3)]
        gmap = phylo.vcf.GenomePosition(genome)
        for c, p in invalids:
            self.assertRaises(Exception, gmap.get_gpos, c, p)

    def test_fail_OOB_get_chr_pos(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        invalids = [0, -1, genome.totlen + 1, genome.totlen * 2]
        gmap = phylo.vcf.GenomePosition(genome)
        for gpos in invalids:
            self.assertRaises(Exception, gmap.get_chr_pos, gpos)

    def test_fail_non_int_pos(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        invalids = [('SDF', 5.3), ('lala', '10'), ('ASDF', None)]
        gmap = phylo.vcf.GenomePosition(genome)
        for c, p in invalids:
            self.assertRaises(Exception, gmap.get_gpos, c, p)
            self.assertRaises(Exception, gmap.get_chr_pos, p)

    def test_spotcheck_edges(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        knowns = [('SDF', 1, 1), ('SDF', 123, 123), ('ASDF', 1, 124), ('ASDF', 256, 379), ('lala', 1, 380),
                  ('lala', 47, 426)]
        gmap = phylo.vcf.GenomePosition(genome)
        for c, p, gpos in knowns:
            self.assertEqual(gpos, gmap.get_gpos(c, p))
            self.assertEqual((c, p), gmap.get_chr_pos(gpos))

    def test_equality_1chrGenome(self):
        genome = StubGenome([('one chr', 10)])
        c = genome.chrs[0]
        gmap = phylo.vcf.GenomePosition(genome)
        for i in range(1, genome.totlen + 1):
            self.assertEqual(i, gmap.get_gpos(c, i))  # c,p -> gpos should produce p=gpos
            self.assertEqual((c, i), gmap.get_chr_pos(i))  # gpos -> c,p should produce p=gpos

    def test_equality_3chrGenome(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        gmap = phylo.vcf.GenomePosition(genome)
        # test gpos -> c,p -> same gpos
        for gpos in range(1, genome.totlen + 1):
            c, p = gmap.get_chr_pos(gpos)
            gpos2 = gmap.get_gpos(c, p)
            self.assertEqual(gpos, gpos2)
        # test c,p -> gpos -> same c,p
        for c, clen in genome.chrlen.items():
            for p in range(1, clen + 1):
                gpos = gmap.get_gpos(c, p)
                c2, p2 = gmap.get_chr_pos(gpos)
                self.assertEqual(c, c2)
                self.assertEqual(p, p2)

    def test_gpos_inbounds(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        gmap = phylo.vcf.GenomePosition(genome)
        for c, clen in genome.chrlen.items():
            for p in range(1, clen + 1):
                gpos = gmap.get_gpos(c, p)
                self.assertIsNotNone(gpos)
                self.assertIsInstance(gpos, int)
                self.assertLessEqual(1, gpos)
                self.assertLessEqual(gpos, genome.totlen)

    def test_chr_pos_inbounds(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        gmap = phylo.vcf.GenomePosition(genome)
        for gpos in range(1, genome.totlen + 1):
            c, p = gmap.get_chr_pos(gpos)
            self.assertIsNotNone(c)
            self.assertIn(c, genome.chrlen)
            self.assertIsNotNone(p)
            self.assertIsInstance(p, int)
            self.assertLessEqual(1, p)
            self.assertLessEqual(p, genome.chrlen[c])

    def test_unique_gpos(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        gmap = phylo.vcf.GenomePosition(genome)
        seen = set()
        for c, clen in genome.chrlen.items():
            for p in range(1, clen + 1):
                gpos = gmap.get_gpos(c, p)
                self.assertNotIn(gpos, seen)
                seen.add(gpos)
        self.assertEqual(len(seen), genome.totlen)

    def test_unique_chr_pos(self):
        genome = StubGenome([('SDF', 123), ('ASDF', 256), ('lala', 47)])
        gmap = phylo.vcf.GenomePosition(genome)
        seen = set()
        for gpos in range(1, genome.totlen + 1):
            c, p = gmap.get_chr_pos(gpos)
            self.assertNotIn((c, p), seen)
            seen.add((c, p))
        self.assertEqual(len(seen), genome.totlen)


class TestVcfReaderPositions(unittest.TestCase):
    ''' Test the OBO errors in the pysam-based VCFReader class (it's prone to such errors) '''

    def setUp(self):
        self.vcf_fname = util.file.get_test_path() + '/input/one_gene.vcf.gz'
        self.vcf_window = ('Pf3D7_13_v3', 1724817, 1726997)
        try:
            self.basestring = basestring
        except NameError:
            self.basestring = str

    def test_sample_names(self):
        expected = ['3D7', 'SenT001.08', 'SenT001.11', 'SenT002.07', 'SenT002.09']
        vcfdb = phylo.vcf.VcfReader(self.vcf_fname)
        self.assertEqual(expected, vcfdb.samples())
        for s in vcfdb.samples():
            self.assertIsInstance(s, self.basestring)

    def test_get_one_base(self):
        vcfdb = phylo.vcf.VcfReader(self.vcf_fname)
        genos = vcfdb.get_snp_genos(self.vcf_window[0], 1726432)
        alleles = [genos[s] for s in vcfdb.samples()]
        expected = ['T', 'T', 'T', 'G', 'T']
        self.assertEqual(expected, alleles)
        for a in alleles:
            self.assertIsInstance(a, self.basestring)

    def test_get_positions_edges(self):
        vcfdb = phylo.vcf.VcfReader(self.vcf_fname)
        out = list(vcfdb.get_positions(self.vcf_window[0]))
        self.assertEqual(out[0][1], self.vcf_window[1])
        self.assertEqual(out[-1][2], self.vcf_window[2])

    def test_get_range_edges(self):
        vcfdb = phylo.vcf.VcfReader(self.vcf_fname)
        out = [p for c, p, alleles, genos in vcfdb.get_range(self.vcf_window[0])]
        self.assertEqual(out[0], self.vcf_window[1])
        self.assertEqual(out[-1], self.vcf_window[2])
