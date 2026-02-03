# Unit tests for interhost.py

__author__ = "irwin@broadinstitute.org"

import viral_ngs.interhost
import viral_ngs.core.file
import unittest
import argparse
import itertools

from tests import TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in interhost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


def makeTempFasta(seqs):
    fn = viral_ngs.core.file.mkstempfname('.fasta')
    with open(fn, 'wt') as outf:
        for line in viral_ngs.core.file.fastaMaker(seqs):
            outf.write(line)
    return fn


class TestCoordMapper(TestCaseWithTmp):

    def setUp(self):
        super(TestCoordMapper, self).setUp()
        self.genomeA = makeTempFasta([('chr1', 'ATGCACGTACGTATGCAAATCGG'), ('chr2', 'AGTCGGTTTTCAG'),])
        self.genomeB = makeTempFasta([('first_chrom', 'GCACGTACGTATTTGCAAATC'), ('second_chr', 'AGTCGGTTTCCAC'),])
        self.cm = interhost.CoordMapper()
        self.cm.align_and_load_sequences([self.genomeA, self.genomeB])

    def test_no_indels(self):
        for pos in range(1, 14):
            self.assertEqual(self.cm.mapAtoB('chr2', pos), ('second_chr', pos))
            self.assertEqual(self.cm.mapBtoA('second_chr', pos), ('chr2', pos))

    def test_map_indels(self):
        expLists = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, [11, 13], 14, 15, 16, 17, 18, 19, 20, 21],
                    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16, 17, 18, 19, 20, 21],]
        for mapper, fromChrom, goodRange, toChrom, expected in [
            [self.cm.mapAtoB, 'chr1', range(3, 22), 'first_chrom', expLists[0]], [self.cm.mapBtoA, 'first_chrom',
                                                                                  range(1, 22), 'chr1', expLists[1]]
        ]:
            result = [mapper(fromChrom, pos) for pos in goodRange]
            for chrom, mappedPos in result:
                self.assertEqual(chrom, toChrom)
            self.assertEqual(expected, [mappedPos for chrom, mappedPos in result])

    def test_side_param(self):
        self.assertEqual(self.cm.mapAtoB('chr1', 13), ('first_chrom', [11, 13]))
        self.assertEqual(self.cm.mapAtoB('chr1', 13, 0), ('first_chrom', [11, 13]))
        self.assertEqual(self.cm.mapAtoB('chr1', 13, -1), ('first_chrom', 11))
        self.assertEqual(self.cm.mapAtoB('chr1', 13, 1), ('first_chrom', 13))
        self.assertEqual(self.cm.mapAtoB('chr1', 12), ('first_chrom', 10))
        self.assertEqual(self.cm.mapAtoB('chr1', 12, 0), ('first_chrom', 10))
        self.assertEqual(self.cm.mapAtoB('chr1', 12, -1), ('first_chrom', 10))
        self.assertEqual(self.cm.mapAtoB('chr1', 12, 1), ('first_chrom', 10))

    def test_oob_errors(self):
        for pos in [-1, 0, 1, 2, 22, 23, 24]:
            self.assertEqual(self.cm.mapAtoB('chr1', pos), ('first_chrom', None))
        for pos in [-1, 0, 14, 15]:
            self.assertEqual(self.cm.mapBtoA('second_chr', pos), ('chr2', None))

    def test_invalid_pos_error(self):
        with self.assertRaises(TypeError):
            self.cm.mapAtoB('chr1', 1.5)
        with self.assertRaises(TypeError):
            self.cm.mapBtoA('second_chr', 4.5)

    def test_invalid_chr_error(self):
        with self.assertRaises(KeyError):
            self.cm.mapAtoB('nonexistentchr', 2)
        with self.assertRaises(KeyError):
            self.cm.mapBtoA('nonexistentchr', 2)

    def test_unequal_genomes_error(self):
        genomeA = makeTempFasta([('chr1', 'ATGCACGTACGTATGCAAATCGG'), ('chr2', 'AGTCGGTTTTCAG'),])
        genomeB = makeTempFasta([('first_chrom', 'GCACGTACGTATTTGCAAATC')])
        with self.assertRaises(viral_ngs.core.file.TranspositionError):
            cm = interhost.CoordMapper()
            cm.align_and_load_sequences([genomeA, genomeB])

    def test_map_chr_only(self):
        self.assertEqual(self.cm.mapAtoB('chr1'), 'first_chrom')
        self.assertEqual(self.cm.mapBtoA('first_chrom'), 'chr1')
        self.assertEqual(self.cm.mapAtoB('chr2'), 'second_chr')
        self.assertEqual(self.cm.mapBtoA('second_chr'), 'chr2')
        with self.assertRaises(KeyError):
            self.cm.mapAtoB('nonexistentchr')


class TestCoordMapperMultipleSeqs(TestCaseWithTmp):

    def setUp(self):
        super(TestCoordMapperMultipleSeqs, self).setUp()
        self.genomeA = makeTempFasta([
            ('chr1', 'ATGCACGTACGTATGCAAATCGG'),
            ('chr2', 'AGTCGGTTTTCAG'),
            ('chr3', 'GACTTTTGGCTGA'),
        ])
        self.genomeB = makeTempFasta([
            ('first_chrom', 'GCACGTACGTATTTGCAAATC'),
            ('second_chr', 'AGTCGGTTTCCAC'),
            ('third_chr', 'CACCTTTGGCTGA'),
        ])
        self.cm = interhost.CoordMapper()
        self.cm.align_and_load_sequences([self.genomeA, self.genomeB])

    def test_legacy_call(self):
        '''
            If mapAtoB or mapBtoA is called on a CoordMapper object with >2 sequences,
            an AssertionError should be raised.
        '''
        self.assertRaises(LookupError, self.cm.mapAtoB, 'chr1', 1)
        self.assertRaises(LookupError, self.cm.mapBtoA, 'chr2', 1)

    def test_no_indels(self):
        for pos in range(1, 14):
            self.assertEqual(self.cm.mapChr('chr2', 'second_chr', pos), ('second_chr', pos))
            self.assertEqual(self.cm.mapChr('second_chr', 'chr2', pos), ('chr2', pos))

    def test_map_indels(self):
        expLists = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, [11, 13], 14, 15, 16, 17, 18, 19, 20, 21],
                    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16, 17, 18, 19, 20, 21],]
        for mapper, fromChrom, goodRange, toChrom, expected in [
            [self.cm.mapChr, 'chr1', range(3, 22), 'first_chrom', expLists[0]], [self.cm.mapChr, 'first_chrom',
                                                                                 range(1, 22), 'chr1', expLists[1]]
        ]:
            result = [mapper(fromChrom, toChrom, pos) for pos in goodRange]
            for chrom, mappedPos in result:
                self.assertEqual(chrom, toChrom)
            self.assertEqual(expected, [mappedPos for chrom, mappedPos in result])

    def test_side_param(self):
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 13), ('first_chrom', [11, 13]))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 13, 0), ('first_chrom', [11, 13]))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 13, -1), ('first_chrom', 11))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 13, 1), ('first_chrom', 13))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 12), ('first_chrom', 10))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 12, 0), ('first_chrom', 10))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 12, -1), ('first_chrom', 10))
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', 12, 1), ('first_chrom', 10))

    def test_oob_errors(self):
        for pos in [-1, 0, 1, 2, 22, 23, 24]:
            self.assertEqual(self.cm.mapChr('chr1', 'first_chrom', pos), ('first_chrom', None))
        for pos in [-1, 0, 14, 15]:
            self.assertEqual(self.cm.mapChr('second_chr', 'chr2', pos), ('chr2', None))

    def test_invalid_pos_error(self):
        with self.assertRaises(TypeError):
            self.cm.mapChr('chr1', 'first_chrom', 1.5)
        with self.assertRaises(TypeError):
            self.cm.mapChr('second_chr', 'chr2', 4.5)

    def test_invalid_chr_error(self):
        self.assertRaises(KeyError, self.cm.mapChr, 'nonexistentchr', 'chr1', 2)
        self.assertRaises(KeyError, self.cm.mapChr, 'chr1', 'nonexistentchr', 2)

    def test_unequal_genomes_error(self):
        genomeA = makeTempFasta([('chr1', 'ATGCACGTACGTATGCAAATCGG'), ('chr2', 'AGTCGGTTTTCAG'),])
        genomeB = makeTempFasta([('first_chrom', 'GCACGTACGTATTTGCAAATC')])
        with self.assertRaises(Exception):
            cm = interhost.CoordMapper()
            cm.align_and_load_sequences([genomeA, genomeB])

    def test_duplicate_chr_names_error(self):
        genomeA = makeTempFasta([('chr1', 'ATGCACGTACGTATGCAAATCGG'),])
        genomeB = makeTempFasta([('chr1', 'GCACGTACGTATTTGCAAATC')])
        with self.assertRaises(Exception):
            cm = interhost.CoordMapper()
            cm.align_and_load_sequences([genomeA, genomeB])

    def test_multiple_input_genomes(self):
        genomeA = makeTempFasta([('chr1', 'ATGCACGTACGTATGCAAATCGG'),])
        genomeB = makeTempFasta([('first_chr', 'ATGCACTACGTATGCAAATCGG')])
        genomeC = makeTempFasta([('chr_one', 'ATGCACGTACGTATGCAATCGG')])
        cm = interhost.CoordMapper()
        cm.align_and_load_sequences([genomeA, genomeB, genomeC])
        # check that toChrom is in the map
        self.assertEqual(cm.mapChr('chr1', 'chr_one'), 'chr_one')

    def test_single_chr_error(self):
        genomeA = makeTempFasta([('chr1', 'ATGCACGTACGTATGCAAATCGG'),])
        genomeB = makeTempFasta([])
        with self.assertRaises(Exception):
            cm = interhost.CoordMapper()
            cm.align_and_load_sequences([genomeA, genomeB])

    def test_map_chr_only(self):
        self.assertEqual(self.cm.mapChr('chr1', 'first_chrom'), 'first_chrom')
        self.assertEqual(self.cm.mapChr('first_chrom', 'chr1'), 'chr1')
        self.assertEqual(self.cm.mapChr('chr2', 'second_chr'), 'second_chr')
        self.assertEqual(self.cm.mapChr('second_chr', 'chr2'), 'chr2')
        self.assertEqual(self.cm.mapChr('chr3', 'third_chr'), 'third_chr')
        self.assertEqual(self.cm.mapChr('third_chr', 'chr3'), 'chr3')
        self.assertRaises(KeyError, self.cm.mapChr, 'nonexistentchr', 'chr1')


class TestSpecificAlignments(TestCaseWithTmp):
    """ For the most part, CoordMapper2Seqs is tested implicitly when
        CoordMapper is tested. Focus here on special cases that are hard
        or impossible to get out of the aligner.
    """

    def test_basic_alignment(self):
        alignment = makeTempFasta([('s1', 'ATCG'), ('s2', 'ACCG'), ('s3', 'AG-T'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment])

    def test_unequal_len(self):
        alignment = makeTempFasta([('s1', 'AA'), ('s2', 'A'),])
        cm = interhost.CoordMapper()
        with self.assertRaises(Exception):
            cm.load_alignments([alignment])

    def test_no_real_bases_in_sample(self):
        alignment1 = makeTempFasta([('s1', 'AA'), ('s2', '--'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment1])
        with self.assertRaises(Exception):
            cm.mapChr('s1', 's2', 1)

        alignment2 = makeTempFasta([('s1', '--'), ('s2', 'AA'), ('s3', 'TT'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment2])
        with self.assertRaises(Exception):
            cm.mapChr('s2', 's1', 1)

    def test_no_real_bases_at_position(self):
        alignment = makeTempFasta([('s1', 'AT-G'), ('s2', 'AC-G'), ('s3', 'AG-T'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment])
        for i in (1, 2, 3):
            self.assertEqual(cm.mapChr('s1', 's2', i), ('s2', i))
            self.assertEqual(cm.mapChr('s2', 's1', i), ('s1', i))
            self.assertEqual(cm.mapChr('s1', 's3', i), ('s3', i))
            self.assertEqual(cm.mapChr('s3', 's1', i), ('s1', i))
            self.assertEqual(cm.mapChr('s2', 's3', i), ('s3', i))
            self.assertEqual(cm.mapChr('s3', 's2', i), ('s2', i))

    def test_aligned_gaps(self):
        alignment = makeTempFasta([('s1', 'ATCG'), ('s2', 'AC-G'), ('s3', 'AG-T'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment])
        for i in (1, 2, 3):
            self.assertEqual(cm.mapChr('s2', 's3', i), ('s3', i))
            self.assertEqual(cm.mapChr('s3', 's2', i), ('s2', i))
        for x, y in ((1, 1), (2, 2), (3, 2), (4, 3)):
            self.assertEqual(cm.mapChr('s1', 's2', x), ('s2', y))
            self.assertEqual(cm.mapChr('s1', 's3', x), ('s3', y))
        for x, y in ((1, 1), (2, [2, 3]), (3, 4)):
            self.assertEqual(cm.mapChr('s2', 's1', x), ('s1', y))
            self.assertEqual(cm.mapChr('s3', 's1', x), ('s1', y))

    def test_adjacent_gaps(self):
        alignment = makeTempFasta([
            ('s1', 'ATCTG'),
            ('s2', 'AC--G'),
            ('s3', 'A-TTG'),
            ('s4', 'A-C-G'),
            ('s5', 'A--CG'),
        ])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment])
        for x, y in ((1, 1), (2, 2), (3, 2), (4, 2), (5, 3)):
            self.assertEqual(cm.mapChr('s1', 's2', x), ('s2', y))
        for x, y in ((1, 1), (2, [2, 4]), (3, 5)):
            self.assertEqual(cm.mapChr('s2', 's1', x), ('s1', y))
        for x, y in ((1, 1), (2, 1), (3, 2), (4, 3), (5, 4)):
            self.assertEqual(cm.mapChr('s1', 's3', x), ('s3', y))
        for x, y in ((1, [1, 2]), (2, 3), (3, 4), (4, 5)):
            self.assertEqual(cm.mapChr('s3', 's1', x), ('s1', y))
        for x, y in ((1, 1), (2, [2, 3]), (3, 4)):
            self.assertEqual(cm.mapChr('s2', 's3', x), ('s3', y))
        for x, y in ((1, 1), (2, 2), (3, 2), (4, 3)):
            self.assertEqual(cm.mapChr('s3', 's2', x), ('s2', y))
        for a, b in itertools.combinations(('s2', 's4', 's5'), 2):
            for i in (1, 2, 3):
                self.assertEqual(cm.mapChr(a, b, i), (b, i))
                self.assertEqual(cm.mapChr(b, a, i), (a, i))

    def test_one_real_base(self):
        alignment = makeTempFasta([('s1', 'AC-'), ('s2', '-CA'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment])
        self.assertEqual(cm.mapChr('s1', 's2', 1), ('s2', None))
        self.assertEqual(cm.mapChr('s1', 's2', 2), ('s2', 1))
        self.assertEqual(cm.mapChr('s2', 's1', 1), ('s1', 2))
        self.assertEqual(cm.mapChr('s2', 's1', 2), ('s1', None))

    def test_exactly_two_pairs(self):
        alignment = makeTempFasta([('s1', 'A--T'), ('s2', 'AGGT'),])
        cm = interhost.CoordMapper()
        cm.load_alignments([alignment])
        self.assertEqual(cm.mapChr('s1', 's2', 1), ('s2', [1, 3]))
        self.assertEqual(cm.mapChr('s1', 's2', 2), ('s2', 4))
        self.assertEqual(cm.mapChr('s2', 's1', 1), ('s1', 1))
        self.assertEqual(cm.mapChr('s2', 's1', 2), ('s1', 1))
        self.assertEqual(cm.mapChr('s2', 's1', 3), ('s1', 1))
        self.assertEqual(cm.mapChr('s2', 's1', 4), ('s1', 2))
