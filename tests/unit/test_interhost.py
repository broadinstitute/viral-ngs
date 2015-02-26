# Unit tests for interhost.py

__author__ = "irwin@broadinstitute.org"

import interhost
import test, util.file
import unittest, shutil, argparse, os

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in interhost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

def makeTempFasta(seqs):
    fn = util.file.mkstempfname('.fasta')
    with open(fn, 'wt') as outf:
        for line in util.file.fastaMaker(seqs):
            outf.write(line)
    return fn

class TestCoordMapper(test.TestCaseWithTmp):
    def setUp(self):
        super(TestCoordMapper, self).setUp()
        self.genomeA = makeTempFasta([
            ('chr1',        'ATGCACGTACGTATGCAAATCGG'),
            ('chr2',        'AGTCGGTTTTCAG'),
            ])
        self.genomeB = makeTempFasta([
            ('first_chrom',   'GCACGTACGTATTTGCAAATC'),
            ('second_chr',  'AGTCGGTTTCCAC'),
            ])
        self.cm = interhost.CoordMapper(self.genomeA, self.genomeB)
    
    def test_no_indels(self):
        for pos in range(1,14):
            self.assertEqual(self.cm.mapAtoB('chr2', pos), ('second_chr', pos))
            self.assertEqual(self.cm.mapBtoA('second_chr', pos), ('chr2', pos))
    
    def test_map_indels(self) :
        expLists = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, [11, 13], 14, 15, 16, 17,
                        18, 19, 20, 21],
                    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16,
                        17, 18, 19, 20, 21],
                    ]
        for mapper, fromChrom, goodRange, toChrom, expected in [
                [self.cm.mapAtoB, 'chr1', range(3, 22), 'first_chrom', expLists[0]],
                [self.cm.mapBtoA, 'first_chrom', range(1, 22), 'chr1', expLists[1]]] :
            result = [mapper(fromChrom, pos) for pos in goodRange]
            for chrom, mappedPos in result :
                self.assertEqual(chrom, toChrom)
            self.assertEqual(expected,
                             [mappedPos for chrom, mappedPos in result])
    
    def test_oob_errors(self):
        for pos in [-1, 0, 1, 2, 22, 23, 24] :
            with self.assertRaises(IndexError) :
                self.cm.mapAtoB('chr1', pos)
        for pos in [-1, 0, 14, 15] :
            with self.assertRaises(IndexError) :
                self.cm.mapBtoA('second_chr', pos)

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
        genomeA = makeTempFasta([
            ('chr1',        'ATGCACGTACGTATGCAAATCGG'),
            ('chr2',        'AGTCGGTTTTCAG'),
            ])
        genomeB = makeTempFasta([
            ('first_chrom',   'GCACGTACGTATTTGCAAATC')
            ])
        with self.assertRaises(Exception):
            cm = interhost.CoordMapper(genomeA, genomeB)

