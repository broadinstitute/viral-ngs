# Unit tests for interhost.py

__author__ = "PLACEHOLDER"

import interhost
import test, util.file
import unittest, shutil, argparse, os

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in interhost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestCoordMapper(test.TestCaseWithTmp):
    def setUp(self):
        super(TestCoordMapper, self).setUp()
        pass
    def test_coord_mapper(self) :
        myInputDir = util.file.get_test_input_path(self)
        
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        cm = interhost.CoordMapper(os.path.join(myInputDir, 'A.fasta'),
                                   os.path.join(myInputDir, 'B.fasta'))

        expLists = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, [11, 13], 14, 15, 16, 17,
                        18, 19, 20, 21],
                    [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16,
                        17, 18, 19, 20, 21],
                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]]

        for mapper, fromChrom, goodRange, toChrom, expected in [
                [cm.mapAtoB, 'chr1', range(3, 22), 'first_chrom', expLists[0]],
                [cm.mapBtoA, 'first_chrom', range(1, 22), 'chr1', expLists[1]],
                [cm.mapAtoB, 'chr2', range(1, 14), 'second_chr', expLists[2]]] :
            result = [mapper(fromChrom, pos) for pos in goodRange]
            for chrom, mappedPos in result :
                self.assertEqual(chrom, toChrom)
            self.assertEqual(expected,
                             [mappedPos for chrom, mappedPos in result])

        # Check that IndexError is raised when past ends of other sequence
        for pos in [-1, 0, 1, 2, 22, 23, 24] :
            with self.assertRaises(IndexError) :
                cm.mapAtoB('chr1', pos)
