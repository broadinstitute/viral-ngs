#!/usr/bin/env python
# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"

import unittest, os, sys, tempfile
import taxon_filter, util.file, tools.last


class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, main_fun, parser_fun in taxon_filter.__commands__:
            parser = parser_fun()
            helpstring = parser.format_help()


class TestTrimmomatic(unittest.TestCase) :

    def setUp(self) :
        util.file.set_tmpDir('TestTrimmomatic')

    def tearDown(self) :
        util.file.destroy_tmpDir()

    def test_trimmomatic(self) :
        commonInputDir = util.file.get_test_input_path()
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        pairedOutFastq1 = util.file.mkstempfname()
        pairedOutFastq2 = util.file.mkstempfname()
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        parser = taxon_filter.parser_trim_trimmomatic()
        args = parser.parse_args([ inFastq1, inFastq2, pairedOutFastq1,
            pairedOutFastq2, clipFasta])
        taxon_filter.main_trim_trimmomatic(args)

        # Check that results match expected
        expected1 = open(os.path.join(myInputDir, 'expected1.fastq'), 'rb')
        expected2 = open(os.path.join(myInputDir, 'expected2.fastq'), 'rb')
        result1 = open(pairedOutFastq1, 'rb')
        result2 = open(pairedOutFastq2, 'rb')

        self.assertEqual(result1.read(), expected1.read())
        self.assertEqual(result2.read(), expected2.read())

        expected1.close(); expected2.close();
        result1.close(); result2.close();


class TestFilterLastal(unittest.TestCase) :

    def setUp(self) :
        util.file.set_tmpDir('TestFilterLastal')

    def tearDown(self) :
        util.file.destroy_tmpDir()

    def test_filter_lastal(self) :
        # Create refDbs
        commonInputDir = util.file.get_test_input_path()
        myInputDir = util.file.get_test_input_path(self)

        refFasta = os.path.join(commonInputDir, 'ebola.fasta')
        dbsDir = tempfile.mkdtemp()
        refDbs = os.path.join(dbsDir, 'ebola')
        lastdbPath = tools.last.Lastdb().install_and_get_path()

        os.system('{lastdbPath} {refDbs} {refFasta}'.format(**locals()))

        # Call main_filter_lastal
        inFastq = os.path.join( myInputDir, 'in.fastq')
        outFastq = util.file.mkstempfname()
        args = taxon_filter.parser_filter_lastal().parse_args([inFastq, refDbs,
            outFastq])
        taxon_filter.main_filter_lastal(args)

        # Check that results match expected
        expectedFastq = open(os.path.join(myInputDir, 'expected.fastq'), 'rb')
        resultFastq = open('{}.fastq'.format(outFastq), 'rb')
        self.assertEqual(resultFastq.read(), expectedFastq.read())

        expectedFastq.close(); resultFastq.close()


if __name__ == '__main__':
    unittest.main()
