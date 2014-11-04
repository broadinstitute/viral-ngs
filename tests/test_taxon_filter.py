# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"

import unittest, os, sys, tempfile, shutil
import taxon_filter, util.file, tools.last, tools.bmtagger
from test import assert_equal_contents


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
        args = parser.parse_args([inFastq1, inFastq2, pairedOutFastq1,
            pairedOutFastq2, clipFasta])
        taxon_filter.main_trim_trimmomatic(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, pairedOutFastq1, expected1Fastq)
        assert_equal_contents(self, pairedOutFastq2, expected2Fastq)

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

        assert not os.system(
            '{lastdbPath} {refDbs} {refFasta}'.format(**locals()))

        # Call main_filter_lastal
        inFastq = os.path.join( myInputDir, 'in.fastq')
        outFastq = util.file.mkstempfname()
        args = taxon_filter.parser_filter_lastal().parse_args([inFastq, refDbs,
            outFastq])
        taxon_filter.main_filter_lastal(args)

        # Check that results match expected
        expectedFastq = os.path.join(myInputDir, 'expected.fastq')
        assert_equal_contents(self, outFastq + '.fastq', expectedFastq)

class TestBmtagger(unittest.TestCase) :
    """
    How test data was created:
      humanChr1Subset.fa has 200 bases from human chr1
      humanChr9Subset.fa has 200 bases from human chr9
      bmtool -d humanChr1Subset.fa -o humanChr1Subset.bitmask -w 8
      bmtool -d humanChr9Subset.fa -o humanChr9Subset.bitmask -w 8
      in[12].fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
          with arbitrary quality scores.
    """
    def setUp(self) :
        util.file.set_tmpDir('TestBmtagger')
    
    def tearDown(self) :
        util.file.destroy_tmpDir()
        
    def test_bmtagger(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        srprismPath = tools.bmtagger.SrprismTool().install_and_get_path()
        for db in ['humanChr1Subset', 'humanChr9Subset'] :
            # .map file is > 100M, so recreate instead of copying
            dbfa = os.path.join(myInputDir, db + '.fa')
            dbsrprism = os.path.join(tempDir, db + '.srprism')
            assert not os.system(
                '{srprismPath} mkindex -i {dbfa} -o {dbsrprism}'.format(
                    **locals()))
            # .bitmask and .srprism.* files must be in same dir, so copy
            shutil.copy(os.path.join(myInputDir, db + '.bitmask'), tempDir)
        
        # Partition the input files
        taxon_filter.partition_bmtagger(
            os.path.join(myInputDir, 'in1.fastq'),
            os.path.join(myInputDir, 'in2.fastq'),
            [os.path.join(tempDir, 'humanChr1Subset'),
             os.path.join(tempDir, 'humanChr9Subset')],
            os.path.join(tempDir, 'outMatch'),
            os.path.join(tempDir, 'outNoMatch'))
            
        # Compare to expected
        for case in ['Match.1', 'Match.2', 'NoMatch.1', 'NoMatch.2'] :
            assert_equal_contents(self,
                os.path.join(tempDir, 'out' + case + '.fastq'),
                os.path.join(myInputDir, 'expected.' + case + '.fastq'))

if __name__ == '__main__':
    unittest.main()
