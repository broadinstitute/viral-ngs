# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"

import unittest, os, sys, tempfile, shutil
import taxon_filter, util.file, tools.last, tools.bmtagger
from test import assert_equal_contents, TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, main_fun, parser_fun in taxon_filter.__commands__:
            parser = parser_fun()
            helpstring = parser.format_help()


class TestTrimmomatic(TestCaseWithTmp) :

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

class TestFilterLastal(TestCaseWithTmp) :

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

class TestBmtagger(TestCaseWithTmp) :
    """
    How test data was created:
      humanChr1Subset.fa has 200 bases from human chr1
      humanChr9Subset.fa has 200 bases from human chr9
      bmtool -d humanChr1Subset.fa -o humanChr1Subset.bitmask -w 8
      bmtool -d humanChr9Subset.fa -o humanChr9Subset.bitmask -w 8
      in[12].fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
          with arbitrary quality scores.
    """
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
        outMatch   = [os.path.join(tempDir,   'outMatch.{}.fastq'.format(n))
                      for n in '12']
        outNoMatch = [os.path.join(tempDir, 'outNoMatch.{}.fastq'.format(n))
                      for n in '12']
        args = taxon_filter.parser_partition_bmtagger().parse_args(
            [os.path.join(myInputDir, 'in1.fastq'),
             os.path.join(myInputDir, 'in2.fastq'),
             os.path.join(tempDir, 'humanChr1Subset'),
             os.path.join(tempDir, 'humanChr9Subset'),
             '--outMatch', outMatch[0], outMatch[1],
             '--outNoMatch', outNoMatch[0], outNoMatch[1]])
        taxon_filter.main_partition_bmtagger(args)
            
        # Compare to expected
        for case in ['Match.1', 'Match.2', 'NoMatch.1', 'NoMatch.2'] :
            assert_equal_contents(self,
                os.path.join(tempDir, 'out' + case + '.fastq'),
                os.path.join(myInputDir, 'expected.' + case + '.fastq'))

class TestMvicuna(TestCaseWithTmp) :
    """
    Input consists of 3 read pairs.
    Second read pair is identical to first.
    Third read pair has same 5' read as first, but different 3' read.
    What Mvicuna did was create paired output files in which the 2nd read
        was deleted. It created an empty unpaired output file. Although
        it initially created the postDupRm pair, it renamed them to the output
        pair.
    [IJ:]I have no idea if this is the correct behavior, but test checks that it
        doesn't change.
    """
    def test_mvicuna(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        
        # Run mvicuna
        inFastq1 = os.path.join(myInputDir, 'in.1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in.2.fastq')
        pairedOutFastq1 = os.path.join(tempDir, 'pairedOut.1.fastq')
        pairedOutFastq2 = os.path.join(tempDir, 'pairedOut.2.fastq')
        unpairedOutFastq = os.path.join(tempDir, 'unpairedOut.fastq')
        args = taxon_filter.parser_dup_remove_mvicuna().parse_args(
            [inFastq1, inFastq2,
             pairedOutFastq1, pairedOutFastq2,
             '--unpairedOutFastq', unpairedOutFastq])
        taxon_filter.main_dup_remove_mvicuna(args)
            
        # Compare to expected
        for filename in ['pairedOut.1.fastq', 'pairedOut.2.fastq',
                         'unpairedOut.fastq'] :
            assert_equal_contents(self,
                os.path.join(tempDir, filename),
                os.path.join(myInputDir, 'expected_' + filename))


if __name__ == '__main__':
    unittest.main()
