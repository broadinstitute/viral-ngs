# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest, os, tempfile
import util, util.file, read_utils
from test import assert_equal_contents, TestCaseWithTmp

class TestPurgeUnmated(TestCaseWithTmp) :
    def test_purge_unmated(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        outFastq1 = util.file.mkstempfname()
        outFastq2 = util.file.mkstempfname()
        parser = read_utils.parser_purge_unmated()
        args = parser.parse_args([inFastq1, inFastq2, outFastq1, outFastq2])
        read_utils.main_purge_unmated(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, outFastq1, expected1Fastq)
        assert_equal_contents(self, outFastq2, expected2Fastq)

class TestFastqToFasta(TestCaseWithTmp) :
    def test_fastq_to_fasta(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outFasta = util.file.mkstempfname()
        parser = read_utils.parser_fastq_to_fasta()
        args = parser.parse_args([inFastq, outFasta])
        read_utils.main_fastq_to_fasta(args)

        # Check that results match expected
        expectedFasta = os.path.join(myInputDir, 'expected.fasta')
        assert_equal_contents(self, outFasta, expectedFasta)


if __name__ == '__main__':
    unittest.main()
