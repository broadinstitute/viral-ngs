# Unit tests for Trimmomatic tool

__author__ = "dpark@broadinstitute.org"

import os
import unittest
import util.file
import tools
import tools.trimmomatic
from test import TestCaseWithTmp, assert_equal_contents

class TestTrimmomatic(TestCaseWithTmp):

    def test_trimmomatic_paired(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        pairedOutFastq1 = util.file.mkstempfname('.out1.fastq')
        pairedOutFastq2 = util.file.mkstempfname('.out2.fastq')
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        tools.trimmomatic.TrimmomaticTool().execute(inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2, clipFasta)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, pairedOutFastq1, expected1Fastq)
        assert_equal_contents(self, pairedOutFastq2, expected2Fastq)

    def test_trimmomatic_single(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        pairedOutFastq1 = util.file.mkstempfname('.out1.fastq')
        pairedOutFastq2 = util.file.mkstempfname('.out2.fastq')
        unpairedOutFastq1 = util.file.mkstempfname('.out3.fastq')
        unpairedOutFastq2 = util.file.mkstempfname('.out4.fastq')
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        tools.trimmomatic.TrimmomaticTool().execute(inFastq1, None, pairedOutFastq1, pairedOutFastq2, clipFasta,
            unpairedOutFastq1=unpairedOutFastq1, unpairedOutFastq2=unpairedOutFastq2)

        # Check that results match expected
        emptyFastq = os.path.join(myInputDir, 'empty.fastq')
        expectedFastq = os.path.join(myInputDir, 'expected1.fastq')
        assert_equal_contents(self, pairedOutFastq1, emptyFastq)
        assert_equal_contents(self, pairedOutFastq2, emptyFastq)
        assert_equal_contents(self, unpairedOutFastq1, expectedFastq)
