# Unit tests for Trimmomatic tool

__author__ = "dpark@broadinstitute.org"

import os
import unittest
import viral_ngs.core
import viral_ngs.core
import viral_ngs.core.trimmomatic
from test import TestCaseWithTmp, assert_equal_contents

class TestTrimmomatic(TestCaseWithTmp):

    def test_trimmomatic_paired(self):
        myInputDir = viral_ngs.core.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        pairedOutFastq1 = viral_ngs.core.file.mkstempfname('.out1.fastq')
        pairedOutFastq2 = viral_ngs.core.file.mkstempfname('.out2.fastq')
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        viral_ngs.core.trimmomatic.TrimmomaticTool().execute(inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2, clipFasta)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, pairedOutFastq1, expected1Fastq)
        assert_equal_contents(self, pairedOutFastq2, expected2Fastq)


    def test_trimmomatic_paired_maxinfo(self):
        myInputDir = viral_ngs.core.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        with viral_ngs.core.file.tempfnames(('.out1.fastq', '.out2.fastq')) as (pairedOutFastq1, pairedOutFastq2):
            viral_ngs.core.trimmomatic.TrimmomaticTool().execute(inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2, clipFasta,
                                                        maxinfo_target_length=30, maxinfo_strictness=.3)

            # Check that results match expected
            expected1Fastq = os.path.join(myInputDir, 'expected1.maxinfo.fastq')
            expected2Fastq = os.path.join(myInputDir, 'expected2.maxinfo.fastq')
            assert_equal_contents(self, pairedOutFastq1, expected1Fastq)
            assert_equal_contents(self, pairedOutFastq2, expected2Fastq)


    def test_trimmomatic_single(self):
        myInputDir = viral_ngs.core.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        pairedOutFastq1 = viral_ngs.core.file.mkstempfname('.out1.fastq')
        pairedOutFastq2 = viral_ngs.core.file.mkstempfname('.out2.fastq')
        unpairedOutFastq1 = viral_ngs.core.file.mkstempfname('.out3.fastq')
        unpairedOutFastq2 = viral_ngs.core.file.mkstempfname('.out4.fastq')
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        viral_ngs.core.trimmomatic.TrimmomaticTool().execute(inFastq1, None, pairedOutFastq1, pairedOutFastq2, clipFasta,
            unpairedOutFastq1=unpairedOutFastq1, unpairedOutFastq2=unpairedOutFastq2)

        # Check that results match expected
        emptyFastq = os.path.join(myInputDir, 'empty.fastq')
        expectedFastq = os.path.join(myInputDir, 'expected1.fastq')
        assert_equal_contents(self, pairedOutFastq1, emptyFastq)
        assert_equal_contents(self, pairedOutFastq2, emptyFastq)
        assert_equal_contents(self, unpairedOutFastq1, expectedFastq)
