# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest
import argparse
import filecmp
import os
import read_utils
import shutil
import tempfile
import tools
import tools.bwa
import tools.samtools
import util
import util.file
from test import TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in read_utils.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestPurgeUnmated(TestCaseWithTmp):

    def test_purge_unmated(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        outFastq1 = util.file.mkstempfname('.fastq')
        outFastq2 = util.file.mkstempfname('.fastq')
        parser = read_utils.parser_purge_unmated(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, outFastq1, outFastq2])
        args.func_main(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        self.assertEqualContents(outFastq1, expected1Fastq)
        self.assertEqualContents(outFastq2, expected2Fastq)

    # test on FASTQs with read IDs in the style of SRA fastq-dump
    def test_purge_unmated_sra(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in_sra1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in_sra2.fastq')
        outFastq1 = util.file.mkstempfname('.fastq')
        outFastq2 = util.file.mkstempfname('.fastq')
        parser = read_utils.parser_purge_unmated(argparse.ArgumentParser())
        args = parser.parse_args(['--regex', '^@(\S+).[1|2] .*', inFastq1, inFastq2, outFastq1, outFastq2])
        args.func_main(args)

        # The expected outputs are identical to the previous case.
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        self.assertEqualContents(outFastq1, expected1Fastq)
        self.assertEqualContents(outFastq2, expected2Fastq)


class TestBwamemIdxstats(TestCaseWithTmp):

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        self.ebolaRef = util.file.mkstempfname('.fasta', directory=self.tempDir)
        shutil.copyfile(os.path.join(util.file.get_test_input_path(), 'G5012.3.fasta'), self.ebolaRef)
        self.bwa = tools.bwa.Bwa()
        self.samtools = tools.samtools.SamtoolsTool()
        self.bwa.index(self.ebolaRef)

    def test_bwamem_idxstats(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outBam = util.file.mkstempfname('.bam', directory=self.tempDir)
        outStats = util.file.mkstempfname('.stats.txt', directory=self.tempDir)
        read_utils.bwamem_idxstats(inBam, self.ebolaRef, outBam, outStats)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertEqual(actual_count, self.samtools.count(outBam, opts=['-F', '4']))
        self.assertGreater(actual_count, 18000)

    def test_bwamem_idxstats_no_bam_output(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = util.file.mkstempfname('.stats.txt', directory=self.tempDir)
        read_utils.bwamem_idxstats(inBam, self.ebolaRef, None, outStats)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertGreater(actual_count, 18000)

class TestFastqToFasta(TestCaseWithTmp):

    def test_fastq_to_fasta(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outFasta = util.file.mkstempfname('.fasta')
        parser = read_utils.parser_fastq_to_fasta(argparse.ArgumentParser())
        args = parser.parse_args([inFastq, outFasta])
        args.func_main(args)

        # Check that results match expected
        expectedFasta = os.path.join(myInputDir, 'expected.fasta')
        self.assertEqualContents(outFasta, expectedFasta)


class TestFastqBam(TestCaseWithTmp):
    'Class for testing fastq <-> bam conversions'

    def test_fastq_bam(self):
        myInputDir = util.file.get_test_input_path(self)

        # Define file names
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        inHeader = os.path.join(myInputDir, 'inHeader.txt')
        expected1_7Sam = os.path.join(myInputDir, 'expected.java1_7.sam')
        expected1_8Sam = os.path.join(myInputDir, 'expected.java1_8.sam')
        expected1_8Sam_v15 = os.path.join(myInputDir, 'expected.java1_8_v1.5.sam')
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq1')
        outBamCmd = util.file.mkstempfname('.bam')
        outBamTxt = util.file.mkstempfname('.bam')
        outSam = util.file.mkstempfname('.sam')
        outFastq1 = util.file.mkstempfname('.fastq')
        outFastq2 = util.file.mkstempfname('.fastq')
        outHeader = util.file.mkstempfname('.txt')
        outHeaderFix = util.file.mkstempfname('.fix.txt')

        # in1.fastq, in2.fastq -> out.bam; header params from command-line
        parser = read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1,
                                  inFastq2,
                                  outBamCmd,
                                  '--sampleName',
                                  'FreeSample',
                                  '--JVMmemory',
                                  '1g',
                                  '--picardOptions',
                                  'LIBRARY_NAME=Alexandria',
                                  'PLATFORM=9.75',
                                  'SEQUENCING_CENTER=KareemAbdul-Jabbar',])
        args.func_main(args)

        # Note for developers: if you're fixing the tests to handle non-bugs
        # (ie our testing here is too brittle), let's just replace a lot of this
        # in the future with code that just reads the header, sorts it, and
        # tests for equality of sorted values in the RG line (and stricter
        # equality in the non-header lines). This is kind of hacky.

        # samtools view for out.sam and compare to expected
        samtools = tools.samtools.SamtoolsTool()
        samtools.view(['-h'], outBamCmd, outSam)
        # picard.sam.FastqToSam outputs header fields in different order for
        #    java version 1.8 vs 1.7/1.6, so compare both
        self.assertTrue(filecmp.cmp(outSam,
                                    expected1_7Sam,
                                    shallow=False) or filecmp.cmp(outSam,
                                                                  expected1_8Sam,
                                                                  shallow=False) or
                                                      filecmp.cmp(outSam,
                                                                  expected1_8Sam_v15,
                                                                  shallow=False))

        # in1.fastq, in2.fastq, inHeader.txt -> out.bam; header from txt
        parser = read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, outBamTxt, '--header', inHeader])
        args.func_main(args)

        # out.bam -> out1.fastq, out2.fastq, outHeader.txt; trim 1 base from 1
        parser = read_utils.parser_bam_to_fastq(argparse.ArgumentParser())
        args = parser.parse_args([outBamTxt,
                                  outFastq1,
                                  outFastq2,
                                  '--outHeader',
                                  outHeader,
                                  '--JVMmemory',
                                  '1g',
                                  '--picardOptions',
                                  'READ1_TRIM=1',])
        args.func_main(args)

        # filter out any "PG" lines from header for testing purposes
        # I don't like this... let's replace later.
        with open(outHeader, 'rt') as inf:
            with open(outHeaderFix, 'wt') as outf:
                for line in inf:
                    if not line.startswith('@PG'):
                        outf.write(line)

        # compare to out1.fastq, out2.fastq, outHeader.txt to in and expected
        self.assertEqualContents(outFastq1, expectedFastq1)  # 1 base trimmed
        self.assertEqualContents(outFastq2, inFastq2)
        self.assertEqualContents(outHeaderFix, inHeader)


class TestRmdupUnaligned(TestCaseWithTmp):
    def test_mvicuna_canned_input(self):
        samtools = tools.samtools.SamtoolsTool()

        input_bam = os.path.join(util.file.get_test_input_path(self), 'input.bam')
        expected_bam = os.path.join(util.file.get_test_input_path(self), 'expected.bam')
        output_bam = util.file.mkstempfname("output.bam")
        read_utils.rmdup_mvicuna_bam(
            input_bam,
            output_bam
        )

        self.assertEqual(samtools.count(output_bam), samtools.count(expected_bam))

    def test_mvicuna_empty_input(self):
        samtools = tools.samtools.SamtoolsTool()
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        output_bam = util.file.mkstempfname("output.bam")
        read_utils.rmdup_mvicuna_bam(
            empty_bam,
            output_bam
        )
        self.assertEqual(samtools.count(output_bam), 0)


class TestSplitReads(TestCaseWithTmp):
    'Test various options of split_reads command.'

    def test_max_reads(self):
        'Test splitting fastq using --maxReads option, with indexLen 1.'
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outPrefix = util.file.mkstempfname()

        # Split
        parser = read_utils.parser_split_reads(argparse.ArgumentParser())
        args = parser.parse_args([inFastq, outPrefix, '--maxReads', '4', '--indexLen', '1'])
        args.func_main(args)

        # Check that results match expected
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq.1')
        expectedFastq2 = os.path.join(myInputDir, 'expected.fastq.2')
        self.assertEqualContents(outPrefix + '1', expectedFastq1)
        self.assertEqualContents(outPrefix + '2', expectedFastq2)

    def test_num_chunks(self):
        'Test spliting fastq.gz using --numChunks option, with default indexLen.'
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq.gz')
        outPrefix = util.file.mkstempfname()

        # Split
        parser = read_utils.parser_split_reads(argparse.ArgumentParser())
        args = parser.parse_args([inFastq, outPrefix, '--numChunks', '3'])
        args.func_main(args)

        # Check that results match expected
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq.01')
        expectedFastq2 = os.path.join(myInputDir, 'expected.fastq.02')
        expectedFastq3 = os.path.join(myInputDir, 'expected.fastq.03')
        self.assertEqualContents(outPrefix + '01', expectedFastq1)
        self.assertEqualContents(outPrefix + '02', expectedFastq2)
        self.assertEqualContents(outPrefix + '03', expectedFastq3)

    def test_fasta(self):
        'Test splitting fasta file.'
        myInputDir = util.file.get_test_input_path(self)
        inFasta = os.path.join(myInputDir, 'in.fasta')
        outPrefix = util.file.mkstempfname()

        # Split
        parser = read_utils.parser_split_reads(argparse.ArgumentParser())
        args = parser.parse_args([inFasta, outPrefix, '--numChunks', '2', '--format', 'fasta'])
        args.func_main(args)

        # Check that results match expected
        expectedFasta1 = os.path.join(myInputDir, 'expected.fasta.01')
        expectedFasta2 = os.path.join(myInputDir, 'expected.fasta.02')
        self.assertEqualContents(outPrefix + '01', expectedFasta1)
        self.assertEqualContents(outPrefix + '02', expectedFasta2)


class TestAlignAndFix(TestCaseWithTmp):
    def setUp(self):
        super(TestAlignAndFix, self).setUp()
        orig_ref = os.path.join(util.file.get_test_input_path(), 'ebov-makona.fasta')
        self.refFasta = util.file.mkstempfname('.ref.fasta')
        shutil.copyfile(orig_ref, self.refFasta)

    def test_novoalign(self):
        self.simple_execution('novoalign')

    def test_bwa(self):
        self.simple_execution('bwa')

    def simple_execution(self, aligner):
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')
        outBamFiltered = util.file.mkstempfname('.outBamFiltered.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll, '--outBamFiltered', outBamFiltered, '--aligner', aligner])
        args.func_main(args)


class TestMvicuna(TestCaseWithTmp):
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

    def test_mvicuna(self):
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)

        # Run mvicuna
        inFastq1 = os.path.join(myInputDir, 'in.1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in.2.fastq')
        pairedOutFastq1 = os.path.join(tempDir, 'pairedOut.1.fastq')
        pairedOutFastq2 = os.path.join(tempDir, 'pairedOut.2.fastq')
        unpairedOutFastq = os.path.join(tempDir, 'unpairedOut.fastq')
        args = read_utils.parser_dup_remove_mvicuna(argparse.ArgumentParser()).parse_args(
            [inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2, '--unpairedOutFastq', unpairedOutFastq])
        args.func_main(args)

        # Compare to expected
        for filename in ['pairedOut.1.fastq', 'pairedOut.2.fastq', 'unpairedOut.fastq']:
            self.assertEqualContents(os.path.join(tempDir, filename), os.path.join(myInputDir, 'expected_' + filename))
