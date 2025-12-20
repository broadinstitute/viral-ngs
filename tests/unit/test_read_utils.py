# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest
import argparse
import filecmp
import os
import glob

import read_utils
import shutil
import tempfile
import tools
import tools.bwa
import tools.samtools
import util
import util.file
from test import TestCaseWithTmp, assert_equal_bam_reads


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in read_utils.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


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

    def test_bwamem_idxstats_with_filtering(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outBam = util.file.mkstempfname('.bam', directory=self.tempDir)
        outStats = util.file.mkstempfname('.stats.txt', directory=self.tempDir)
        read_utils.bwamem_idxstats(inBam, self.ebolaRef, outBam, outStats, filterReadsAfterAlignment=True)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertEqual(actual_count, self.samtools.count(outBam, opts=['-F', '4']))
        self.assertLess(actual_count, 18000)

        outBamNoFiltering = util.file.mkstempfname('.bam', directory=self.tempDir)
        outStatsNoFiltering = util.file.mkstempfname('.stats.txt', directory=self.tempDir)
        read_utils.bwamem_idxstats(inBam, self.ebolaRef, outBamNoFiltering, outStatsNoFiltering, filterReadsAfterAlignment=False)
        with open(outStatsNoFiltering, 'rt') as inf:
            count_without_filtering = int(inf.readline().strip().split('\t')[2])

        # the filtered count should be less than the count without filtering
        self.assertLess(actual_count, count_without_filtering)

    def test_bwamem_idxstats_no_bam_output(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = util.file.mkstempfname('.stats.txt', directory=self.tempDir)
        read_utils.bwamem_idxstats(inBam, self.ebolaRef, None, outStats)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertGreater(actual_count, 18000)


class TestMinimap2Idxstats(TestCaseWithTmp):
    """Tests for the minimap2_idxstats command with new PAF-based implementation."""

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        self.ebolaRef = util.file.mkstempfname('.fasta', directory=self.tempDir)
        shutil.copyfile(os.path.join(util.file.get_test_input_path(), 'G5012.3.fasta'), self.ebolaRef)

    def test_minimap2_idxstats(self):
        """Test basic minimap2_idxstats with new signature (no outBam parameter)."""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = util.file.mkstempfname('.stats.txt', directory=self.tempDir)

        # New signature: minimap2_idxstats(inBam, refFasta, outStats, outReadlist=None, threads=None)
        read_utils.minimap2_idxstats(inBam, self.ebolaRef, outStats)

        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertGreater(actual_count, 18000)

    def test_minimap2_idxstats_with_readlist(self):
        """Test minimap2_idxstats with optional readlist output."""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = util.file.mkstempfname('.stats.txt', directory=self.tempDir)
        outReadlist = util.file.mkstempfname('.readlist.txt', directory=self.tempDir)

        read_utils.minimap2_idxstats(inBam, self.ebolaRef, outStats, outReadlist=outReadlist)

        # Verify stats file
        with open(outStats, 'rt') as inf:
            stats_count = int(inf.readline().strip().split('\t')[2])
        self.assertGreater(stats_count, 18000)

        # Verify readlist file
        with open(outReadlist, 'rt') as inf:
            read_ids = [line.strip() for line in inf if line.strip()]
        self.assertGreater(len(read_ids), 0)

        # All read IDs should be unique
        self.assertEqual(len(read_ids), len(set(read_ids)))

        # Readlist contains unique read pairs (~half of stats count for paired-end data)
        self.assertGreater(len(read_ids), 9000)


class TestFastqBam(TestCaseWithTmp):
    'Class for testing fastq <-> bam conversions'

    def test_fastq_bam(self):
        import pysam

        myInputDir = util.file.get_test_input_path(self)

        # Define file names
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        inHeader = os.path.join(myInputDir, 'inHeader.txt')
        expected1_7Sam = os.path.join(myInputDir, 'expected.java1_7.sam')
        outBamCmd = util.file.mkstempfname('.bam')
        outBamTxt = util.file.mkstempfname('.bam')

        # in1.fastq, in2.fastq -> out.bam; header params from command-line
        # Note: --JVMmemory is kept for backwards compatibility but ignored with samtools
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

        # Verify BAM was created with reads
        samtools = tools.samtools.SamtoolsTool()
        self.assertEqual(samtools.count(outBamCmd), 2)

        # Verify RG tags match expected values
        with pysam.AlignmentFile(outBamCmd, 'rb', check_sq=False) as bam:
            rg_list = bam.header.get('RG', [])
            self.assertEqual(len(rg_list), 1)
            rg = rg_list[0]
            self.assertEqual(rg.get('SM'), 'FreeSample')
            self.assertEqual(rg.get('LB'), 'Alexandria')
            self.assertEqual(rg.get('PL'), '9.75')
            self.assertEqual(rg.get('CN'), 'KareemAbdul-Jabbar')

        # Verify reads match expected (flags and sequences)
        assert_equal_bam_reads(self, outBamCmd, expected1_7Sam)

        # Test with header file
        parser = read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, outBamTxt, '--header', inHeader])
        args.func_main(args)

        # Verify header was replaced correctly
        with pysam.AlignmentFile(outBamTxt, 'rb', check_sq=False) as bam:
            rg_list = bam.header.get('RG', [])
            self.assertEqual(len(rg_list), 1)
            rg = rg_list[0]
            # Values should come from inHeader.txt
            self.assertEqual(rg.get('SM'), 'txtSample')
            self.assertEqual(rg.get('LB'), 'txtLib')
            self.assertEqual(rg.get('PL'), 'txtPlatform')
            self.assertEqual(rg.get('CN'), 'txtCenter')
            self.assertEqual(rg.get('DT'), '2014-11-10')

    def test_fastq_to_bam_empty_inputs(self):
        """Test that fastq_to_bam handles empty FASTQ files correctly.

        Empty FASTQ inputs should produce a valid BAM file with:
        - A proper BAM header (non-zero file size)
        - Zero reads (SamtoolsTool.isEmpty should return True)
        - Readable by samtools (no corruption)
        """
        # Create empty FASTQ files
        emptyFastq1 = util.file.mkstempfname('.fastq')
        emptyFastq2 = util.file.mkstempfname('.fastq')
        outBam = util.file.mkstempfname('.bam')

        # Create zero-byte FASTQ files
        open(emptyFastq1, 'w').close()
        open(emptyFastq2, 'w').close()

        # Convert empty FASTQs to BAM (should now succeed with defensive code)
        read_utils.fastq_to_bam(
            emptyFastq1,
            emptyFastq2,
            outBam,
            sampleName='EmptySample',
            picardOptions=['LIBRARY_NAME=EmptyLibrary']
        )

        # Verify the BAM file was created
        self.assertTrue(os.path.exists(outBam), "Output BAM file should exist")

        # Verify the BAM file is non-zero (contains header)
        bam_size = os.path.getsize(outBam)
        self.assertGreater(bam_size, 0, "BAM file should be non-zero (contains header)")

        # Verify the BAM file is empty (no reads) using SamtoolsTool
        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(samtools.isEmpty(outBam), "BAM should be empty (no reads)")

        # Verify the BAM file has zero reads via count
        read_count = samtools.count(outBam)
        self.assertEqual(read_count, 0, "BAM should contain zero reads")

        # Verify the BAM file is readable (has valid header)
        header_file = util.file.mkstempfname('.txt')
        samtools.dumpHeader(outBam, header_file)
        header_size = os.path.getsize(header_file)
        self.assertGreater(header_size, 0, "BAM header should be non-empty")


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

    def test_cdhit_canned_input(self):
        samtools = tools.samtools.SamtoolsTool()

        input_bam = os.path.join(util.file.get_test_input_path(self), 'input.bam')
        expected_bam = os.path.join(util.file.get_test_input_path(self), 'expected.bam')
        output_bam = util.file.mkstempfname("output.bam")
        read_utils.rmdup_cdhit_bam(
            input_bam,
            output_bam
        )

        self.assertEqual(samtools.count(output_bam), 1772)

    def test_cdhit_empty_input(self):
        samtools = tools.samtools.SamtoolsTool()
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        output_bam = util.file.mkstempfname("output.bam")
        read_utils.rmdup_cdhit_bam(
            empty_bam,
            output_bam
        )
        self.assertEqual(samtools.count(output_bam), 0)


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
        tools.mvicuna.MvicunaTool().rmdup(
            (inFastq1, inFastq2), (pairedOutFastq1, pairedOutFastq2),
            outUnpaired=unpairedOutFastq)

        # Compare to expected
        for filename in ['pairedOut.1.fastq', 'pairedOut.2.fastq', 'unpairedOut.fastq']:
            self.assertEqualContents(os.path.join(tempDir, filename), os.path.join(myInputDir, 'expected_' + filename))



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

    def test_minimap2(self):
        self.simple_execution('minimap2')

    def simple_execution(self, aligner):
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')
        outBamFiltered = util.file.mkstempfname('.outBamFiltered.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll, '--outBamFiltered', outBamFiltered, '--aligner', aligner])
        args.func_main(args)

    def test_empty_reads(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta,  '--aligner', 'minimap2',
            '--outBamAll', util.file.mkstempfname('.outBamAll.bam'),
            '--outBamFiltered', util.file.mkstempfname('.outBamFiltered.bam')])
        args.func_main(args)


class TestDownsampleBams(TestCaseWithTmp):
    def setUp(self):
        super(TestDownsampleBams, self).setUp()
        orig_larger_bam  = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        orig_smaller_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.mini.bam')
        orig_with_dup    = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned','input.bam')
        self.larger_bam  = util.file.mkstempfname('.larger.bam')
        self.smaller_bam = util.file.mkstempfname('.smaller.bam')
        self.with_dups   = util.file.mkstempfname('.with_dups.bam')
        shutil.copyfile(orig_larger_bam, self.larger_bam)
        shutil.copyfile(orig_smaller_bam, self.smaller_bam)
        shutil.copyfile(orig_with_dup, self.with_dups)

        self.samtools = tools.samtools.SamtoolsTool()

    def test_normalization_to_lowest_cardinality(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = self.samtools.count(self.smaller_bam)
        # target count not passed in since we are checking that the count of the smaller file is used
        read_utils.main_downsample_bams([self.larger_bam, self.smaller_bam], temp_dir, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_to_target_count(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = 4000
        read_utils.main_downsample_bams([self.larger_bam, self.smaller_bam], temp_dir, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_to_target_count_without_subdir(self):
        target_count = 4000
        read_utils.main_downsample_bams([self.larger_bam], out_path=None, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(os.path.dirname(self.larger_bam), '*downsampled-*.bam')))
        
        print(output_bams)
        self.assertGreater(len(output_bams), 0, msg="No output files matching *downsampled-*.bam found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_with_dedup_after(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = 1500
        read_utils.main_downsample_bams([self.with_dups], temp_dir, deduplicate_after=True, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertLess(self.samtools.count(out_bam), target_count, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_with_dedup_before(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = 1500
        read_utils.main_downsample_bams([self.with_dups], temp_dir, deduplicate_before=True, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_to_too_large_target_count(self):
        """ Should fail """
        temp_dir = tempfile.mkdtemp()

        target_count = 20000

        with self.assertRaises(ValueError):
            read_utils.main_downsample_bams([self.larger_bam, self.smaller_bam], temp_dir, specified_read_count=target_count, JVMmemory="1g")
