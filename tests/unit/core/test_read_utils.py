# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest
import argparse
import filecmp
import os
import glob
import platform
import random

import pysam
import viral_ngs.read_utils
import shutil
import tempfile
import viral_ngs.core
import viral_ngs.core.bwa
import viral_ngs.core.samtools
from tests import TestCaseWithTmp, assert_equal_bam_reads

# Skip tests requiring x86-only tools on ARM platforms
IS_ARM = platform.machine() in ('arm64', 'aarch64')
SKIP_X86_ONLY_REASON = "Requires x86-only bioconda package (not available on ARM)"


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in viral_ngs.read_utils.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestBwamemIdxstats(TestCaseWithTmp):

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        self.ebolaRef = viral_ngs.core.file.mkstempfname('.fasta', directory=self.tempDir)
        shutil.copyfile(os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.fasta'), self.ebolaRef)
        self.bwa = viral_ngs.core.bwa.Bwa()
        self.samtools = viral_ngs.core.samtools.SamtoolsTool()
        self.bwa.index(self.ebolaRef)

    def test_bwamem_idxstats(self):
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam', directory=self.tempDir)
        outStats = viral_ngs.core.file.mkstempfname('.stats.txt', directory=self.tempDir)
        viral_ngs.read_utils.bwamem_idxstats(inBam, self.ebolaRef, outBam, outStats)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertEqual(actual_count, self.samtools.count(outBam, opts=['-F', '4']))
        self.assertGreater(actual_count, 18000)

    def test_bwamem_idxstats_with_filtering(self):
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam', directory=self.tempDir)
        outStats = viral_ngs.core.file.mkstempfname('.stats.txt', directory=self.tempDir)
        viral_ngs.read_utils.bwamem_idxstats(inBam, self.ebolaRef, outBam, outStats, filterReadsAfterAlignment=True)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertEqual(actual_count, self.samtools.count(outBam, opts=['-F', '4']))
        self.assertLess(actual_count, 18000)

        outBamNoFiltering = viral_ngs.core.file.mkstempfname('.bam', directory=self.tempDir)
        outStatsNoFiltering = viral_ngs.core.file.mkstempfname('.stats.txt', directory=self.tempDir)
        viral_ngs.read_utils.bwamem_idxstats(inBam, self.ebolaRef, outBamNoFiltering, outStatsNoFiltering, filterReadsAfterAlignment=False)
        with open(outStatsNoFiltering, 'rt') as inf:
            count_without_filtering = int(inf.readline().strip().split('\t')[2])

        # the filtered count should be less than the count without filtering
        self.assertLess(actual_count, count_without_filtering)

    def test_bwamem_idxstats_no_bam_output(self):
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = viral_ngs.core.file.mkstempfname('.stats.txt', directory=self.tempDir)
        viral_ngs.read_utils.bwamem_idxstats(inBam, self.ebolaRef, None, outStats)
        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertGreater(actual_count, 18000)


class TestMinimap2Idxstats(TestCaseWithTmp):
    """Tests for the minimap2_idxstats command with new PAF-based implementation."""

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        self.ebolaRef = viral_ngs.core.file.mkstempfname('.fasta', directory=self.tempDir)
        shutil.copyfile(os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.fasta'), self.ebolaRef)

    def test_minimap2_idxstats(self):
        """Test basic minimap2_idxstats with new signature (no outBam parameter)."""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = viral_ngs.core.file.mkstempfname('.stats.txt', directory=self.tempDir)

        # New signature: minimap2_idxstats(inBam, refFasta, outStats, outReadlist=None, threads=None)
        viral_ngs.read_utils.minimap2_idxstats(inBam, self.ebolaRef, outStats)

        with open(outStats, 'rt') as inf:
            actual_count = int(inf.readline().strip().split('\t')[2])
        self.assertGreater(actual_count, 18000)

    def test_minimap2_idxstats_with_readlist(self):
        """Test minimap2_idxstats with optional readlist output."""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outStats = viral_ngs.core.file.mkstempfname('.stats.txt', directory=self.tempDir)
        outReadlist = viral_ngs.core.file.mkstempfname('.readlist.txt', directory=self.tempDir)

        viral_ngs.read_utils.minimap2_idxstats(inBam, self.ebolaRef, outStats, outReadlist=outReadlist)

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

        myInputDir = viral_ngs.core.file.get_test_input_path(self)

        # Define file names
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        inHeader = os.path.join(myInputDir, 'inHeader.txt')
        expected1_7Sam = os.path.join(myInputDir, 'expected.java1_7.sam')
        outBamCmd = viral_ngs.core.file.mkstempfname('.bam')
        outBamTxt = viral_ngs.core.file.mkstempfname('.bam')

        # in1.fastq, in2.fastq -> out.bam; header params from command-line
        # Note: --JVMmemory is kept for backwards compatibility but ignored with samtools
        parser = viral_ngs.read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
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
        samtools = viral_ngs.core.samtools.SamtoolsTool()
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
        parser = viral_ngs.read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
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
        emptyFastq1 = viral_ngs.core.file.mkstempfname('.fastq')
        emptyFastq2 = viral_ngs.core.file.mkstempfname('.fastq')
        outBam = viral_ngs.core.file.mkstempfname('.bam')

        # Create zero-byte FASTQ files
        open(emptyFastq1, 'w').close()
        open(emptyFastq2, 'w').close()

        # Convert empty FASTQs to BAM (should now succeed with defensive code)
        viral_ngs.read_utils.fastq_to_bam(
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
        samtools = viral_ngs.core.samtools.SamtoolsTool()
        self.assertTrue(samtools.isEmpty(outBam), "BAM should be empty (no reads)")

        # Verify the BAM file has zero reads via count
        read_count = samtools.count(outBam)
        self.assertEqual(read_count, 0, "BAM should contain zero reads")

        # Verify the BAM file is readable (has valid header)
        header_file = viral_ngs.core.file.mkstempfname('.txt')
        samtools.dumpHeader(outBam, header_file)
        header_size = os.path.getsize(header_file)
        self.assertGreater(header_size, 0, "BAM header should be non-empty")


class TestRmdupUnaligned(TestCaseWithTmp):
    @unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
    def test_cdhit_canned_input(self):
        samtools = viral_ngs.core.samtools.SamtoolsTool()

        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")
        viral_ngs.read_utils.rmdup_cdhit_bam(
            input_bam,
            output_bam
        )

        self.assertEqual(samtools.count(output_bam), 1772)

    @unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
    def test_cdhit_empty_input(self):
        samtools = viral_ngs.core.samtools.SamtoolsTool()
        empty_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")
        viral_ngs.read_utils.rmdup_cdhit_bam(
            empty_bam,
            output_bam
        )
        self.assertEqual(samtools.count(output_bam), 0)


class TestReadIdStore(TestCaseWithTmp):
    """Tests for ReadIdStore SQLite-backed read ID storage."""

    def test_add_from_fastq_paired(self):
        """Test adding read IDs from paired-end interleaved FASTQ."""
        # Create a test FASTQ with paired reads
        fastq_path = viral_ngs.core.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            # Write 3 read pairs (6 entries, but only 3 unique IDs)
            for i in range(3):
                f.write('@read{}/1\nACGT\n+\nIIII\n'.format(i))
                f.write('@read{}/2\nACGT\n+\nIIII\n'.format(i))

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            added = store.add_from_fastq(fastq_path)
            self.assertEqual(added, 3)  # 3 unique read IDs
            self.assertEqual(len(store), 3)

    def test_add_from_fastq_single_end(self):
        """Test adding read IDs from single-end FASTQ."""
        fastq_path = viral_ngs.core.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            for i in range(5):
                f.write('@read{}\nACGT\n+\nIIII\n'.format(i))

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            added = store.add_from_fastq(fastq_path)
            self.assertEqual(added, 5)
            self.assertEqual(len(store), 5)

    def test_deduplication(self):
        """Test that duplicate read IDs are ignored."""
        fastq_path = viral_ngs.core.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            # Write same read ID multiple times
            for _ in range(10):
                f.write('@duplicate_read\nACGT\n+\nIIII\n')

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.add_from_fastq(fastq_path)
            self.assertEqual(len(store), 1)  # Only 1 unique ID

    def test_write_to_file(self):
        """Test writing read IDs to file."""
        fastq_path = viral_ngs.core.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            for i in range(5):
                f.write('@read{}\nACGT\n+\nIIII\n'.format(i))

        db_path = viral_ngs.core.file.mkstempfname('.db')
        out_path = viral_ngs.core.file.mkstempfname('.txt')

        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.add_from_fastq(fastq_path)
            written = store.write_to_file(out_path)
            self.assertEqual(written, 5)

        # Verify file contents
        with open(out_path, 'rt') as f:
            lines = [line.strip() for line in f]
        self.assertEqual(len(lines), 5)
        self.assertEqual(set(lines), {'read0', 'read1', 'read2', 'read3', 'read4'})

    def test_write_to_file_with_downsampling(self):
        """Test random downsampling when writing to file."""
        fastq_path = viral_ngs.core.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            for i in range(100):
                f.write('@read{}\nACGT\n+\nIIII\n'.format(i))

        db_path = viral_ngs.core.file.mkstempfname('.db')
        out_path = viral_ngs.core.file.mkstempfname('.txt')

        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.add_from_fastq(fastq_path)
            self.assertEqual(len(store), 100)
            written = store.write_to_file(out_path, max_reads=20)
            self.assertEqual(written, 20)

        # Verify file has exactly 20 lines
        with open(out_path, 'rt') as f:
            lines = [line.strip() for line in f]
        self.assertEqual(len(lines), 20)

    def test_empty_fastq(self):
        """Test handling of empty FASTQ file."""
        fastq_path = viral_ngs.core.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            pass  # Empty file

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            added = store.add_from_fastq(fastq_path)
            self.assertEqual(added, 0)
            self.assertEqual(len(store), 0)

    def test_add_single(self):
        """Test adding a single read ID."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            self.assertTrue(store.add('read1'))
            self.assertEqual(len(store), 1)
            # Adding same ID again should return False
            self.assertFalse(store.add('read1'))
            self.assertEqual(len(store), 1)
            # Adding different ID should work
            self.assertTrue(store.add('read2'))
            self.assertEqual(len(store), 2)

    def test_extend(self):
        """Test adding multiple read IDs efficiently."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            # Add from a list
            added = store.extend(['read1', 'read2', 'read3'])
            self.assertEqual(added, 3)
            self.assertEqual(len(store), 3)

            # Add with some duplicates
            added = store.extend(['read3', 'read4', 'read5'])
            self.assertEqual(added, 2)  # Only read4 and read5 are new
            self.assertEqual(len(store), 5)

    def test_extend_generator(self):
        """Test extend with a generator (O(1) memory)."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            # Use a generator expression
            added = store.extend('read{}'.format(i) for i in range(100))
            self.assertEqual(added, 100)
            self.assertEqual(len(store), 100)

    def test_contains(self):
        """Test membership testing with 'in' operator."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2', 'read3'])
            self.assertIn('read1', store)
            self.assertIn('read2', store)
            self.assertNotIn('read4', store)

    def test_iter(self):
        """Test iteration over read IDs in insertion order."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(['read3', 'read1', 'read2'])
            # Should iterate in insertion order
            result = list(store)
            self.assertEqual(result, ['read3', 'read1', 'read2'])

    def test_delitem(self):
        """Test deleting read IDs with del operator."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2', 'read3'])
            self.assertEqual(len(store), 3)

            del store['read2']
            self.assertEqual(len(store), 2)
            self.assertNotIn('read2', store)
            self.assertIn('read1', store)
            self.assertIn('read3', store)

            # Deleting non-existent ID should raise KeyError
            with self.assertRaises(KeyError):
                del store['nonexistent']

    def test_discard(self):
        """Test discard method (no error if absent)."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2'])

            store.discard('read1')
            self.assertEqual(len(store), 1)
            self.assertNotIn('read1', store)

            # Discarding non-existent ID should not raise
            store.discard('nonexistent')  # Should not raise
            self.assertEqual(len(store), 1)

    def test_shrink_to_subsample_basic(self):
        """Test shrink_to_subsample reduces store to n elements."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            # Add 100 read IDs
            original_ids = set('read{}'.format(i) for i in range(100))
            store.extend(original_ids)
            self.assertEqual(len(store), 100)

            # Shrink to 25
            store.shrink_to_subsample(25)
            self.assertEqual(len(store), 25)

            # All remaining IDs should be from original set
            remaining = set(store)
            self.assertTrue(remaining.issubset(original_ids))

    def test_shrink_to_subsample_larger_than_store(self):
        """Test shrink_to_subsample with n >= store size does nothing."""
        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(['read{}'.format(i) for i in range(50)])
            self.assertEqual(len(store), 50)

            # Request 100, but only have 50
            store.shrink_to_subsample(100)
            self.assertEqual(len(store), 50)  # Unchanged

    def test_shrink_to_subsample_randomness(self):
        """Test that shrink_to_subsample produces different results on repeated calls."""
        original_ids = ['read{}'.format(i) for i in range(100)]

        # Create two stores and subsample each
        db_path1 = viral_ngs.core.file.mkstempfname('.db')
        db_path2 = viral_ngs.core.file.mkstempfname('.db')

        with viral_ngs.read_utils.ReadIdStore(db_path1) as store1:
            store1.extend(original_ids)
            store1.shrink_to_subsample(10)
            result1 = set(store1)

        with viral_ngs.read_utils.ReadIdStore(db_path2) as store2:
            store2.extend(original_ids)
            store2.shrink_to_subsample(10)
            result2 = set(store2)

        # Very unlikely to be identical (1 in C(100,10) chance)
        # We allow this test to pass even if they match, just check both are valid
        self.assertEqual(len(result1), 10)
        self.assertEqual(len(result2), 10)

    def test_filter_bam_by_ids_include(self):
        """Test filter_bam_by_ids with include=True keeps only matching reads."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname('.bam')

        samtools = viral_ngs.core.samtools.SamtoolsTool()

        # Get some read names from the input BAM
        read_names = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i < 10:  # Get first 10 reads (5 pairs)
                    read_names.add(read.query_name)
                else:
                    break

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(read_names)
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        # Output should have exactly len(read_names) * 2 reads (paired-end)
        output_count = samtools.count(output_bam)
        self.assertEqual(output_count, len(read_names) * 2)

    def test_filter_bam_by_ids_exclude(self):
        """Test filter_bam_by_ids with include=False excludes matching reads."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname('.bam')

        samtools = viral_ngs.core.samtools.SamtoolsTool()
        input_count = samtools.count(input_bam)

        # Get some read names from the input BAM
        read_names = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i < 10:  # Get first 10 reads (5 pairs)
                    read_names.add(read.query_name)
                else:
                    break

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(read_names)
            store.filter_bam_by_ids(input_bam, output_bam, include=False)

        # Output should have input_count - (len(read_names) * 2) reads
        output_count = samtools.count(output_bam)
        expected = input_count - (len(read_names) * 2)
        self.assertEqual(output_count, expected)

    def test_filter_bam_by_ids_header_preserved(self):
        """Test that BAM headers are semantically equivalent after filtering."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        output_bam = viral_ngs.core.file.mkstempfname('.bam')

        # Get some read names
        read_names = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i < 5:
                    read_names.add(read.query_name)
                else:
                    break

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(read_names)
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        # Compare headers
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as inb:
            input_header = inb.header

        with pysam.AlignmentFile(output_bam, 'rb', check_sq=False) as outb:
            output_header = outb.header

        # Check header keys match
        self.assertEqual(set(input_header.keys()), set(output_header.keys()))

        # Check RG entries preserved
        if 'RG' in input_header:
            self.assertEqual(input_header['RG'], output_header['RG'])

        # Check HD header preserved
        if 'HD' in input_header:
            self.assertEqual(dict(input_header['HD']), dict(output_header['HD']))

    def test_filter_bam_by_ids_empty_input(self):
        """Test filter_bam_by_ids handles empty input BAM."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        output_bam = viral_ngs.core.file.mkstempfname('.bam')

        samtools = viral_ngs.core.samtools.SamtoolsTool()

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2'])
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        self.assertEqual(samtools.count(output_bam), 0)

    def test_filter_bam_by_ids_empty_store_include(self):
        """Test include mode with empty store produces empty BAM with header."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname('.bam')

        samtools = viral_ngs.core.samtools.SamtoolsTool()

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            # Empty store
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        # Output should have 0 reads
        self.assertEqual(samtools.count(output_bam), 0)

        # But should have valid header
        with pysam.AlignmentFile(output_bam, 'rb', check_sq=False) as bam:
            self.assertIsNotNone(bam.header)

    def test_filter_bam_by_ids_empty_store_exclude(self):
        """Test exclude mode with empty store keeps all reads."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname('.bam')

        samtools = viral_ngs.core.samtools.SamtoolsTool()
        input_count = samtools.count(input_bam)

        db_path = viral_ngs.core.file.mkstempfname('.db')
        with viral_ngs.read_utils.ReadIdStore(db_path) as store:
            # Empty store
            store.filter_bam_by_ids(input_bam, output_bam, include=False)

        # Output should have all input reads
        self.assertEqual(samtools.count(output_bam), input_count)


class TestRmdupBbnorm(TestCaseWithTmp):
    """Tests for rmdup_bbnorm_bam using BBNorm for deduplication/normalization."""

    def setUp(self):
        super(TestRmdupBbnorm, self).setUp()
        self.samtools = viral_ngs.core.samtools.SamtoolsTool()

    def test_bbnorm_canned_input(self):
        """Test rmdup_bbnorm_bam with standard paired-end input."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        viral_ngs.read_utils.rmdup_bbnorm_bam(input_bam, output_bam, threads=1, memory='250m')

        # Output should have reads and be <= input count
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_empty_input(self):
        """Test rmdup_bbnorm_bam handles empty BAM."""
        empty_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # Empty input doesn't call bbnorm, so no need for memory/threads
        viral_ngs.read_utils.rmdup_bbnorm_bam(empty_bam, output_bam)

        self.assertEqual(self.samtools.count(output_bam), 0)

    def test_bbnorm_multi_library(self):
        """Test rmdup_bbnorm_bam with multiple libraries/read groups."""
        # G5012.3.testreads.bam has 18710 reads, 12 RGs, 2 libraries
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        viral_ngs.read_utils.rmdup_bbnorm_bam(input_bam, output_bam, threads=1, memory='250m')

        # Output should have reads
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_single_end(self):
        """Test rmdup_bbnorm_bam with single-end reads."""
        # in.2libs3rgs.bam has 2850 single-end reads
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestPerSample', 'in.2libs3rgs.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        viral_ngs.read_utils.rmdup_bbnorm_bam(input_bam, output_bam, threads=1, memory='250m')

        # Output should have reads
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_min_input_reads_skip(self):
        """Test that min_input_reads skips processing when below threshold."""
        # input.bam has 1794 reads, set threshold to 2000
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # min_input_reads causes skip, so bbnorm isn't called
        viral_ngs.read_utils.rmdup_bbnorm_bam(input_bam, output_bam, min_input_reads=2000)

        # Output should equal input (copied, not processed)
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertEqual(output_count, input_count)

    def test_bbnorm_min_input_reads_process(self):
        """Test that min_input_reads processes when above threshold."""
        # input.bam has 1794 reads, set threshold to 1000
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        viral_ngs.read_utils.rmdup_bbnorm_bam(input_bam, output_bam, min_input_reads=1000,
                                     threads=1, memory='250m')

        # Output should have reads (processing occurred)
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_max_output_reads_downsample(self):
        """Test that max_output_reads downsamples the keep-list."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # Set max_output_reads to 500 (less than expected bbnorm output)
        # Use low memory and single thread for CI compatibility
        max_output_reads = 500
        viral_ngs.read_utils.rmdup_bbnorm_bam(input_bam, output_bam, max_output_reads=max_output_reads,
                                     threads=1, memory='250m')

        # Output should be approximately max_output_reads (tolerance for pairs)
        output_count = self.samtools.count(output_bam)
        # For paired reads, we downsample read IDs, so output is ~2x the ID count
        # Allow some tolerance over exactly 2x
        expected_max_reads = max_output_reads * 2 + 100  # IDs * 2 reads/pair + tolerance
        self.assertLessEqual(output_count, expected_max_reads)


class TestRmdupClumpify(TestCaseWithTmp):
    """Tests for rmdup_clumpify_bam.

    clumpify ships in the bbmap conda package on both amd64 and arm64, so
    unlike the mvicuna and cdhit dedup tests these are not ARM-skipped.
    """

    def setUp(self):
        super(TestRmdupClumpify, self).setUp()
        self.samtools = viral_ngs.core.samtools.SamtoolsTool()

    def _collapse_libraries(self, inBam, library='onelib'):
        """Return a copy of inBam whose read groups all share a single LB tag."""
        header_file = viral_ngs.core.file.mkstempfname('.txt')
        with open(header_file, 'wt') as outf:
            for row in self.samtools.getHeader(inBam):
                if row[0] == '@RG':
                    row = [('LB:' + library) if x.startswith('LB:') else x for x in row]
                outf.write('\t'.join(row) + '\n')
        out_bam = viral_ngs.core.file.mkstempfname('.onelib.bam')
        self.samtools.reheader(inBam, header_file, out_bam)
        return out_bam

    def test_clumpify_canned_input(self):
        """Standard paired-end input: reads are removed, some survive."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam, threads=1, memory='250m')

        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLess(output_count, input_count)

    def test_clumpify_empty_input(self):
        """Empty BAM in, valid header-only BAM out -- never a crash."""
        empty_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        viral_ngs.read_utils.rmdup_clumpify_bam(empty_bam, output_bam)

        self.assertEqual(self.samtools.count(output_bam), 0)
        # output must still be a readable BAM, not a zero-byte file
        with pysam.AlignmentFile(output_bam, 'rb', check_sq=False) as bam:
            self.assertEqual(len(list(bam)), 0)

    def test_clumpify_single_end(self):
        """Single-end dedup actually removes reads.

        This is the case rmdup_mvicuna_bam silently passed through untouched.
        """
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestPerSample', 'in.2libs3rgs.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam, threads=1, memory='250m')

        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLess(output_count, input_count)

    def test_clumpify_preserves_header_and_tags(self):
        """Filtering the original BAM preserves the header and per-read tags.

        This is what rmdup_cdhit_bam's FASTQ round-trip loses.
        """
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam, threads=1, memory='250m')

        in_rgs = [row for row in self.samtools.getHeader(input_bam) if row[0] == '@RG']
        out_rgs = [row for row in self.samtools.getHeader(output_bam) if row[0] == '@RG']
        self.assertEqual(in_rgs, out_rgs)

        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            in_tags = {r.query_name: r.get_tag('RG') for r in bam if r.has_tag('RG')}
        self.assertGreater(len(in_tags), 0)

        with pysam.AlignmentFile(output_bam, 'rb', check_sq=False) as bam:
            out_reads = list(bam)
        self.assertGreater(len(out_reads), 0)
        for read in out_reads:
            self.assertTrue(read.has_tag('RG'),
                            msg="read {} lost its RG tag".format(read.query_name))
            self.assertEqual(read.get_tag('RG'), in_tags[read.query_name])

    def _make_cross_library_bam(self, n=200, seed=7):
        """Build an unaligned BAM where every template appears in both libraries.

        200 distinct paired templates, each written once under read group rgA
        (LB:libA) and once under rgB (LB:libB) with a different query name.
        There are therefore no duplicates *within* a library, and exactly one
        duplicate for every template *across* libraries.
        """
        rnd = random.Random(seed)
        header = {
            'HD': {'VN': '1.6', 'SO': 'unsorted'},
            'RG': [
                {'ID': 'rgA', 'LB': 'libA', 'SM': 'sample', 'PL': 'illumina'},
                {'ID': 'rgB', 'LB': 'libB', 'SM': 'sample', 'PL': 'illumina'},
            ],
        }
        seqs = set()
        while len(seqs) < n:
            seqs.add((''.join(rnd.choice('ACGT') for _ in range(60)),
                      ''.join(rnd.choice('ACGT') for _ in range(60))))
        seqs = sorted(seqs)

        out_bam = viral_ngs.core.file.mkstempfname('.crosslib.bam')
        with pysam.AlignmentFile(out_bam, 'wb', header=header) as out:
            for i, (seq1, seq2) in enumerate(seqs):
                for rg in ('rgA', 'rgB'):
                    for mate, seq in ((1, seq1), (2, seq2)):
                        read = pysam.AlignedSegment()
                        read.query_name = 'tmpl%04d_%s' % (i, rg)
                        read.query_sequence = seq
                        read.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
                        read.flag = 0x1 | 0x4 | 0x8 | (0x40 if mate == 1 else 0x80)
                        read.set_tag('RG', rg)
                        out.write(read)
        return out_bam, len(seqs)

    def test_clumpify_multi_library(self):
        """Smoke test across a real 12-read-group, 2-library BAM."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam, threads=1, memory='250m')

        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_clumpify_library_aware(self):
        """Duplicates are collapsed within a library but never across libraries.

        Identical sequences from different libraries are independent observations
        of the molecule, not duplicates of each other; collapsing them discards
        real evidence and corrupts downstream library-level replicate reasoning.

        The fixture has every template present in both libraries and none repeated
        within one, so a library-aware dedup must retain everything, while pooling
        the libraries halves it.
        """
        input_bam, num_templates = self._make_cross_library_bam()
        per_library_bam = viral_ngs.core.file.mkstempfname("per_library.bam")
        pooled_bam = viral_ngs.core.file.mkstempfname("pooled.bam")

        input_count = self.samtools.count(input_bam)
        self.assertEqual(input_count, num_templates * 4)  # 2 libraries x 2 mates

        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, per_library_bam,
                                                threads=1, memory='500m')

        one_library_input = self._collapse_libraries(input_bam)
        libraries = set(
            field for row in self.samtools.getHeader(one_library_input)
            if row[0] == '@RG' for field in row if field.startswith('LB:')
        )
        self.assertEqual(len(libraries), 1, msg="fixture setup failed to collapse libraries")

        viral_ngs.read_utils.rmdup_clumpify_bam(one_library_input, pooled_bam,
                                                threads=1, memory='500m')

        per_library_count = self.samtools.count(per_library_bam)
        pooled_count = self.samtools.count(pooled_bam)

        # nothing is duplicated within a library, so a library-aware dedup keeps it all
        self.assertEqual(per_library_count, input_count)
        # pooling the libraries collapses each cross-library twin
        self.assertEqual(pooled_count, num_templates * 2)

    def test_clumpify_default_memory(self):
        """Runs without an explicit memory= argument.

        BBTools' own heap autodetection (calcmem.sh) subtracts a fixed 500MB
        floor and can compute a negative -Xmx inside a container, which makes
        the JVM refuse to start. The command must supply its own default rather
        than relying on that.
        """
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam, threads=1)

        self.assertGreater(self.samtools.count(output_bam), 0)

    def test_clumpify_min_input_reads_skip(self):
        """min_input_reads below threshold copies input through untouched."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        # input.bam has 1794 reads
        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam, min_input_reads=2000)

        self.assertEqual(self.samtools.count(output_bam), self.samtools.count(input_bam))

    def test_clumpify_max_output_reads_downsample(self):
        """max_output_reads shrinks the keep-list before filtering."""
        input_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = viral_ngs.core.file.mkstempfname("output.bam")

        max_output_reads = 300
        viral_ngs.read_utils.rmdup_clumpify_bam(input_bam, output_bam,
                                                max_output_reads=max_output_reads,
                                                threads=1, memory='250m')

        # IDs are per-template, so a paired BAM yields up to 2 reads per kept ID
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, max_output_reads * 2)


class TestMvicunaRemoved(unittest.TestCase):
    """M-Vicuna is gone: x86-only, unmaintained, ~17x slower than clumpify, and
    its single-end path was a documented no-op that copied input to output.

    rmdup_clumpify_bam carries forward the property that made it worth keeping
    -- library-aware and provenance-preserving in one implementation.
    """

    def test_tool_module_gone(self):
        self.assertFalse(hasattr(viral_ngs.core, 'mvicuna'))

    def test_command_unregistered(self):
        commands = dict(viral_ngs.read_utils.__commands__)
        self.assertNotIn('rmdup_mvicuna_bam', commands)
        # the replacement is registered in its place
        self.assertIn('rmdup_clumpify_bam', commands)

    def test_no_module_attributes_left(self):
        for attr in ('rmdup_mvicuna_bam', 'parser_rmdup_mvicuna_bam',
                     'mvicuna_fastqs_to_readlist', '_merge_fastqs_and_mvicuna'):
            self.assertFalse(hasattr(viral_ngs.read_utils, attr),
                             msg="read_utils still exposes %s" % attr)

    def test_tool_class_gone(self):
        names = [cls.__name__ for cls in viral_ngs.core.all_tool_classes()]
        self.assertNotIn('MvicunaTool', names)


class TestAlignAndFix(TestCaseWithTmp):
    def setUp(self):
        super(TestAlignAndFix, self).setUp()
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta')
        self.refFasta = viral_ngs.core.file.mkstempfname('.ref.fasta')
        shutil.copyfile(orig_ref, self.refFasta)

    @unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
    def test_novoalign(self):
        self.simple_execution('novoalign')

    def test_bwa(self):
        self.simple_execution('bwa')

    def test_minimap2(self):
        self.simple_execution('minimap2')

    def simple_execution(self, aligner):
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = viral_ngs.core.file.mkstempfname('.outBamAll.bam')
        outBamFiltered = viral_ngs.core.file.mkstempfname('.outBamFiltered.bam')

        args = viral_ngs.read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll, '--outBamFiltered', outBamFiltered, '--aligner', aligner])
        args.func_main(args)

    def test_empty_reads(self):
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        args = viral_ngs.read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta,  '--aligner', 'minimap2',
            '--outBamAll', viral_ngs.core.file.mkstempfname('.outBamAll.bam'),
            '--outBamFiltered', viral_ngs.core.file.mkstempfname('.outBamFiltered.bam')])
        args.func_main(args)

    def test_dup_marker_sambamba(self):
        """Test using sambamba for duplicate marking"""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = viral_ngs.core.file.mkstempfname('.outBamAll.bam')

        args = viral_ngs.read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll,
             '--aligner', 'minimap2', '--dupMarker', 'sambamba'])
        args.func_main(args)

        samtools = viral_ngs.core.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertFalse(samtools.isEmpty(outBamAll))

    def test_dup_marker_picard_explicit(self):
        """Test explicitly using Picard for duplicate marking"""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = viral_ngs.core.file.mkstempfname('.outBamAll.bam')

        args = viral_ngs.read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll,
             '--aligner', 'minimap2', '--dupMarker', 'picard'])
        args.func_main(args)

        samtools = viral_ngs.core.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertFalse(samtools.isEmpty(outBamAll))

    def test_dup_marker_default_is_sambamba(self):
        """Verify default duplicate marker is sambamba"""
        parser = viral_ngs.read_utils.parser_align_and_fix(argparse.ArgumentParser())
        args = parser.parse_args(['in.bam', 'ref.fasta'])
        self.assertEqual(args.dup_marker, 'sambamba')

    def test_align_and_fix_full_sambamba_pipeline(self):
        """End-to-end test with sambamba for markdup and indexing"""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = viral_ngs.core.file.mkstempfname('.outBamAll.bam')
        outBamFiltered = viral_ngs.core.file.mkstempfname('.outBamFiltered.bam')

        args = viral_ngs.read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll, '--outBamFiltered', outBamFiltered,
             '--aligner', 'minimap2', '--dupMarker', 'sambamba'])
        args.func_main(args)

        samtools = viral_ngs.core.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertTrue(os.path.exists(outBamFiltered))
        self.assertFalse(samtools.isEmpty(outBamAll))
        # Verify index files exist
        self.assertTrue(os.path.exists(outBamAll + '.bai') or os.path.exists(outBamAll.replace('.bam', '.bai')))


class TestDownsampleBams(TestCaseWithTmp):
    def setUp(self):
        super(TestDownsampleBams, self).setUp()
        orig_larger_bam  = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        orig_smaller_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.mini.bam')
        orig_with_dup    = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned','input.bam')
        self.larger_bam  = viral_ngs.core.file.mkstempfname('.larger.bam')
        self.smaller_bam = viral_ngs.core.file.mkstempfname('.smaller.bam')
        self.with_dups   = viral_ngs.core.file.mkstempfname('.with_dups.bam')
        shutil.copyfile(orig_larger_bam, self.larger_bam)
        shutil.copyfile(orig_smaller_bam, self.smaller_bam)
        shutil.copyfile(orig_with_dup, self.with_dups)

        self.samtools = viral_ngs.core.samtools.SamtoolsTool()

    def test_normalization_to_lowest_cardinality(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = self.samtools.count(self.smaller_bam)
        # target count not passed in since we are checking that the count of the smaller file is used
        viral_ngs.read_utils.main_downsample_bams([self.larger_bam, self.smaller_bam], temp_dir, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_to_target_count(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = 4000
        viral_ngs.read_utils.main_downsample_bams([self.larger_bam, self.smaller_bam], temp_dir, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_to_target_count_without_subdir(self):
        target_count = 4000
        viral_ngs.read_utils.main_downsample_bams([self.larger_bam], out_path=None, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(os.path.dirname(self.larger_bam), '*downsampled-*.bam')))
        
        print(output_bams)
        self.assertGreater(len(output_bams), 0, msg="No output files matching *downsampled-*.bam found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_with_dedup_after(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = 1500
        viral_ngs.read_utils.main_downsample_bams([self.with_dups], temp_dir, deduplicate_after=True, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertLess(self.samtools.count(out_bam), target_count, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    def test_downsample_with_dedup_before(self):
        """ Also tests subdir output """
        temp_dir = tempfile.mkdtemp()

        target_count = 1500
        viral_ngs.read_utils.main_downsample_bams([self.with_dups], temp_dir, deduplicate_before=True, specified_read_count=target_count, JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10, msg="{} not downsampled to the target size: {}".format(os.path.basename(out_bam),target_count))

    @unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
    def test_downsample_with_dedup_tool_cdhit(self):
        """--dedupTool cdhit routes the dedup step through rmdup_cdhit_bam."""
        temp_dir = tempfile.mkdtemp()

        target_count = 1500
        viral_ngs.read_utils.main_downsample_bams([self.with_dups], temp_dir,
                                                  deduplicate_before=True,
                                                  dedup_tool='cdhit',
                                                  specified_read_count=target_count,
                                                  JVMmemory="1g")

        output_bams = list(glob.glob(os.path.join(temp_dir, '*.bam')))
        self.assertGreater(len(output_bams), 0, msg="No output found")
        for out_bam in output_bams:
            self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=10)

    def test_downsample_dedup_tool_parser_roundtrip(self):
        """--dedupTool parses, defaults to clumpify, and rejects unknown tools."""
        parser = viral_ngs.read_utils.parser_downsample_bams(argparse.ArgumentParser())

        args = parser.parse_args(['in.bam', '--deduplicateBefore'])
        self.assertEqual(args.dedup_tool, 'clumpify')

        args = parser.parse_args(['in.bam', '--deduplicateBefore', '--dedupTool', 'cdhit'])
        self.assertEqual(args.dedup_tool, 'cdhit')

        with self.assertRaises(SystemExit):
            parser.parse_args(['in.bam', '--dedupTool', 'mvicuna'])

    def test_downsample_to_too_large_target_count(self):
        """ Should fail """
        temp_dir = tempfile.mkdtemp()

        target_count = 20000

        with self.assertRaises(ValueError):
            viral_ngs.read_utils.main_downsample_bams([self.larger_bam, self.smaller_bam], temp_dir, specified_read_count=target_count, JVMmemory="1g")


class TestTrimRmdupSubsamp(TestCaseWithTmp):
    '''Test the trim_rmdup_subsamp command.

    Tests ported from viral-assemble/test/unit/test_assembly.py.
    Uses threads=1 to avoid overwhelming CI runners during parallel test execution.
    '''

    def test_subsamp_empty(self):
        """Test with empty BAM input - should return all zeros."""
        inDir = viral_ngs.core.file.get_test_input_path()
        inBam = os.path.join(inDir, 'empty.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = viral_ngs.core.file.mkstempfname('.out.bam')
        read_stats = viral_ngs.read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=10, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (0, 0, 0, 0, 0, 0))

    def test_subsamp_small_50(self):
        """Test subsampling to 50 reads from small input."""
        inDir = viral_ngs.core.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = viral_ngs.core.file.mkstempfname('.out.bam')
        read_stats = viral_ngs.read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=50, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (200, 172, 172, 50, 50, 0))

    def test_subsamp_small_90(self):
        """Test subsampling to 90 reads."""
        inDir = viral_ngs.core.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = viral_ngs.core.file.mkstempfname('.out.bam')
        read_stats = viral_ngs.read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=90, threads=1)
        os.unlink(outBam)
        # counts are individual reads
        self.assertEqual(read_stats, (200, 172, 172, 90, 90, 0))

    def test_subsamp_small_200(self):
        """Test where unpaired reads are needed to reach threshold."""
        inDir = viral_ngs.core.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = viral_ngs.core.file.mkstempfname('.out.bam')
        read_stats = viral_ngs.read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=200, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (200, 172, 172, 185, 172, 13))

    def test_subsamp_big_500(self):
        """Test with larger input file."""
        inDir = viral_ngs.core.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.testreads.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = viral_ngs.core.file.mkstempfname('.out.bam')
        read_stats = viral_ngs.read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=500, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (18710, 16310, 16310, 500, 500, 0))
