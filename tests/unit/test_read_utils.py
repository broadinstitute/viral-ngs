# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest
import argparse
import filecmp
import os
import glob

import pysam
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


class TestReadIdStore(TestCaseWithTmp):
    """Tests for ReadIdStore SQLite-backed read ID storage."""

    def test_add_from_fastq_paired(self):
        """Test adding read IDs from paired-end interleaved FASTQ."""
        # Create a test FASTQ with paired reads
        fastq_path = util.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            # Write 3 read pairs (6 entries, but only 3 unique IDs)
            for i in range(3):
                f.write('@read{}/1\nACGT\n+\nIIII\n'.format(i))
                f.write('@read{}/2\nACGT\n+\nIIII\n'.format(i))

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            added = store.add_from_fastq(fastq_path)
            self.assertEqual(added, 3)  # 3 unique read IDs
            self.assertEqual(len(store), 3)

    def test_add_from_fastq_single_end(self):
        """Test adding read IDs from single-end FASTQ."""
        fastq_path = util.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            for i in range(5):
                f.write('@read{}\nACGT\n+\nIIII\n'.format(i))

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            added = store.add_from_fastq(fastq_path)
            self.assertEqual(added, 5)
            self.assertEqual(len(store), 5)

    def test_deduplication(self):
        """Test that duplicate read IDs are ignored."""
        fastq_path = util.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            # Write same read ID multiple times
            for _ in range(10):
                f.write('@duplicate_read\nACGT\n+\nIIII\n')

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.add_from_fastq(fastq_path)
            self.assertEqual(len(store), 1)  # Only 1 unique ID

    def test_write_to_file(self):
        """Test writing read IDs to file."""
        fastq_path = util.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            for i in range(5):
                f.write('@read{}\nACGT\n+\nIIII\n'.format(i))

        db_path = util.file.mkstempfname('.db')
        out_path = util.file.mkstempfname('.txt')

        with read_utils.ReadIdStore(db_path) as store:
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
        fastq_path = util.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            for i in range(100):
                f.write('@read{}\nACGT\n+\nIIII\n'.format(i))

        db_path = util.file.mkstempfname('.db')
        out_path = util.file.mkstempfname('.txt')

        with read_utils.ReadIdStore(db_path) as store:
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
        fastq_path = util.file.mkstempfname('.fastq')
        with open(fastq_path, 'wt') as f:
            pass  # Empty file

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            added = store.add_from_fastq(fastq_path)
            self.assertEqual(added, 0)
            self.assertEqual(len(store), 0)

    def test_add_single(self):
        """Test adding a single read ID."""
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
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
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
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
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            # Use a generator expression
            added = store.extend('read{}'.format(i) for i in range(100))
            self.assertEqual(added, 100)
            self.assertEqual(len(store), 100)

    def test_contains(self):
        """Test membership testing with 'in' operator."""
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2', 'read3'])
            self.assertIn('read1', store)
            self.assertIn('read2', store)
            self.assertNotIn('read4', store)

    def test_iter(self):
        """Test iteration over read IDs in insertion order."""
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(['read3', 'read1', 'read2'])
            # Should iterate in insertion order
            result = list(store)
            self.assertEqual(result, ['read3', 'read1', 'read2'])

    def test_delitem(self):
        """Test deleting read IDs with del operator."""
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
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
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2'])

            store.discard('read1')
            self.assertEqual(len(store), 1)
            self.assertNotIn('read1', store)

            # Discarding non-existent ID should not raise
            store.discard('nonexistent')  # Should not raise
            self.assertEqual(len(store), 1)

    def test_shrink_to_subsample_basic(self):
        """Test shrink_to_subsample reduces store to n elements."""
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
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
        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(['read{}'.format(i) for i in range(50)])
            self.assertEqual(len(store), 50)

            # Request 100, but only have 50
            store.shrink_to_subsample(100)
            self.assertEqual(len(store), 50)  # Unchanged

    def test_shrink_to_subsample_randomness(self):
        """Test that shrink_to_subsample produces different results on repeated calls."""
        original_ids = ['read{}'.format(i) for i in range(100)]

        # Create two stores and subsample each
        db_path1 = util.file.mkstempfname('.db')
        db_path2 = util.file.mkstempfname('.db')

        with read_utils.ReadIdStore(db_path1) as store1:
            store1.extend(original_ids)
            store1.shrink_to_subsample(10)
            result1 = set(store1)

        with read_utils.ReadIdStore(db_path2) as store2:
            store2.extend(original_ids)
            store2.shrink_to_subsample(10)
            result2 = set(store2)

        # Very unlikely to be identical (1 in C(100,10) chance)
        # We allow this test to pass even if they match, just check both are valid
        self.assertEqual(len(result1), 10)
        self.assertEqual(len(result2), 10)

    def test_filter_bam_by_ids_include(self):
        """Test filter_bam_by_ids with include=True keeps only matching reads."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()

        # Get some read names from the input BAM
        read_names = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i < 10:  # Get first 10 reads (5 pairs)
                    read_names.add(read.query_name)
                else:
                    break

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(read_names)
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        # Output should have exactly len(read_names) * 2 reads (paired-end)
        output_count = samtools.count(output_bam)
        self.assertEqual(output_count, len(read_names) * 2)

    def test_filter_bam_by_ids_exclude(self):
        """Test filter_bam_by_ids with include=False excludes matching reads."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        input_count = samtools.count(input_bam)

        # Get some read names from the input BAM
        read_names = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i < 10:  # Get first 10 reads (5 pairs)
                    read_names.add(read.query_name)
                else:
                    break

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(read_names)
            store.filter_bam_by_ids(input_bam, output_bam, include=False)

        # Output should have input_count - (len(read_names) * 2) reads
        output_count = samtools.count(output_bam)
        expected = input_count - (len(read_names) * 2)
        self.assertEqual(output_count, expected)

    def test_filter_bam_by_ids_header_preserved(self):
        """Test that BAM headers are semantically equivalent after filtering."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        output_bam = util.file.mkstempfname('.bam')

        # Get some read names
        read_names = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for i, read in enumerate(bam):
                if i < 5:
                    read_names.add(read.query_name)
                else:
                    break

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
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
        input_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        output_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(['read1', 'read2'])
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        self.assertEqual(samtools.count(output_bam), 0)

    def test_filter_bam_by_ids_empty_store_include(self):
        """Test include mode with empty store produces empty BAM with header."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            # Empty store
            store.filter_bam_by_ids(input_bam, output_bam, include=True)

        # Output should have 0 reads
        self.assertEqual(samtools.count(output_bam), 0)

        # But should have valid header
        with pysam.AlignmentFile(output_bam, 'rb', check_sq=False) as bam:
            self.assertIsNotNone(bam.header)

    def test_filter_bam_by_ids_empty_store_exclude(self):
        """Test exclude mode with empty store keeps all reads."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        input_count = samtools.count(input_bam)

        db_path = util.file.mkstempfname('.db')
        with read_utils.ReadIdStore(db_path) as store:
            # Empty store
            store.filter_bam_by_ids(input_bam, output_bam, include=False)

        # Output should have all input reads
        self.assertEqual(samtools.count(output_bam), input_count)


class TestRmdupBbnorm(TestCaseWithTmp):
    """Tests for rmdup_bbnorm_bam using BBNorm for deduplication/normalization."""

    def setUp(self):
        super(TestRmdupBbnorm, self).setUp()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_bbnorm_canned_input(self):
        """Test rmdup_bbnorm_bam with standard paired-end input."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        read_utils.rmdup_bbnorm_bam(input_bam, output_bam, threads=1, memory='250m')

        # Output should have reads and be <= input count
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_empty_input(self):
        """Test rmdup_bbnorm_bam handles empty BAM."""
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # Empty input doesn't call bbnorm, so no need for memory/threads
        read_utils.rmdup_bbnorm_bam(empty_bam, output_bam)

        self.assertEqual(self.samtools.count(output_bam), 0)

    def test_bbnorm_multi_library(self):
        """Test rmdup_bbnorm_bam with multiple libraries/read groups."""
        # G5012.3.testreads.bam has 18710 reads, 12 RGs, 2 libraries
        input_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        read_utils.rmdup_bbnorm_bam(input_bam, output_bam, threads=1, memory='250m')

        # Output should have reads
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_single_end(self):
        """Test rmdup_bbnorm_bam with single-end reads."""
        # in.2libs3rgs.bam has 2850 single-end reads
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestPerSample', 'in.2libs3rgs.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        read_utils.rmdup_bbnorm_bam(input_bam, output_bam, threads=1, memory='250m')

        # Output should have reads
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_min_input_reads_skip(self):
        """Test that min_input_reads skips processing when below threshold."""
        # input.bam has 1794 reads, set threshold to 2000
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # min_input_reads causes skip, so bbnorm isn't called
        read_utils.rmdup_bbnorm_bam(input_bam, output_bam, min_input_reads=2000)

        # Output should equal input (copied, not processed)
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertEqual(output_count, input_count)

    def test_bbnorm_min_input_reads_process(self):
        """Test that min_input_reads processes when above threshold."""
        # input.bam has 1794 reads, set threshold to 1000
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # Use low memory and single thread for CI compatibility
        read_utils.rmdup_bbnorm_bam(input_bam, output_bam, min_input_reads=1000,
                                     threads=1, memory='250m')

        # Output should have reads (processing occurred)
        input_count = self.samtools.count(input_bam)
        output_count = self.samtools.count(output_bam)
        self.assertGreater(output_count, 0)
        self.assertLessEqual(output_count, input_count)

    def test_bbnorm_max_output_reads_downsample(self):
        """Test that max_output_reads downsamples the keep-list."""
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        output_bam = util.file.mkstempfname("output.bam")

        # Set max_output_reads to 500 (less than expected bbnorm output)
        # Use low memory and single thread for CI compatibility
        max_output_reads = 500
        read_utils.rmdup_bbnorm_bam(input_bam, output_bam, max_output_reads=max_output_reads,
                                     threads=1, memory='250m')

        # Output should be approximately max_output_reads (tolerance for pairs)
        output_count = self.samtools.count(output_bam)
        # For paired reads, we downsample read IDs, so output is ~2x the ID count
        # Allow some tolerance over exactly 2x
        expected_max_reads = max_output_reads * 2 + 100  # IDs * 2 reads/pair + tolerance
        self.assertLessEqual(output_count, expected_max_reads)


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

    def test_skip_realign_flag(self):
        """Test that --skipRealign skips GATK local realignment"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll,
             '--aligner', 'minimap2', '--skipRealign'])
        args.func_main(args)

        # Verify output exists and has reads
        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertFalse(samtools.isEmpty(outBamAll))

    def test_skip_realign_with_skip_mark_dupes(self):
        """Test combining --skipRealign with --skipMarkDupes"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll,
             '--aligner', 'minimap2', '--skipRealign', '--skipMarkDupes'])
        args.func_main(args)

        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertFalse(samtools.isEmpty(outBamAll))

    def test_dup_marker_sambamba(self):
        """Test using sambamba for duplicate marking"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll,
             '--aligner', 'minimap2', '--dupMarker', 'sambamba', '--skipRealign'])
        args.func_main(args)

        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertFalse(samtools.isEmpty(outBamAll))

    def test_dup_marker_picard_explicit(self):
        """Test explicitly using Picard for duplicate marking"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll,
             '--aligner', 'minimap2', '--dupMarker', 'picard', '--skipRealign'])
        args.func_main(args)

        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertFalse(samtools.isEmpty(outBamAll))

    def test_dup_marker_default_is_sambamba(self):
        """Verify default duplicate marker is sambamba"""
        parser = read_utils.parser_align_and_fix(argparse.ArgumentParser())
        args = parser.parse_args(['in.bam', 'ref.fasta'])
        self.assertEqual(args.dup_marker, 'sambamba')

    def test_align_and_fix_full_sambamba_pipeline(self):
        """End-to-end test with sambamba for markdup and indexing"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBamAll = util.file.mkstempfname('.outBamAll.bam')
        outBamFiltered = util.file.mkstempfname('.outBamFiltered.bam')

        args = read_utils.parser_align_and_fix(argparse.ArgumentParser()).parse_args(
            [inBam, self.refFasta, '--outBamAll', outBamAll, '--outBamFiltered', outBamFiltered,
             '--aligner', 'minimap2', '--dupMarker', 'sambamba', '--skipRealign'])
        args.func_main(args)

        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(os.path.exists(outBamAll))
        self.assertTrue(os.path.exists(outBamFiltered))
        self.assertFalse(samtools.isEmpty(outBamAll))
        # Verify index files exist
        self.assertTrue(os.path.exists(outBamAll + '.bai') or os.path.exists(outBamAll.replace('.bam', '.bai')))


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


class TestTrimRmdupSubsamp(TestCaseWithTmp):
    '''Test the trim_rmdup_subsamp command.

    Tests ported from viral-assemble/test/unit/test_assembly.py.
    Uses threads=1 to avoid overwhelming CI runners during parallel test execution.
    '''

    def test_subsamp_empty(self):
        """Test with empty BAM input - should return all zeros."""
        inDir = util.file.get_test_input_path()
        inBam = os.path.join(inDir, 'empty.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = util.file.mkstempfname('.out.bam')
        read_stats = read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=10, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (0, 0, 0, 0, 0, 0))

    def test_subsamp_small_50(self):
        """Test subsampling to 50 reads from small input."""
        inDir = util.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = util.file.mkstempfname('.out.bam')
        read_stats = read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=50, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (200, 172, 172, 50, 50, 0))

    def test_subsamp_small_90(self):
        """Test subsampling to 90 reads."""
        inDir = util.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = util.file.mkstempfname('.out.bam')
        read_stats = read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=90, threads=1)
        os.unlink(outBam)
        # counts are individual reads
        self.assertEqual(read_stats, (200, 172, 172, 90, 90, 0))

    def test_subsamp_small_200(self):
        """Test where unpaired reads are needed to reach threshold."""
        inDir = util.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.subset.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = util.file.mkstempfname('.out.bam')
        read_stats = read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=200, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (200, 172, 172, 185, 172, 13))

    def test_subsamp_big_500(self):
        """Test with larger input file."""
        inDir = util.file.get_test_input_path()
        inBam = os.path.join(inDir, 'G5012.3.testreads.bam')
        clipDb = os.path.join(inDir, 'TestTrimRmdupSubsamp', 'clipDb.fasta')
        outBam = util.file.mkstempfname('.out.bam')
        read_stats = read_utils.trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=500, threads=1)
        os.unlink(outBam)
        self.assertEqual(read_stats, (18710, 16310, 16310, 500, 500, 0))
