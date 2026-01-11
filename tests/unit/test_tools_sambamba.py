# Unit tests for tools.sambamba

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import tempfile
import shutil
import platform

import pytest

import util.file
import tools
import tools.sambamba
import tools.samtools
from test import TestCaseWithTmp


# Skip merge tests on ARM platforms where x86 emulation causes segfaults
IS_ARM = platform.machine() in ('arm64', 'aarch64')
SKIP_MERGE_REASON = "sambamba merge segfaults under x86 emulation on ARM"


class TestToolSambamba(TestCaseWithTmp):
    """Tests for SambambaTool - high-performance BAM processing"""

    def setUp(self):
        super().setUp()
        self.sambamba = tools.sambamba.SambambaTool()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_sambamba_installed(self):
        """Verify sambamba is installed and accessible"""
        self.sambamba.install()
        self.assertTrue(os.path.isfile(self.sambamba.install_and_get_path()))

    def test_version(self):
        """Verify version detection works"""
        self.sambamba.install()
        version = self.sambamba.version()
        self.assertIsNotNone(version)
        self.assertIn('.', version)  # Version should contain a dot (e.g., "1.0.1")


class TestSambambaSort(TestCaseWithTmp):
    """Tests for sambamba sort functionality"""

    def setUp(self):
        super().setUp()
        self.sambamba = tools.sambamba.SambambaTool()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_sort_coordinate(self):
        """Test coordinate sort (default)"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBam = util.file.mkstempfname('.sorted.bam')

        self.sambamba.sort(inBam, outBam)

        # Verify output exists and is non-empty
        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

        # Verify read count matches input
        expected_count = self.samtools.count(inBam)
        actual_count = self.samtools.count(outBam)
        self.assertEqual(actual_count, expected_count)

    def test_sort_queryname(self):
        """Test queryname sort"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBam = util.file.mkstempfname('.namesorted.bam')

        self.sambamba.sort(inBam, outBam, sort_order='queryname')

        # Verify output exists
        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

        # Verify read count is preserved
        expected_count = self.samtools.count(inBam)
        actual_count = self.samtools.count(outBam)
        self.assertEqual(actual_count, expected_count)

    def test_sort_with_threads(self):
        """Test multi-threaded sorting"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBam = util.file.mkstempfname('.sorted.bam')

        # Should work with explicit thread count
        self.sambamba.sort(inBam, outBam, threads=2)

        self.assertTrue(os.path.exists(outBam))
        expected_count = self.samtools.count(inBam)
        actual_count = self.samtools.count(outBam)
        self.assertEqual(actual_count, expected_count)


class TestSambambaIndex(TestCaseWithTmp):
    """Tests for sambamba index functionality"""

    def setUp(self):
        super().setUp()
        self.sambamba = tools.sambamba.SambambaTool()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_index_basic(self):
        """Test creating .bai index file"""
        # First, sort the BAM (indexing requires coordinate-sorted BAM)
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        # Now index it
        self.sambamba.index(sortedBam)

        # Verify index file was created
        indexFile = sortedBam + '.bai'
        self.assertTrue(os.path.exists(indexFile))
        self.assertGreater(os.path.getsize(indexFile), 0)

    def test_index_with_threads(self):
        """Test multi-threaded indexing"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        self.sambamba.index(sortedBam, threads=2)

        indexFile = sortedBam + '.bai'
        self.assertTrue(os.path.exists(indexFile))

    def test_index_file_exists(self):
        """Verify .bai file is created in expected location"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        # Remove any existing index
        indexFile = sortedBam + '.bai'
        if os.path.exists(indexFile):
            os.unlink(indexFile)

        self.sambamba.index(sortedBam)

        self.assertTrue(os.path.exists(indexFile))


class TestSambambaMerge(TestCaseWithTmp):
    """Tests for sambamba merge functionality"""

    def setUp(self):
        super().setUp()
        self.sambamba = tools.sambamba.SambambaTool()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_merge_two_bams(self):
        """Test merging two BAM files"""
        # Sort both input BAMs first (merge requires sorted input)
        inBam1 = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        inBam2 = os.path.join(util.file.get_test_input_path(), 'G5012.3.mini.bam')

        sortedBam1 = util.file.mkstempfname('.sorted1.bam')
        sortedBam2 = util.file.mkstempfname('.sorted2.bam')
        self.sambamba.sort(inBam1, sortedBam1)
        self.sambamba.sort(inBam2, sortedBam2)

        outBam = util.file.mkstempfname('.merged.bam')
        self.sambamba.merge([sortedBam1, sortedBam2], outBam)

        # Verify output exists
        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

    def test_merge_preserves_reads(self):
        """Verify read count equals sum of inputs"""
        inBam1 = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        inBam2 = os.path.join(util.file.get_test_input_path(), 'G5012.3.mini.bam')

        sortedBam1 = util.file.mkstempfname('.sorted1.bam')
        sortedBam2 = util.file.mkstempfname('.sorted2.bam')
        self.sambamba.sort(inBam1, sortedBam1)
        self.sambamba.sort(inBam2, sortedBam2)

        count1 = self.samtools.count(sortedBam1)
        count2 = self.samtools.count(sortedBam2)

        outBam = util.file.mkstempfname('.merged.bam')
        self.sambamba.merge([sortedBam1, sortedBam2], outBam)

        merged_count = self.samtools.count(outBam)
        self.assertEqual(merged_count, count1 + count2)


class TestSambambaFlagstat(TestCaseWithTmp):
    """Tests for sambamba flagstat functionality"""

    def setUp(self):
        super().setUp()
        self.sambamba = tools.sambamba.SambambaTool()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_flagstat_basic(self):
        """Test getting alignment statistics"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')

        stats = self.sambamba.flagstat(inBam)

        # Should return a dictionary with statistics
        self.assertIsInstance(stats, dict)

    def test_flagstat_returns_dict(self):
        """Verify parsed output structure"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')

        stats = self.sambamba.flagstat(inBam)

        # Should contain standard flagstat keys
        self.assertIn('total', stats)

    def test_flagstat_empty_bam(self):
        """Test flagstat on empty BAM file"""
        inBam = os.path.join(util.file.get_test_input_path(), 'empty.bam')

        stats = self.sambamba.flagstat(inBam)

        self.assertIsInstance(stats, dict)
        self.assertEqual(stats.get('total', 0), 0)


class TestSambambaMarkdup(TestCaseWithTmp):
    """Tests for sambamba markdup functionality"""

    def setUp(self):
        super().setUp()
        self.sambamba = tools.sambamba.SambambaTool()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_markdup_basic(self):
        """Test basic duplicate marking on aligned BAM"""
        # Use the larger test file which should have duplicates
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')

        # Sort first (markdup requires coordinate-sorted input)
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        outBam = util.file.mkstempfname('.deduped.bam')
        self.sambamba.markdup(sortedBam, outBam)

        # Verify output exists
        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

        # Read count should be preserved (duplicates are marked, not removed)
        in_count = self.samtools.count(sortedBam)
        out_count = self.samtools.count(outBam)
        self.assertEqual(out_count, in_count)

    def test_markdup_with_threads(self):
        """Test multi-threaded duplicate marking"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        outBam = util.file.mkstempfname('.deduped.bam')
        self.sambamba.markdup(sortedBam, outBam, threads=2)

        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

    def test_markdup_empty_input(self):
        """Test handling of empty BAM file"""
        inBam = os.path.join(util.file.get_test_input_path(), 'empty.bam')

        # Sort the empty BAM (should still work)
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        outBam = util.file.mkstempfname('.deduped.bam')
        self.sambamba.markdup(sortedBam, outBam)

        # Output should exist and be empty
        self.assertTrue(os.path.exists(outBam))
        self.assertTrue(self.samtools.isEmpty(outBam))

    def test_markdup_removes_duplicates(self):
        """Test remove_duplicates=True flag"""
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        sortedBam = util.file.mkstempfname('.sorted.bam')
        self.sambamba.sort(inBam, sortedBam)

        outBam = util.file.mkstempfname('.deduped.bam')
        self.sambamba.markdup(sortedBam, outBam, remove_duplicates=True)

        # Output should exist
        self.assertTrue(os.path.exists(outBam))

        # When removing duplicates, output should have fewer or equal reads
        in_count = self.samtools.count(sortedBam)
        out_count = self.samtools.count(outBam)
        self.assertLessEqual(out_count, in_count)
