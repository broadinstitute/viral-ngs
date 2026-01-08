# Unit tests for bbmap aligner and related tools

__author__ = "ilya@broadinstitute.org"

import unittest
import os.path
import shutil
import util.file
import tools.bbmap
import tools.samtools
import pysam
from test import TestCaseWithTmp, assert_md5_equal_to_line_in_file


class TestToolBBMap(TestCaseWithTmp):

    def setUp(self):
        super(TestToolBBMap, self).setUp()
        self.bbmap = tools.bbmap.BBMapTool()
        self.bbmap.install()
        self.samtools = tools.samtools.SamtoolsTool()

    def test_align(self):
        orig_ref = os.path.join(util.file.get_test_input_path(), 'ebola.fasta')
        inRef = util.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        reads = os.path.join(util.file.get_test_input_path(self), 'ebov_reads.bam')
        outBam = util.file.mkstempfname('.bam')
        # Use conservative memory and single thread for CI compatibility
        # BBMap alignment needs more memory than bbnorm for index building
        self.bbmap.align(inBam=reads, refFasta=inRef, outBam=outBam,
                         threads=1, Xmx='1g')
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))

    def test_bbnorm_paired_interleaved(self):
        """Test bbnorm on paired-end interleaved FASTQ file."""
        # Create interleaved FASTQ from TestRmdupUnaligned BAM
        inBam = os.path.join(util.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        inFastq = util.file.mkstempfname('.fastq')
        outFastq = util.file.mkstempfname('.fastq')

        # Convert BAM to interleaved FASTQ
        self.samtools.bam2fq(inBam, inFastq)
        self.assertGreater(os.path.getsize(inFastq), 0)

        # Count input reads
        with open(inFastq, 'rt') as f:
            input_lines = sum(1 for _ in f)
        input_reads = input_lines // 4

        # Run bbnorm (auto-detects interleaved format)
        # Use low memory and single thread for CI compatibility
        with util.file.tmp_dir('_bbnorm_test') as tmpdir:
            self.bbmap.bbnorm(inFastq, outFastq, tmpdir=tmpdir, target=50,
                              threads=1, memory='250m')

        # Verify output exists and has reads
        self.assertTrue(os.path.isfile(outFastq))
        self.assertGreater(os.path.getsize(outFastq), 0)

        # Count output reads
        with open(outFastq, 'rt') as f:
            output_lines = sum(1 for _ in f)
        output_reads = output_lines // 4

        # Output should have reads and be <= input
        self.assertGreater(output_reads, 0)
        self.assertLessEqual(output_reads, input_reads)

    def test_bbnorm_single_end(self):
        """Test bbnorm on single-end FASTQ file."""
        # Create single-end FASTQ from TestPerSample/in.2libs3rgs.bam
        inBam = os.path.join(util.file.get_test_input_path(), 'TestPerSample', 'in.2libs3rgs.bam')
        inFastq = util.file.mkstempfname('.fastq')
        outFastq = util.file.mkstempfname('.fastq')

        # Convert BAM to FASTQ (single-end)
        self.samtools.bam2fq(inBam, inFastq)
        self.assertGreater(os.path.getsize(inFastq), 0)

        # Count input reads
        with open(inFastq, 'rt') as f:
            input_lines = sum(1 for _ in f)
        input_reads = input_lines // 4

        # Run bbnorm (auto-detects single-end format)
        # Use low memory and single thread for CI compatibility
        with util.file.tmp_dir('_bbnorm_test') as tmpdir:
            self.bbmap.bbnorm(inFastq, outFastq, tmpdir=tmpdir, target=50,
                              threads=1, memory='250m')

        # Verify output exists and has reads
        self.assertTrue(os.path.isfile(outFastq))
        self.assertGreater(os.path.getsize(outFastq), 0)

        # Count output reads
        with open(outFastq, 'rt') as f:
            output_lines = sum(1 for _ in f)
        output_reads = output_lines // 4

        # Output should have reads and be <= input
        self.assertGreater(output_reads, 0)
        self.assertLessEqual(output_reads, input_reads)
