# Unit tests for tools.samtools

__author__ = "dpark@broadinstitute.org"

import unittest
import os, os.path
import tempfile
import shutil
import Bio.SeqIO, Bio.SeqRecord, Bio.Seq
import util
import util.file
import tools
import tools.samtools
from test import TestCaseWithTmp


class TestToolSamtools(TestCaseWithTmp):

    def test_count_bam(self):
        sam = os.path.join(util.file.get_test_input_path(self), 'simple.sam')
        n = tools.samtools.SamtoolsTool().count(sam, ['-S'])
        self.assertEqual(n, 2)

    def test_fasta_index(self):
        orig_ref = os.path.join(util.file.get_test_input_path(self), 'in.fasta')
        expected_fai = os.path.join(util.file.get_test_input_path(self), 'in.fasta.fai')
        samtools = tools.samtools.SamtoolsTool()
        for ext in ('.fasta', '.fa'):
            inRef = util.file.mkstempfname(ext)
            shutil.copyfile(orig_ref, inRef)
            outFai = inRef + '.fai'
            samtools.faidx(inRef)
            self.assertEqualContents(outFai, expected_fai)

    def test_messy_fasta_index(self):
        orig_ref = os.path.join(util.file.get_test_input_path(), 'TestToolPicard', 'messy-headers.fasta')
        samtools = tools.samtools.SamtoolsTool()
        with util.file.tempfname('.fasta') as inRef:
            shutil.copyfile(orig_ref, inRef)
            samtools.faidx(inRef, overwrite=True)
        with open(inRef + '.fai', 'rt') as inf:
            seqnames = set()
            for line in inf:
                seq_name = line.strip().split('\t')[0]
                # old versions of code cut this off at "Influenza"
                self.assertGreater(len(seq_name), 50)
                seqnames.add(seq_name)
            # require that all sequence names are unique
            self.assertEqual(len(seqnames), 8)

    def test_isEmpty(self):
        samtools = tools.samtools.SamtoolsTool()
        self.assertTrue(samtools.isEmpty(os.path.join(util.file.get_test_input_path(), 'empty.bam')))
        self.assertFalse(samtools.isEmpty(os.path.join(util.file.get_test_input_path(), 'almost-empty.bam')))
        self.assertFalse(samtools.isEmpty(os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')))
        self.assertFalse(samtools.isEmpty(os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')))

    def test_sam_downsample(self):
        desired_count = 100
        tolerance = 0.1

        in_sam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        out_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()

        samtools.downsample_to_approx_count(in_sam, out_bam, desired_count)

        assert samtools.count(out_bam) in range(
            int(desired_count - (desired_count * tolerance)), int(desired_count + (desired_count * tolerance))+1
        ), "Downsampled bam file does not contain the expected number of reads within tolerance: %s" % tolerance

    def test_filterByCigarString(self):
        # The test input contains three reads to remove; one each: 
        #   leading indel, trailing indel, both leading and trailing
        # It also has a cigar string with an indel between alignment matches
        in_sam = os.path.join(util.file.get_test_input_path(self), 'indel_cigar.sam')
        out_bam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()

        # We'll use the default regex, which matches leading or trailing indels.
        # It is reproduced here in case the default changes:
        # '^((?:[0-9]+[ID]){1}(?:[0-9]+[MNIDSHPX=])+)|((?:[0-9]+[MNIDSHPX=])+(?:[0-9]+[ID]){1})$'
        samtools.filterByCigarString(in_sam, out_bam)

        assert samtools.count(out_bam)==39, "Output read count does not match the expected count."

    def test_bam2fa(self):
        samtools = tools.samtools.SamtoolsTool()
        sam = os.path.join(util.file.get_test_input_path(self), 'simple.sam')

        with samtools.bam2fa_tmp(sam) as (fa1, fa2):
            for fa in (fa1, fa2):
                assert len(list(Bio.SeqIO.parse(fa, 'fasta')))==1

        assert not os.path.isfile(fa1) and not os.path.isfile(fa2)


class TestSamtoolsImport(TestCaseWithTmp):
    """Tests for samtools import (FASTQ to BAM conversion)"""

    def test_import_paired_fastq_basic(self):
        """Test basic paired FASTQ to BAM conversion with sample name only"""
        # Use existing test FASTQ files from TestFastqBam
        inFastq1 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in1.fastq')
        inFastq2 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in2.fastq')
        outBam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        samtools.import_fastq(inFastq1, inFastq2, outBam, sample_name='TestSample')

        # Verify BAM file was created and is non-empty
        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

        # Verify read count matches input (1 read pair = 2 reads)
        self.assertEqual(samtools.count(outBam), 2)

        # Verify BAM is not flagged as empty
        self.assertFalse(samtools.isEmpty(outBam))

    def test_import_with_full_read_group(self):
        """Test that all RG tags are correctly set in BAM header"""
        import pysam

        inFastq1 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in1.fastq')
        inFastq2 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in2.fastq')
        outBam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        samtools.import_fastq(
            inFastq1, inFastq2, outBam,
            sample_name='FreeSample',
            library_name='Alexandria',
            platform='illumina',
            read_group_id='MyRG',
            platform_unit='FC001.1.ATCG',
            sequencing_center='BroadInstitute',
            run_date='2024-12-14'
        )

        # Read BAM header and verify RG tags
        with pysam.AlignmentFile(outBam, 'rb', check_sq=False) as bam:
            rg_list = bam.header.get('RG', [])
            self.assertEqual(len(rg_list), 1)
            rg = rg_list[0]

            self.assertEqual(rg.get('ID'), 'MyRG')
            self.assertEqual(rg.get('SM'), 'FreeSample')
            self.assertEqual(rg.get('LB'), 'Alexandria')
            self.assertEqual(rg.get('PL'), 'illumina')
            self.assertEqual(rg.get('PU'), 'FC001.1.ATCG')
            self.assertEqual(rg.get('CN'), 'BroadInstitute')
            self.assertEqual(rg.get('DT'), '2024-12-14')

    def test_import_rg_defaults(self):
        """Test that RG defaults are applied correctly when optional params omitted"""
        import pysam

        inFastq1 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in1.fastq')
        inFastq2 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in2.fastq')
        outBam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        samtools.import_fastq(inFastq1, inFastq2, outBam, sample_name='TestSample')

        with pysam.AlignmentFile(outBam, 'rb', check_sq=False) as bam:
            rg_list = bam.header.get('RG', [])
            self.assertEqual(len(rg_list), 1)
            rg = rg_list[0]

            # ID should default to 'A' (matching Picard default)
            self.assertEqual(rg.get('ID'), 'A')
            self.assertEqual(rg.get('SM'), 'TestSample')
            # LB should default to sample_name if not specified
            self.assertEqual(rg.get('LB'), 'TestSample')
            # PL should default to 'illumina'
            self.assertEqual(rg.get('PL'), 'illumina')
            # Optional tags should not be present
            self.assertIsNone(rg.get('PU'))
            self.assertIsNone(rg.get('CN'))
            self.assertIsNone(rg.get('DT'))

    def test_import_read_flags(self):
        """Test that paired-end read flags are set correctly (77/141 for unmapped pairs)"""
        import pysam

        inFastq1 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in1.fastq')
        inFastq2 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in2.fastq')
        outBam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        samtools.import_fastq(inFastq1, inFastq2, outBam, sample_name='TestSample')

        with pysam.AlignmentFile(outBam, 'rb', check_sq=False) as bam:
            reads = list(bam.fetch(until_eof=True))
            self.assertEqual(len(reads), 2)

            # First read should have flag 77 (paired, unmapped, mate unmapped, first in pair)
            # Second read should have flag 141 (paired, unmapped, mate unmapped, second in pair)
            flags = sorted([r.flag for r in reads])
            self.assertEqual(flags, [77, 141])

    def test_import_read_rg_tag(self):
        """Test that reads have the RG tag set correctly"""
        import pysam

        inFastq1 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in1.fastq')
        inFastq2 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in2.fastq')
        outBam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        samtools.import_fastq(
            inFastq1, inFastq2, outBam,
            sample_name='TestSample',
            read_group_id='CustomRG'
        )

        with pysam.AlignmentFile(outBam, 'rb', check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                self.assertTrue(read.has_tag('RG'))
                self.assertEqual(read.get_tag('RG'), 'CustomRG')

    def test_import_empty_fastq(self):
        """Test handling of empty FASTQ inputs - should create empty BAM with header"""
        import pysam

        # Create empty FASTQ files
        emptyFastq1 = util.file.mkstempfname('.fastq')
        emptyFastq2 = util.file.mkstempfname('.fastq')
        outBam = util.file.mkstempfname('.bam')

        open(emptyFastq1, 'w').close()
        open(emptyFastq2, 'w').close()

        samtools = tools.samtools.SamtoolsTool()
        samtools.import_fastq(
            emptyFastq1, emptyFastq2, outBam,
            sample_name='EmptySample',
            library_name='EmptyLib'
        )

        # Verify BAM was created
        self.assertTrue(os.path.exists(outBam))
        self.assertGreater(os.path.getsize(outBam), 0)

        # Verify BAM is empty (no reads)
        self.assertTrue(samtools.isEmpty(outBam))
        self.assertEqual(samtools.count(outBam), 0)

        # Verify header has RG
        with pysam.AlignmentFile(outBam, 'rb', check_sq=False) as bam:
            rg_list = bam.header.get('RG', [])
            self.assertEqual(len(rg_list), 1)
            self.assertEqual(rg_list[0].get('SM'), 'EmptySample')
            self.assertEqual(rg_list[0].get('LB'), 'EmptyLib')

    def test_import_multithreaded(self):
        """Test that threads parameter works without error"""
        inFastq1 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in1.fastq')
        inFastq2 = os.path.join(util.file.get_test_input_path(), 'TestFastqBam', 'in2.fastq')
        outBam = util.file.mkstempfname('.bam')

        samtools = tools.samtools.SamtoolsTool()
        # Should work with explicit thread count
        samtools.import_fastq(
            inFastq1, inFastq2, outBam,
            sample_name='TestSample',
            threads=2
        )

        self.assertTrue(os.path.exists(outBam))
        self.assertEqual(samtools.count(outBam), 2)

