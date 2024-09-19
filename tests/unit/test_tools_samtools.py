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


