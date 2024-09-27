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
        in_ref = util.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, in_ref)
        reads = os.path.join(util.file.get_test_input_path(self), 'ebov_reads.bam')
        out_bam = util.file.mkstempfname('.bam')
        self.bbmap.align(in_bam=reads, ref_fasta=in_ref, out_bam=out_bam)
        self.assertTrue(os.path.isfile(out_bam))
        self.assertTrue(os.path.getsize(out_bam))

    def test_dedup_clumpify(self):
        reads = os.path.join(util.file.get_test_input_path(self), 'ebov_reads.bam')
        expected_bam = os.path.join(util.file.get_test_input_path(self), 'ebov_reads_clumpify_dedup_expected.bam')
        out_bam = util.file.mkstempfname('.bam')
        self.bbmap.dedup_clumpify(in_bam=reads, out_bam=out_bam)

        target_count = self.samtools.count(expected_bam)

        self.assertTrue(os.path.isfile(out_bam))
        self.assertTrue(os.path.getsize(out_bam))
        # check that the target count is within 1% of the expected count
        self.assertAlmostEqual(self.samtools.count(out_bam), target_count, delta=target_count*0.01, msg="{} not deduplicated to the target size: {}".format(os.path.basename(out_bam),target_count))
