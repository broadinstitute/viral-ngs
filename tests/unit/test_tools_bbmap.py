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
        self.bbmap.align(inBam=reads, refFasta=inRef, outBam=outBam)
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))
