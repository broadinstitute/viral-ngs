# Unit tests for bwa aligner

__author__ = "dpark@broadinstitute.org"

import os.path
import tempfile
import subprocess
import util.file
import tools.minimap2
import tools.samtools
from test import TestCaseWithTmp

class TestToolMinimap2(TestCaseWithTmp):

    def setUp(self):
        super(TestToolMinimap2, self).setUp()
        self.mm2 = tools.minimap2.Minimap2()
        self.samtools = tools.samtools.SamtoolsTool()
        self.tempDir = tempfile.mkdtemp()
        root_input_dir = util.file.get_test_input_path()
        self.ref_fasta = os.path.join(root_input_dir, '5kb_human_from_chr6.fasta')

    def test_human_bam(self):
        in_bam = os.path.join(util.file.get_test_input_path(), "TestDepleteHuman", 'test-reads-human.bam')
        with util.file.tempfname('.bam') as outfile:
            self.mm2.align_bam(in_bam, self.ref_fasta, outfile, options=['-x', 'sr'])
            self.assertEqual(self.samtools.count(outfile), 20)

    def test_ebola_bam(self):
        in_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        ref_fasta = os.path.join(util.file.get_test_input_path(), 'ebov-makona.fasta')
        with util.file.tempfname('.bam') as outfile:
            self.mm2.align_bam(in_bam, ref_fasta, outfile)
            self.assertEqual(self.samtools.count(outfile), 18843)

    def test_corrupt_bam(self):
        # pipe.poll() should raise an exception
        in_bam = os.path.join(util.file.get_test_input_path(), 'broken.bam')
        with util.file.tempfname('.bam') as outfile:
            self.assertRaises(subprocess.CalledProcessError, self.mm2.align_bam, in_bam, self.ref_fasta, outfile)
