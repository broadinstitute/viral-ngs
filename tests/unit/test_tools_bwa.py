# Unit tests for bwa aligner

__author__ = "tomkinsc@broadinstitute.org"

import os.path
import tempfile
import subprocess
import util.file
import tools.bwa
import tools.samtools
from test import TestCaseWithTmp

class TestToolBwa(TestCaseWithTmp):

    def setUp(self):
        super(TestToolBwa, self).setUp()
        self.bwa = tools.bwa.Bwa()
        self.bwa.install()
        self.samtools = tools.samtools.SamtoolsTool()

        self.tempDir = tempfile.mkdtemp()
        root_input_dir = util.file.get_test_input_path()
        ref_fasta = os.path.join(root_input_dir, '5kb_human_from_chr6.fasta')

        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")

        # create bwa db
        self.bwadb_path = self.bwa.index(ref_fasta, output=self.database_prefix_path)

    def test_working_bam(self):
        in_bam = os.path.join(util.file.get_test_input_path(), "TestDepleteHuman", 'test-reads-human.bam')

        outfile = util.file.mkstempfname('.bam')
        self.bwa.mem(in_bam, self.bwadb_path, outfile, options=['-a'])

        self.assertEqual(self.samtools.count(outfile), 20)

        os.unlink(outfile)

    def test_corrupt_bam(self):
        in_bam = os.path.join(util.file.get_test_input_path(), 'broken.bam')

        outfile = util.file.mkstempfname('.bam')

        # pipe.poll() should raise an exception
        self.assertRaises(subprocess.CalledProcessError, self.bwa.mem, in_bam, self.bwadb_path, outfile, options=['-a'])

        os.unlink(outfile)
