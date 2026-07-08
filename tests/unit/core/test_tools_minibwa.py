# Unit tests for minibwa aligner

import os
import os.path
import shutil
import tempfile
import subprocess
import viral_ngs.core
import viral_ngs.core.minibwa
import viral_ngs.core.samtools
from tests import TestCaseWithTmp


class TestToolMinibwa(TestCaseWithTmp):

    def setUp(self):
        super(TestToolMinibwa, self).setUp()
        self.minibwa = viral_ngs.core.minibwa.Minibwa()
        self.samtools = viral_ngs.core.samtools.SamtoolsTool()
        self.tempDir = tempfile.mkdtemp()

    def _index(self, ref_fasta):
        # minibwa indexes in-place next to the reference; copy to a writable temp location first
        ref_copy = viral_ngs.core.file.mkstempfname('.ref.fasta')
        shutil.copyfile(ref_fasta, ref_copy)
        self.minibwa.index(ref_copy)
        return ref_copy

    def test_ebola_bam(self):
        in_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.testreads.bam')
        ref_fasta = self._index(os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta'))
        with viral_ngs.core.file.tempfname('.bam') as outfile:
            self.minibwa.align_bam(in_bam, ref_fasta, outfile)
            # illumina reads should align to the matching reference; exact count depends on
            # the minibwa version and should be tightened once verified against the installed binary
            self.assertGreater(self.samtools.count(outfile), 10000)

    def test_empty_input(self):
        in_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        ref_fasta = self._index(os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta'))
        with viral_ngs.core.file.tempfname('.bam') as outfile:
            self.minibwa.align_bam(in_bam, ref_fasta, outfile)
            # empty input should produce a valid, header-only BAM with zero reads
            self.assertEqual(self.samtools.count(outfile), 0)

    def test_corrupt_bam(self):
        # the samtools bam2fq pipe should fail and raise
        in_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'broken.bam')
        ref_fasta = self._index(os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta'))
        with viral_ngs.core.file.tempfname('.bam') as outfile:
            self.assertRaises(subprocess.CalledProcessError, self.minibwa.align_bam, in_bam, ref_fasta, outfile)
