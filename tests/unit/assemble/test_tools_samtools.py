# Unit tests for tools.samtools

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import tempfile
import shutil
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

if __name__ == '__main__':
    unittest.main()
