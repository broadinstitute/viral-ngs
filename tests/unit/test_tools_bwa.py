# Unit tests for bwa tool

__author__ = "hlevitin@broadinstitute.org"

import unittest, os.path, shutil
import util.file, tools.bwa
from test import TestCaseWithTmp

class TestToolBwa(TestCaseWithTmp) :

    def setUp(self) :
        super(TestToolBwa, self).setUp()
        self.bwa = tools.bwa.Bwa()
        self.bwa.install()

    def test_index(self) :
        orig_ref = os.path.join(util.file.get_test_input_path(),
            'ebola.fasta')
        inRef = util.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)

        expected_fasta = os.path.join(
            util.file.get_test_input_path(self),
            'ebola_expected.fasta')

        self.bwa.execute('index', [inRef])
        for ext in ('amb', 'ann', 'bwt', 'pac', 'sa'):
            self.assertEqualContents(inRef+'.'+ext, expected_fasta+'.'+ext)

    def test_aln(self) :
        expectedDir = util.file.get_test_input_path(self)

        # can used expected out for index as input
        reference = os.path.join(expectedDir, 'ebola_expected.fasta')
        fastq = os.path.join(expectedDir, 'ebola_aln_input.fastq')
        output = util.file.mkstempfname('.sai')
        expect = os.path.join(expectedDir, 'ebola_aln_expected.sai')

        self.bwa.execute('aln', [reference, fastq],
            options={'-q': 5, '-t': 4},
            post_cmd=" > {}".format(output))
        
        self.assertEqualContents(output, expect)


if __name__ == '__main__':
    unittest.main()
