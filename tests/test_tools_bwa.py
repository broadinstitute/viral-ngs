# Unit tests for bwa tool

__author__ = "hlevitin@broadinstitute.org"

import unittest, os, sys, tempfile
import util.file, tools.bwa

class TestToolBwa(unittest.TestCase) :

    def setUp(self) :
        util.file.set_tmpDir('TestToolBwa')
        self.bwa = tools.bwa.Bwa()
        self.bwa.install()


    def tearDown(self) :
        util.file.destroy_tmpDir()

    def test_tool_bwa_index(self) :
        referenceDir = util.file.get_test_input_path()
        expectedDir = util.file.get_test_input_path(self)

        fasta = os.path.join(referenceDir, 'ebola.fasta')
        self.bwa.execute('index', [fasta])

        expected_fasta = os.path.join(expectedDir, 'ebola_expected.fasta')
        extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
        for ext in extensions:
            result = open('{}.{}'.format(fasta, ext), 'rb')
            expect = open('{}.{}'.format(expected_fasta, ext), 'rb')

            self.assertEqual(result.read(),
                             expect.read() )

            result.close(); expect.close()

    def test_tool_bwa_aln(self) :
        # referenceDir = util.file.get_test_input_path()
        expectedDir = util.file.get_test_input_path(self)

        # can used expected out for index as input
        reference = os.path.join(expectedDir, 'ebola_expected.fasta')
        fastq = os.path.join(expectedDir, 'ebola_aln_input.fastq')
        output = "{}.sai".format(util.file.mkstempfname())
        expect = os.path.join(expectedDir, 'ebola_aln_expected.sai')


        bwa = tools.bwa.Bwa()
        bwa.install()

        bwa.execute('aln', [reference, fastq], options={'-q': 5, '-t': 4},
                post_cmd=" > {}".format(output))

        result = open(output, 'rb')
        gold = open(expect, 'rb')
        self.assertEqual(result.read(),
                         gold.read())
        result.close(); gold.close()

if __name__ == '__main__':
    unittest.main()
