#!/usr/bin/env python
# Unit tests for taxon_filter.py

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
			result = '{}.{}'.format(fasta, ext)
			expect = '{}.{}'.format(expected_fasta, ext)
			self.assertEqual(open(result).read(),
							open(expect).read())

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

		self.assertEqual(open(output).read(),
							open(expect).read())
		# extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
		# for ext in extensions:
			# result = '{}.{}'.format(fasta, ext)
			# expect = '{}.{}'.format(expected_fasta, ext)

if __name__ == '__main__':
    unittest.main()
