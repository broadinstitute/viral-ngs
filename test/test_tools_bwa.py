#!/usr/bin/env python
# Unit tests for taxon_filter.py

__author__ = "hlevitin@broadinstitute.org"

import unittest, os, sys, tempfile
import util.file, tools.bwa

class TestToolBwa(unittest.TestCase) :

	def setUp(self) :
		util.file.set_tmpDir('TestToolBwa')

	def tearDown(self) :
		util.file.destroy_tmpDir()

	def test_tool_bwa_index(self) :
		referenceDir = util.file.get_test_input_path()
		expectedDir = util.file.get_test_input_path(self)

		fasta = os.path.join(referenceDir, 'ebola.fasta')
		bwa = tools.bwa.Bwa()
		bwa.install()

		bwa.execute('index', [fasta])

		expected_fasta = os.path.join(expectedDir, 'ebola_expected.fasta')
		extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
		for ext in extensions:
			result = '{}.{}'.format(fasta, ext)
			expect = '{}.{}'.format(expected_fasta, ext)
			self.assertEqual(open(result).read(),
							open(expect).read())

if __name__ == '__main__':
    unittest.main()
