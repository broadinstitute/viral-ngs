#!/usr/bin/env python
# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org"

import unittest, os, sys, tempfile
import taxon_filter, util.file, tools.last

class TestCommandHelp(unittest.TestCase):
	def test_help_parser_for_each_command(self):
		for cmd_name, main_fun, parser_fun in taxon_filter.__commands__:
			parser = parser_fun()
			helpstring = parser.format_help()

class TestTrimmomatic(unittest.TestCase) :
	def setUp(self) :
		util.file.set_tmpDir('TestTrimmomatic')
	def tearDown(self) :
		util.file.destroy_tmpDir()
	def test_trimmomatic(self) :
		inputDir = util.file.get_test_input_path(self)
		inFastq1 = os.path.join(inputDir, 'in1.fastq')
		inFastq2 = os.path.join(inputDir, 'in2.fastq')
		pairedOutFastq1 = util.file.mkstempfname()
		pairedOutFastq2 = util.file.mkstempfname()
		clipFasta = os.path.join(inputDir, 'clip.fasta')
		parser = taxon_filter.parser_trim_trimmomatic()
		args = parser.parse_args([inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2,
								 clipFasta])
		taxon_filter.main_trim_trimmomatic(args)

		# Check that results match expected
		expected1Fastq = os.path.join(inputDir, 'expected1.fastq')
		expected2Fastq = os.path.join(inputDir, 'expected2.fastq')
		self.assertEqual(open(pairedOutFastq1).read(), open(expected1Fastq).read())
		self.assertEqual(open(pairedOutFastq2).read(), open(expected2Fastq).read())

class TestFilterLastal(unittest.TestCase) :
	def setUp(self) :
		util.file.set_tmpDir('TestFilterLastal')
	def tearDown(self) :
		util.file.destroy_tmpDir()
	def test_filter_lastal(self) :
		# Create refDbs
		inputDir = util.file.get_test_input_path(self)
		refFasta = os.path.join(inputDir, 'ebola.fasta')
		dbsDir = tempfile.mkdtemp()
		refDbs = os.path.join(dbsDir, 'ebola')
		lastdbPath = tools.last.Lastdb().install_and_get_path()
		os.system('{lastdbPath} {refDbs} {refFasta}'.format(lastdbPath = lastdbPath,
															refDbs = refDbs,
															refFasta = refFasta))
		# Call main_filter_lastal
		inFastq = os.path.join(inputDir, 'in.fastq')
		outFastq = util.file.mkstempfname()
		args = taxon_filter.parser_filter_lastal().parse_args([inFastq, refDbs, outFastq])
		taxon_filter.main_filter_lastal(args)

		# Check that results match expected
		expectedFastq = os.path.join(inputDir, 'expected.fastq')
		self.assertEqual(open(outFastq + '.fastq').read(), open(expectedFastq).read())

if __name__ == '__main__':
    unittest.main()
