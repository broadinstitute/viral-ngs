# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org"

import taxon_filter, util
from test_consensus import set_tmpDir, destroy_tmpDir
import unittest, os

class TestCommandHelp(unittest.TestCase):
	def test_help_parser_for_each_command(self):
		for cmd_name, main_fun, parser_fun in taxon_filter.__commands__:
			parser = parser_fun()
			helpstring = parser.format_help()

class TestTaxonFilter(unittest.TestCase):
	def testNothingAtAll(self):
		'''here we test nothing at all and this should pass'''
		pass

class TestFilterLastal(unittest.TestCase):
	def setUp(self):
		set_tmpDir('TestFilterLastal')
	def tearDown(self):
		destroy_tmpDir()
	def test_filter_lastal(self) :
		filterLastalInput = """@fakeRead
		AGTACATGCAGAGCAAGGACTGATACAATATCCAACAGCTTGGCAATCAGTAGGACACATGATGGTGA
		+
		CCCFFFFFHHHHHJJJJJJJJJJJJJJJHIIIIJJJJHIJIIJJIJJFHGIIJJGHHHBDFDDDDDDD
		"""
		filterLastalExpected = """@fakeRead
		AGTACATGCAGAGCAAGGACTGATACAATATCCAACAGCTTGGCAATCAGTAGGACACATGATGGTGA
		+fakeRead
		CCCFFFFFHHHHHJJJJJJJJJJJJJJJHIIIIJJJJHIJIIJJIJJFHGIIJJGHHHBDFDDDDDDD
		"""
		refDbs = os.path.join(os.path.dirname(__file__), 'input', 'ebola')
		inFastq = util.file.mkstempfname()
		outFastq = util.file.mkstempfname()
		open(inFastq, 'w').write(filterLastalInput)
		args = taxon_filter.parser_filter_lastal().parse_args([inFastq, refDbs, outFastq])
		taxon_filter.main_filter_lastal(args)
		self.assertEqual(open(outFastq + '.fastq').read(), filterLastalExpected)

if __name__ == '__main__':
    unittest.main()
