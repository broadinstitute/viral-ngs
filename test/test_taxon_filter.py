#!/usr/bin/env python
# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org"

import unittest, os, sys, tempfile
# The following line is needed to access taxon_filter and util when running from shell
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import taxon_filter, util, tools.last
# The following line is needed to access test_consensus in python 3.4 
sys.path.append(os.path.dirname(__file__))
from test_consensus import set_tmpDir, destroy_tmpDir

class TestCommandHelp(unittest.TestCase):
	def test_help_parser_for_each_command(self):
		for cmd_name, main_fun, parser_fun in taxon_filter.__commands__:
			parser = parser_fun()
			helpstring = parser.format_help()


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

class TestFilterLastal(unittest.TestCase):
	def setUp(self):
		set_tmpDir('TestFilterLastal')
	def tearDown(self):
		destroy_tmpDir()
	def test_filter_lastal(self) :
		# Create refDbs
		refFasta = os.path.join(os.path.dirname(__file__), 'input', 'ebola.fasta')
		dbsDir = tempfile.mkdtemp()
		refDbs = os.path.join(dbsDir, 'ebola')
		lastdbPath = tools.last.Lastdb().install_and_get_path()
		os.system('{lastdbPath} {refDbs} {refFasta}'.format(lastdbPath = lastdbPath,
															refDbs = refDbs,
															refFasta = refFasta))
		inFastq = util.file.mkstempfname()
		outFastq = util.file.mkstempfname()
		open(inFastq, 'w').write(filterLastalInput)
		args = taxon_filter.parser_filter_lastal().parse_args([inFastq, refDbs, outFastq])
		taxon_filter.main_filter_lastal(args)
		self.assertEqual(open(outFastq + '.fastq').read(), filterLastalExpected)

if __name__ == '__main__':
    unittest.main()
