'''Unit tests for util.annot.py'''

__author__ = "dpark@broadinstitute.org"

import util.annot
import unittest

class TestGeneDb(unittest.TestCase):
	def setUp(self):
		pass
		
	def testNothingAtAll(self):
		'''here we test nothing at all and this should pass'''
		pass
		
	def testTautology(self):
		'''here we test 1 = 1'''
		self.assertEqual(1, 1)

''' some tests to write

- GeneDb.getByGeneId returns empty list for IDs not in DB
- GeneDb.getByGeneId returns always returns 0 or 1 hits in list if asked for a 'gene'
- GeneDb.getByGeneId always returns features in sorted order (chr,start) if asked for a child feature type
- GeneDb.getFeaturesInRegion properly handles edge cases
- GeneDb.getGeneNameById returns None for IDs not in DB
'''
