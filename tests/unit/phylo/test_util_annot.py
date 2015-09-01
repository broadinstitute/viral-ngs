# Unit tests for util.annot.py

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
