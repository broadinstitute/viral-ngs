# Unit tests for MEGAN

#from builtins import super
import unittest
import tools.megan
from test import TestCaseWithTmp


class TestToolMegan(TestCaseWithTmp):

    def setUp(self):
        #super().setUp()
        super(TestToolMegan, self).setUp()
        self.megan = tools.megan.Megan()
        self.megan.install()

    # Just a simple test that MEGAN can open/close.
    @unittest.skip('Requires installed MEGAN')
    def test_megan_runs(self):
        self.assertEqual(0, self.megan.execute('quit;\n')[0])


if __name__ == '__main__':
    unittest.main()
