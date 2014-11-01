'''utilities for tests'''

__author__ = "irwin@broadinstitute.org"

import filecmp

def assert_equal_contents(testCase, filename1, filename2) :
    'Assert contents of two files are equal for a unittest.TestCase'
    testCase.assertTrue(filecmp.cmp(filename1, filename2, shallow=False))
