'''utilities for tests'''

__author__ = "irwin@broadinstitute.org"

def assert_equal_contents(testCase, file1, file2) :
    'Assert contents of two files are equal for a unittest.TestCase'
    testCase.assertEqual(open(file1).read(), open(file2).read())