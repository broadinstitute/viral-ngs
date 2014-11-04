'''utilities for tests'''

__author__ = "irwin@broadinstitute.org"

import filecmp, os

"""
When "nose" executes python scripts for automated testing, it excludes ones with
the executable bit set (in case they aren't import safe). To prevent any of the 
tests in this folder from being silently excluded, assure this bit is not set.
"""
def assert_none_executable() :
    testDir = os.path.dirname(__file__)
    assert all(not os.access(os.path.join(testDir, filename), os.X_OK)
               for filename in os.listdir(testDir)
               if filename.endswith('.py'))
assert_none_executable()

def assert_equal_contents(testCase, filename1, filename2) :
    'Assert contents of two files are equal for a unittest.TestCase'
    testCase.assertTrue(filecmp.cmp(filename1, filename2, shallow=False))
