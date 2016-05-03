'''utilities for tests'''

__author__ = "irwin@broadinstitute.org"

# built-ins
import filecmp
import os
import unittest

# third-party
import pysam

# intra-project
import util.file

def assert_equal_contents(testCase, filename1, filename2):
    'Assert contents of two files are equal for a unittest.TestCase'
    testCase.assertTrue(filecmp.cmp(filename1, filename2, shallow=False))

def assert_equal_bam_reads(testCase, bam_filename1, bam_filename2):
    ''' Assert that two bam files are equivalent

        This function converts each file to sam (plaintext) format,
        without header, since the header can be variable.

        Test data should be stored in bam format rather than sam
        to save space, and to ensure the bam->sam conversion
        is the same for both files.
    '''

    sam_one = util.file.mkstempfname(".sam")
    sam_two = util.file.mkstempfname(".sam")

    sorted_sam_one = util.file.mkstempfname("sorted.sam")
    sorted_sam_two = util.file.mkstempfname("sorted.sam")

    # write the bam files to sam format, without header (no -h)
    pysam.view('-o', sam_one, bam_filename1, catch_stdout=False)
    pysam.view('-o', sam_two, bam_filename2, catch_stdout=False)

    try:
        testCase.assertTrue(filecmp.cmp(sorted_sam_one, sorted_sam_two, shallow=False))
    finally:
        for fname in [sam_one, sam_two, sorted_sam_one, sorted_sam_two]:
            if os.path.exists(fname):
                os.unlink(fname)

class TestCaseWithTmp(unittest.TestCase):
    'Base class for tests that use tempDir'

    @classmethod
    def setUpClass(cls):
        cls._class_tempdir = util.file.set_tmp_dir(cls.__name__)

    def setUp(self):
        util.file.set_tmp_dir(type(self).__name__)

    @classmethod
    def tearDownClass(cls):
        util.file.destroy_tmp_dir(cls._class_tempdir)

    def tearDown(self):
        util.file.destroy_tmp_dir()

    def assertEqualContents(self, f1, f2):
        assert_equal_contents(self, f1, f2)


"""
When "nose" executes python scripts for automated testing, it excludes ones with
the executable bit set (in case they aren't import safe). To prevent any of the
tests in this folder from being silently excluded, assure this bit is not set.
"""


def assert_none_executable():
    testDir = os.path.dirname(__file__)
    assert all(not os.access(os.path.join(testDir, filename), os.X_OK) for filename in os.listdir(testDir)
               if filename.endswith('.py'))


assert_none_executable()
