'''utilities for tests'''

__author__ = "irwin@broadinstitute.org"

# built-ins
import filecmp
import os
import re
import unittest
import hashlib
import logging

# third-party
import Bio.SeqIO
import pytest

# intra-project
import util.file
from util.misc import available_cpu_count
from tools.samtools import SamtoolsTool

logging.getLogger('botocore').setLevel(logging.WARNING)
logging.getLogger('boto3').setLevel(logging.WARNING)


if 'PYTEST_XDIST_WORKER_COUNT' in os.environ:
    _CPUS = 1
else:
    _CPUS = available_cpu_count()


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

    samtools = SamtoolsTool()

    sam_one = util.file.mkstempfname(".sam")
    sam_two = util.file.mkstempfname(".sam")

    # write the bam files to sam format, without header (no -h)
    samtools.view(args=[], inFile=bam_filename1, outFile=sam_one)
    samtools.view(args=[], inFile=bam_filename2, outFile=sam_two)

    try:
        testCase.assertTrue(filecmp.cmp(sam_one, sam_two, shallow=False))
    finally:
        for fname in [sam_one, sam_two]:
            if os.path.exists(fname):
                os.unlink(fname)

def assert_md5_equal_to_line_in_file(testCase, filename, checksum_file, msg=None):
    ''' Compare the checksum of a test file with the expected checksum
        stored in a second file
          compare md5(test_output.bam) with expected_output.bam.md5:1
    '''

    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    expected_checksum = ""
    with open(checksum_file, "rb") as f:
        expected_checksum = str(f.readline().decode("utf-8"))

    expected_checksum = expected_checksum.replace("\r","").replace("\n","")

    assert len(expected_checksum) > 0

    testCase.assertEqual(hash_md5.hexdigest(), expected_checksum, msg=msg)

@pytest.mark.usefixtures('tmpdir_class')
class TestCaseWithTmp(unittest.TestCase):
    'Base class for tests that use tempDir'

    def assertEqualContents(self, f1, f2):
        assert_equal_contents(self, f1, f2)

    def assertEqualFasta(self, f1, f2):
        '''Check that two fasta files have the same sequence ids and bases'''
        def seqIdPairs(f):
            return [(rec.id, rec.seq) for rec in Bio.SeqIO.parse(f, 'fasta')]
        self.assertEqual(seqIdPairs(f1), seqIdPairs(f2))

    def assertEqualFastaSeqs(self, f1, f2):
        '''Check that two fasta files have the same sequence bases (sequence ids may differ)'''
        def seqs(f):
            return [rec.seq for rec in Bio.SeqIO.parse(f, 'fasta')]
        self.assertEqual(seqs(f1), seqs(f2))

    def input(self, fname):
        '''Return the full filename for a file in the test input directory for this test class'''
        return os.path.join(util.file.get_test_input_path(self), fname)

    def inputs(self, *fnames):
        '''Return the full filenames for files in the test input directory for this test class'''
        return [self.input(fname) for fname in fnames]

"""
When "nose" executes python scripts for automated testing, it excludes ones with
the executable bit set (in case they aren't import safe). To prevent any of the
tests in this folder from being silently excluded, assure this bit is not set.
"""


def assert_none_executable():
    testDir = os.path.dirname(__file__)
    assert all(not os.access(os.path.join(testDir, filename), os.X_OK) for filename in os.listdir(testDir)
               if filename.endswith('.py'))
