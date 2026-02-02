'''utilities for tests'''

__author__ = "irwin@broadinstitute.org"

# built-ins
import filecmp
import os
import re
import unittest
import hashlib
import logging
import copy

# third-party
import Bio.SeqIO
import pytest
import pysam

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
        if not filecmp.cmp(sam_one, sam_two, shallow=False):
            # Files are different - show debug output
            print("\n" + "="*80)
            print(f"DEBUG: BAM files are not equal")
            print(f"  File 1: {bam_filename1}")
            print(f"  File 2: {bam_filename2}")
            print("\nFirst 5 reads from each file:")
            print("-" * 80)
            print("File 1:")
            with open(sam_one, 'r') as f:
                for i, line in enumerate(f):
                    if i >= 5:
                        break
                    print(f"  {line.rstrip()}")
            print("-" * 80)
            print("File 2:")
            with open(sam_two, 'r') as f:
                for i, line in enumerate(f):
                    if i >= 5:
                        break
                    print(f"  {line.rstrip()}")
            print("="*80 + "\n")
            testCase.fail("BAM files are not equal (see debug output above)")
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

    def assertEqualSamHeaders(self, tested_samfile, expected_samfile, other_allowed_values=None):
        '''
            other_allowed_values is a dict that maps a key to a list of other accepted values
        '''
        other_allowed_values = other_allowed_values or {}

        test_sam = pysam.AlignmentFile(tested_samfile, "rb", check_sq=False)
        expected_sam = pysam.AlignmentFile(expected_samfile, "rb", check_sq=False)

        # check that the two sams contain the same types of header lines
        # note that pysam returns lists of dicts for
        # header lines types that appear more than once (ex. multiple 'RG' lines, etc.)
        self.assertEqual(sorted(test_sam.header.keys()),sorted(expected_sam.header.keys()))

        lines_tag_to_check = test_sam.header.keys()

        for line_tag in lines_tag_to_check:
            line_type = type(test_sam.header[line_tag])
            # for header lines that occur more than once
            if line_type == list: # if list of dicts
                # compose alterates for each line present, and check to see if exists in list of expected
                for test_dict in test_sam.header[line_tag]:
                    alternate_allowed_test_dicts = []
                    for key,val in test_dict.items():
                        alterate_test_dict = copy.deepcopy(test_dict)
                        if key in other_allowed_values.keys():
                            for alternate_val in other_allowed_values[key]:
                                alterate_test_dict[key] = alternate_val
                                alternate_allowed_test_dicts.append(copy.deepcopy(alterate_test_dict))
                    for expected_d in expected_sam.header[line_tag]:
                        self.assertTrue(any(d==expected_d for d in alternate_allowed_test_dicts+[test_dict]), msg="{} expected but not seen in {}".format(expected_d, alternate_allowed_test_dicts+[test_dict]))
            # for header lines that occur only once
            elif line_type in (dict,): # may need to change object type for pysam >=0.18
                test_dict = dict(test_sam.header[line_tag])
                alternate_allowed_test_dicts = []
                for key,val in test_dict.items():
                    alterate_test_dict = copy.deepcopy(test_dict)
                    if key in other_allowed_values.keys():
                        for alternate_val in other_allowed_values[key]:
                            alterate_test_dict[key] = alternate_val
                            alternate_allowed_test_dicts.append(copy.deepcopy(alterate_test_dict))
                self.assertTrue(any(d==expected_sam.header[line_tag] for d in alternate_allowed_test_dicts+[test_dict]), msg="{} expected but not seen in {}".format(expected_sam.header[line_tag],alternate_allowed_test_dicts+[test_dict]))

"""
When "nose" executes python scripts for automated testing, it excludes ones with
the executable bit set (in case they aren't import safe). To prevent any of the
tests in this folder from being silently excluded, assure this bit is not set.
"""


def assert_none_executable():
    testDir = os.path.dirname(__file__)
    assert all(not os.access(os.path.join(testDir, filename), os.X_OK) for filename in os.listdir(testDir)
               if filename.endswith('.py'))
