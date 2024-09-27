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
import math
import csv
import itertools

# third-party
import Bio.SeqIO
import pytest
import pysam

# intra-project
import util.file
from util.misc import available_cpu_count, zip_dicts, is_number
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

    def assertApproxEqualValuesInDelimitedFiles(self, file_one, file_two, dialect="tsv", numeric_rel_tol=1e-5, header_lines_to_skip=0, use_first_processed_line_for_fieldnames=False):
        '''
            This test checks whether two delimited (tsv, csv, or custom dialect) files
            are approximately the same, by comparing the values present at each line
            with a configurable tolerance for comparison of numeric values, and exact
            comparison of string values. A specified number of header lines can be skipped 
            before the comparison begins to account for lines written by some tools with
            for metadata, invokation params, etc.
            If use_first_processed_line_for_fieldnames is specified, the columns can be
            in any order since the rows are parsed as dicts, otherwise the columns are
            assumed to be in the same order in both of the files.
        '''
        header_fieldnames=None

        csv.register_dialect('tsv', quoting=csv.QUOTE_MINIMAL, delimiter="\t")
        csv.register_dialect('csv', quoting=csv.QUOTE_MINIMAL, delimiter=",")
        
        with util.file.open_or_gzopen(file_one, 'rU') as inf1, util.file.open_or_gzopen(file_two, 'rU') as inf2:
            report_type=None
            for line_num, (line1,line2) in enumerate(itertools.zip_longest(inf1,inf2)):
                self.assertIsNotNone(line1, msg="%s appears to be shorter than %s" % (inf1, inf2))
                self.assertIsNotNone(line2, msg="%s appears to be shorter than %s" % (inf2, inf1))

                # continue to this next pair of lines until we have 
                # skipped past the header
                if line_num < header_lines_to_skip:
                    continue
                else:
                    if header_fieldnames is None and use_first_processed_line_for_fieldnames:
                        inf1_row = next(csv.reader([line1.strip().rstrip('\n')], dialect=dialect))
                        inf2_row = next(csv.reader([line1.strip().rstrip('\n')], dialect=dialect))
                        self.assertEqual(inf1_row,inf2_row, msg="header lines are not the same")
                        header_fieldnames=inf1_row
                        continue

                    # if the header names are defined
                    if header_fieldnames is not None:
                        inf1_row = next(csv.DictReader([line1.strip().rstrip('\n')], dialect=dialect, fieldnames=header_fieldnames))
                        inf2_row = next(csv.DictReader([line2.strip().rstrip('\n')], dialect=dialect, fieldnames=header_fieldnames))
                    else:
                        inf1_row = next(csv.reader([line1.strip().rstrip('\n')], dialect=dialect))
                        inf2_row = next(csv.reader([line2.strip().rstrip('\n')], dialect=dialect))

                    # assume the rows have the same number of elements
                    self.assertTrue(len(inf1_row) == len(inf2_row), msg="Files have lines of different length on line %s %s %s" % (line_num, inf1_row, inf2_row))

                    if header_fieldnames is not None:
                        def _dict_values_only(dict1_in, dict2_in):
                            for key, values in zip_dicts(dict1_in,dict2_in):
                                yield values
                        items_to_compare=_dict_values_only(inf1_row,inf2_row)
                    else:
                        items_to_compare=zip(inf1_row,inf2_row)

                    for inf1_row_item, inf2_row_item in items_to_compare:
                        # we assume at the item from the same position is a number 
                        # or not a number in both files
                        self.assertTrue(is_number(inf1_row_item)==is_number(inf2_row_item))
                        
                        # if we're dealing with numbers, check that they're approximately equal
                        if is_number(inf1_row_item):
                            assert float(inf2_row_item) == pytest.approx(float(inf1_row_item), rel=numeric_rel_tol)
                        else:
                            # otherwise we're probably dealing with a string
                            self.assertEqual(inf2_row_item, inf1_row_item)


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
