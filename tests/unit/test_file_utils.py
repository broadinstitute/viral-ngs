# Unit tests for ncbi.py

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import os
import tempfile
import unittest
import argparse

# module-specific
import file_utils
import util.file
from test import TestCaseWithTmp

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_func in file_utils.__commands__:
            parser = parser_func(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestTarballMerger(TestCaseWithTmp):
    def setUp(self):
        super(TestTarballMerger, self).setUp()
        self.input_dir = util.file.get_test_input_path(self)
        self.raw_files = ["file{}".format(x) for x in range(1,5)]
        self.input_mixed_files = list(os.path.join(self.input_dir, "mixed-compressed-input", fn)
                for fn in sorted(os.listdir(os.path.join(self.input_dir, "mixed-compressed-input"))))
    
    def test_simple_merge(self):
        """
            Simple execution
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        file_utils.merge_tarballs( out_tarball_file,
                                    self.input_mixed_files
                                    )

class TestTsvJoin(TestCaseWithTmp):
   
    def test_join(self):
        infiles = list(os.path.join(util.file.get_test_input_path(self), x) for x in ('tab-1.txt', 'tab-2.txt'))
        expected = os.path.join(util.file.get_test_input_path(self), 'expected-out.txt')
        actual = util.file.mkstempfname('.joined.txt')
        file_utils.tsv_join(infiles, actual, join_id='h1')
        self.assertEqualContents(expected, actual)
