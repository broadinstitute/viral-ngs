# Unit tests for ncbi.py

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import os
import tempfile
import unittest
import argparse
import tarfile

# module-specific
import file_utils
import util.file
from test import assert_equal_bam_reads, TestCaseWithTmp, assert_equal_contents, assert_md5_equal_to_line_in_file


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_func in file_utils.__commands__:
            parser = parser_func(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestTarballMerger(TestCaseWithTmp):
    def setUp(self):
        super(TestTarballMerger, self).setUp()
        self.input_dir = util.file.get_test_input_path(self)

    def test_simple_merge(self):
        """
            Simple round-trip test
        """
        temp_dir = tempfile.gettempdir()
        raw_files = ["file{}".format(x) for x in range(1,5)]

        input_files = [os.path.join(self.input_dir, "compressed-input", x+".tar.gz") for x in raw_files]
        

        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        file_utils.merge_tarballs( out_tarball_file,
                                    input_files                                    
                                    )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",raw_files[i])
            outf = os.path.join(temp_dir,raw_files[i])

            assert_equal_contents(self, inf, outf)
