# Unit tests for ncbi.py

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import os
import sys
import tempfile
import unittest
import argparse
import tarfile
import subprocess

# third-party
import pytest
from mock import patch

# module-specific
import file_utils
import util.file
from test import TestCaseWithTmp, assert_equal_contents


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
        self.input_tgz_files = [os.path.join(self.input_dir, "compressed-input", x+".tar.gz") for x in self.raw_files]
    
    def test_simple_merge(self):
        """
            Simple repack test
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        file_utils.merge_tarballs( out_tarball_file,
                                    self.input_tgz_files
                                    )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            assert_equal_contents(self, inf, outf)


    def test_merge_with_extract(self):
        """
            Test streaming repack with intermediate extraction to disk
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")
        out_extracted_path = os.path.join(temp_dir,"extracted")

        file_utils.merge_tarballs( out_tarball_file,
                                    self.input_tgz_files,
                                    extract_to_disk_path=out_extracted_path
                                    )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        # inspect merged
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            assert_equal_contents(self, inf, outf)

        # inspect extracted
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(out_extracted_path,self.raw_files[i])

            assert_equal_contents(self, inf, outf)

    def test_merge_with_extract_repack_from_disk(self):
        """
            Test with repack from disk source after extraction
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")
        out_extracted_path = os.path.join(temp_dir,"extracted")

        util.file.repack_tarballs( out_tarball_file,
                                    self.input_tgz_files,
                                    extract_to_disk_path=out_extracted_path,
                                    avoid_disk_roundtrip=False
                                    )

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        # inspect merged
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            assert_equal_contents(self, inf, outf)

        # inspect extracted
        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(out_extracted_path,self.raw_files[i])

            assert_equal_contents(self, inf, outf)


    def test_piped_in_merge(self):
        """
            Test with streamed input
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        ps = subprocess.Popen("cat {files}".format(files=' '.join(self.input_tgz_files)).split(), stdout=subprocess.PIPE)
        with patch('sys.stdin', ps.stdout):
            file_utils.merge_tarballs( out_tarball_file,
                                        ["-"],
                                        pipe_hint_in="gz" )
        ps.wait()

        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            assert_equal_contents(self, inf, outf)

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_piped_out_merge(self):
        """
            Test with streamed output
        """
        temp_dir = tempfile.gettempdir()
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        with open(out_tarball_file, "wb", 0) as outf:
            # temporarily disable pytest's capture of sys.stdout
            with self.capsys.disabled():
                with patch('sys.stdout', outf):
                    file_utils.merge_tarballs( "-",
                                                self.input_tgz_files,
                                                pipe_hint_out="gz" )
        
        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            assert_equal_contents(self, inf, outf)

    def test_merge_piped_in_and_out(self):
        """
            Test with streamed input and output
        """
        temp_dir = tempfile.gettempdir()
        
        out_tarball_file = os.path.join(temp_dir,"out.tar.gz")

        ps = subprocess.Popen("cat {files}".format(files=' '.join(self.input_tgz_files)).split(), stdout=subprocess.PIPE)
        with patch('sys.stdin', ps.stdout):
            with open(out_tarball_file, "wb", 0) as outf:
                # temporarily disable pytest's capture of sys.stdout
                with self.capsys.disabled():
                    with patch('sys.stdout', outf):
                        file_utils.merge_tarballs( "-",
                                                  ["-"],
                                                  pipe_hint_out="gz",
                                                  pipe_hint_in="gz" )
        ps.wait()
        
        tb = tarfile.open(out_tarball_file)
        tb.extractall(path=temp_dir)

        for i in range(len(self.raw_files)):
            inf = os.path.join(self.input_dir,"raw-input",self.raw_files[i])
            outf = os.path.join(temp_dir,self.raw_files[i])

            assert_equal_contents(self, inf, outf)

