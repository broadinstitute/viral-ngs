# Unit tests for interhost.py

__author__ = "PLACEHOLDER"

import interhost
import test
import unittest, shutil, argparse

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in interhost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


@unittest.skip('not implemented')
class TestCoordMapper(test.TestCaseWithTmp):
    def setUp(self):
        super(TestCoordMapper, self).setUp()
        pass
    def test_one_chr_no_gaps(self):
        pass
    def test_one_chr_gaps(self):
        pass
    def test_two_chr_no_gaps(self):
        pass
    
    
    
