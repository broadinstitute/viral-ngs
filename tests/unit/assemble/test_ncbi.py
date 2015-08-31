# Unit tests for ncbi.py

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import unittest, argparse

#module-specific
import ncbi

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in ncbi.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


