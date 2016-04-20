# Unit tests for metagenomics.py

import unittest
import argparse
import argparse
import metagenomics

class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in metagenomics.__commands__:
            print(cmd_name)
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()
