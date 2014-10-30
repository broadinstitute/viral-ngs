# Unit tests for intrahost.py

__author__ = "PLACEHOLDER"

import intrahost
import unittest

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, main_fun, parser_fun in intrahost.__commands__:
            parser = parser_fun()
            helpstring = parser.format_help()
