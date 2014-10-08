# Unit tests for ncbi.py

__author__ = "PLACEHOLDER"

import ncbi
import unittest


class TestCommandHelp(unittest.TestCase):
	def test_help_parser_for_each_command(self):
		for cmd_name, main_fun, parser_fun in ncbi.__commands__:
			parser = parser_fun()
			helpstring = parser.format_help()

class TestNCBI(unittest.TestCase):
	def testNothingAtAll(self):
		'''here we test nothing at all and this should pass'''
		pass

