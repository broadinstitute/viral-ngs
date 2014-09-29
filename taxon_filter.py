#!/usr/bin/env python
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.

Requires python >= 2.7 and BioPython.  On the Broad cluster, it is known
to work with the Python-2.7 and Python-3.4 dotkits.'''

__author__ = "PLACEHOLDER"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging
import util.cmd, util.file, util.vcf, util.misc

log = logging.getLogger(__name__)
global_tool_paths = {}


def parser_trim_trimmomatic():
	parser = argparse.ArgumentParser(
		description='''Trim read sequences with Trimmomatic. Perhaps move this to
		a separate script of general bam/alignment utility functions?''')
	parser.add_argument("inBam", help="Input BAM file")
	parser.add_argument("outBam", help="Output BAM file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_trim_trimmomatic(args):
	raise ("not yet implemented")
	return 0
__commands__.append(('trim_trimmomatic', main_trim_trimmomatic, parser_trim_trimmomatic))


def parser_filter_lastal():
	parser = argparse.ArgumentParser(
		description='''Restrict input reads to those that align to the given
		reference databases using LASTAL.''')
	parser.add_argument("inBam", help="Input BAM file")
	parser.add_argument("refDbs", nargs='+',
		help="""Reference databases (one or more) to retain from input""")
	parser.add_argument("outBam", help="Output BAM file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_filter_lastal(args):
	raise ("not yet implemented")
	return 0
__commands__.append(('filter_lastal', main_filter_lastal, parser_filter_lastal))


def parser_deplete_bmtagger():
	parser = argparse.ArgumentParser(
		description='''Deplete human reads and other contaminants using bmtagger''')
	parser.add_argument("inBam", help="Input BAM file")
	parser.add_argument("refDbs", nargs='+',
		help="""Reference databases (one or more) to deplete from input""")
	parser.add_argument("outBam", help="Output BAM file")
	util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
	return parser
def main_deplete_bmtagger(args):
	raise ("not yet implemented")
	return 0
__commands__.append(('deplete_bmtagger', main_deplete_bmtagger, parser_deplete_bmtagger))


if __name__ == '__main__':
	util.cmd.main_argparse(__commands__, __doc__)
