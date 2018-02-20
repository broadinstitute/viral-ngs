#!/usr/bin/env python
"""Dummy CLI defining simple commands for testing util/cmd.py and util/metadata.py"""

__commands__ = []

import os
import os.path
import argparse

import util.file
from util.metadata import InFile, OutFile

def get_file_size(in_fname, size_fname):
    """Writes size of input file to output file"""

    util.file.dump_file(size_fname, os.path.getsize(in_fname))

def parser_get_file_size(parser=argparse.ArgumentParser()):
    parser.add_argument('in_fname', type=InFile, help='File the size of which we want')
    parser.add_argument('size_fname', type=OutFile, help='File to which to write size of file `in_fname`')
    util.cmd.attach_main(parser, get_file_size, split_args=True)
    return parser

__commands__.append(('get_file_size', parser_get_file_size))

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
