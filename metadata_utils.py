#!/usr/bin/env python
''' Utilities for managing metadata.  See util.metadata .
'''

__author__ = "ilya@broadinstitute.org"
__commands__ = []

# built-ins
import argparse
import logging

# intra-module
import util.cmd
from util.metadata import InFile, OutFile, ProvenanceGraph

# third-party

log = logging.getLogger(__name__)

def print_provenance(fname):
    """Print provenance of a given file"""

    G = ProvenanceGraph()
    G.load()

    G.print_provenance(fname)
            
def parser_print_provenance(parser=argparse.ArgumentParser()):
    parser.add_argument('fname', type=InFile, help='File for which to print provenance')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, print_provenance, split_args=True)
    return parser

__commands__.append(('print_provenance', parser_print_provenance))

###############################################

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
