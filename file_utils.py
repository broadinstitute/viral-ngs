#!/usr/bin/env python
"""
Utilities for dealing with files.
"""

__author__ = "tomkinsc@broadinstitute.org"
__commands__ = []

import argparse
import logging

import util.cmd
import util.file

log = logging.getLogger(__name__)


# ==============================
# ***  merge_tarballs   ***
# ==============================

def merge_tarballs(out_tarball, in_tarballs, threads=None, extract_to_disk_path=None, pipe_hint=None):
    ''' Merges separate tarballs into one tarball
        data can be piped in and/or out
    '''
    util.file.repack_tarballs(out_tarball, in_tarballs, threads=threads, extract_to_disk_path=extract_to_disk_path, pipe_hint=pipe_hint)
    return 0
def parser_merge_tarballs(parser=argparse.ArgumentParser()):
    parser.add_argument(
        'out_tarball', 
        help='''output tarball (*.tar.gz|*.tar.lz4|*.tar.bz2|-) 
        Note: if "-" is used, a gzip-compressed tarball 
        will be written to stdout''')
    parser.add_argument(
        'in_tarballs', nargs='+',
        help=('input tarballs (*.tar.gz|*.tar.lz4|*.tar.bz2)')
    )
    parser.add_argument('--extractToDiskPath',
                        dest="extract_to_disk_path",
                        help='If specified, the tar contents will also be extracted to a local directory.')
    parser.add_argument('--pipeHint',
                        dest="pipe_hint",
                        default=".gz",
                        help='If specified, the compression type specifed (as a file extension) is used for input piped in.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, merge_tarballs, split_args=True)
    return parser
__commands__.append(('merge_tarballs', parser_merge_tarballs))


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)