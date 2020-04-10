#!/usr/bin/env python3
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

def merge_tarballs(out_tarball, in_tarballs, threads=None, extract_to_disk_path=None, pipe_hint_in=None, pipe_hint_out=None):
    ''' Merges separate tarballs into one tarball
        data can be piped in and/or out
    '''
    util.file.repack_tarballs(out_tarball, in_tarballs, threads=threads, extract_to_disk_path=extract_to_disk_path, pipe_hint_in=pipe_hint_in, pipe_hint_out=pipe_hint_out)
    return 0
def parser_merge_tarballs(parser=argparse.ArgumentParser()):
    parser.add_argument(
        'out_tarball', 
        help='''output tarball (*.tar.gz|*.tar.lz4|*.tar.bz2|*.tar.zst|-);
                compression is inferred by the file extension.
        Note: if "-" is used, output will be written to stdout and
         --pipeOutHint must be provided to indicate compression type
         when compression type is not gzip (gzip is used by default).
        ''')
    parser.add_argument(
        'in_tarballs', nargs='+',
        help=('input tarballs (*.tar.gz|*.tar.lz4|*.tar.bz2|*.tar.zst)')
    )
    parser.add_argument('--extractToDiskPath',
                        dest="extract_to_disk_path",
                        help='If specified, the tar contents will also be extracted to a local directory.')
    parser.add_argument('--pipeInHint',
                        dest="pipe_hint_in",
                        default="gz",
                        help='If specified, the compression type used is used for piped input.')
    parser.add_argument('--pipeOutHint',
                        dest="pipe_hint_out",
                        default="gz",
                        help='If specified, the compression type used is used for piped output.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, merge_tarballs, split_args=True)
    return parser
__commands__.append(('merge_tarballs', parser_merge_tarballs))


def parser_rename_fasta_sequences(parser=argparse.ArgumentParser()):
    parser.add_argument("in_fasta", help="input fasta sequences")
    parser.add_argument("out_fasta", help="output (renamed) fasta sequences")
    parser.add_argument("new_name", help="new sequence base name")
    parser.add_argument(
        "--suffix_always",
        help="append numeric index '-1' to <new_name> if only one sequence exists in <input> (default: %(default)s)",
        default=False,
        action="store_true",
        dest="suffix_always"
    )

    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_rename_fasta_sequences)
    return parser
def main_rename_fasta_sequences(args):
    ''' Renames the sequences in a fasta file. Behavior modes:
        1. If input file has exactly one sequence and suffix_always is False,
            then the output file's sequence is named new_name.
        2. In all other cases,
            the output file's sequences are named <new_name>-<i> where <i> is an increasing number from 1..<# of sequences>
    '''
    n_seqs = util.file.fasta_length(args.in_fasta)
    with open(args.in_fasta, 'rt') as inf:
      with open(args.out_fasta, 'wt') as outf:
        if (n_seqs == 1) and not args.suffix_always:
          inf.readline()
          outf.write('>' + args.new_name + '\n')
          for line in inf:
            outf.write(line)
        else:
          i = 1
          for line in inf:
            if line.startswith('>'):
              line = args.new_name + '-' + str(i) + '\n'
            outf.write(line)

    return 0
__commands__.append(('rename_fasta_sequences', parser_rename_fasta_sequences))


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
