#!/usr/bin/env python

"""Commands for working with sets of kmers"""


from __future__ import print_function
__author__ = "ilya@broadinstitute.org"

__commands__ = []

import argparse
import glob
import logging
import os
import tempfile
import shutil
import functools

from Bio import SeqIO
import pysam

import util.cmd
import util.file
import util.misc
import tools
import tools.picard
import tools.samtools
import tools.kmc
import read_utils

log = logging.getLogger(__name__)

# =================================

def build_kmer_db(seq_files, kmer_db, kmer_size=tools.kmc.DEFAULT_KMER_SIZE, min_occs=None, max_occs=None,
                  counter_cap=tools.kmc.DEFAULT_COUNTER_CAP, mem_limit_gb=8, threads=None):
    """Build a database of kmers occurring in given sequences."""
    tools.kmc.KmcTool().build_kmer_db(seq_files=seq_files, kmer_size=kmer_size, min_occs=min_occs, max_occs=max_occs, counter_cap=counter_cap,
                                      kmer_db=kmer_db, mem_limit_gb=mem_limit_gb, threads=threads)

def parser_build_kmer_db(parser=argparse.ArgumentParser()):
    parser.add_argument('seq_files', nargs='+', help='Files from which to extract kmers (fasta/fastq/bam, fasta/fastq may be .gz or .bz2)')
    parser.add_argument('kmer_db', help='kmer database (with or without .kmc_pre/.kmc_suf suffix)')
    parser.add_argument('--kmerSize', '-k', dest='kmer_size', type=int, help='kmer size')
    parser.add_argument('--minOccs', '-ci', dest='min_occs', type=int, help='drop kmers with fewer than this many occurrences')
    parser.add_argument('--maxOccs', '-cx', dest='max_occs', type=int, help='drop kmers with more than this many occurrences')
    parser.add_argument('--counterCap', '-cs', dest='counter_cap', type=int, default=tools.kmc.DEFAULT_COUNTER_CAP, help='cap kmer counts at this value')
    parser.add_argument('--memLimitGb', dest='mem_limit_gb', default=8, type=int, help='Max memory to use, in GB')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, build_kmer_db, split_args=True)
    return parser

__commands__.append(('build_kmer_db', parser_build_kmer_db))

# =========================

def dump_kmers(kmer_db, out_kmers, min_occs=None, max_occs=None, threads=None):
    """Dump kmers from kmer database to a text file"""
    tools.kmc.KmcTool().dump_kmers(kmer_db=kmer_db, out_kmers=out_kmers, min_occs=min_occs, max_occs=max_occs, threads=threads)

def parser_dump_kmers(parser=argparse.ArgumentParser()):
    parser.add_argument('kmer_db', help='kmer database (with or without .kmc_pre/.kmc_suf suffix)')
    parser.add_argument('out_kmers', help='text file to which to write the kmers')
    parser.add_argument('--minOccs', '-ci', dest='min_occs', type=int, help='drop kmers with fewer than this many occurrences')
    parser.add_argument('--maxOccs', '-cx', dest='max_occs', type=int, help='drop kmers with more than this many occurrences')

    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, dump_kmers, split_args=True)
    return parser

__commands__.append(('dump_kmers', parser_dump_kmers))

# =========================

def filter_by_kmers(kmer_db, in_reads, out_reads, db_min_occs=None, db_max_occs=None, 
                    read_min_occs=None, read_max_occs=None, hard_mask=False, threads=None):
    """Filter sequences based on their kmer contents."""
    print('filter_by_kmers: starting run')
    tools.kmc.KmcTool().filter_reads(kmer_db=kmer_db, in_reads=in_reads, out_reads=out_reads, db_min_occs=db_min_occs, db_max_occs=db_max_occs,
                                     read_min_occs=read_min_occs, read_max_occs=read_max_occs, hard_mask=hard_mask,
                                     threads=threads)

def parser_filter_by_kmers(parser=argparse.ArgumentParser()):
    parser.add_argument('kmer_db', help='kmer database (with or without .kmc_pre/.kmc_suf suffix)')
    parser.add_argument('in_reads', help='input reads, as fasta/fastq/bam')
    parser.add_argument('out_reads', help='output reads')
    parser.add_argument('--dbMinOccs', dest='db_min_occs', type=int, help='ignore datatbase kmers with count below this')
    parser.add_argument('--dbMaxOccs', dest='db_max_occs', type=int, help='ignore datatbase kmers with count above this')
    int_or_float = functools.partial(util.misc.as_type, **dict(types=(int, float)))
    parser.add_argument('--readMinOccs', dest='read_min_occs', type=int_or_float,
                        help='filter out reads with fewer than this many db kmers; if a float, interpreted as fraction of read length')
    parser.add_argument('--readMaxOccs', dest='read_max_occs', type=int_or_float,
                        help='filter out reads with more than this many db kmers; if a float, interpreted as fraction of read length')
    parser.add_argument('--hardMask', dest='hard_mask', default=False, action='store_true', help='In the output reads, mask the invalid kmers')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_by_kmers, split_args=True)
    return parser

__commands__.append(('filter_by_kmers', parser_filter_by_kmers))

# =========================

def kmers_binary_op(op, kmer_db1, kmer_db2, kmer_db_out, threads=None):
    """Perform a simple binary operation on kmer sets."""

    tools.kmc.KmcTool().kmers_binary_op(op, kmer_db1, kmer_db2, kmer_db_out, threads=threads)

def parser_kmers_binary_op(parser=argparse.ArgumentParser()):
    parser.add_argument('op', choices=('intersect', 'union', 'kmers_subtract', 'counters_subtract'), help='binary operation to perform')
    parser.add_argument('kmer_db1', help='first kmer set')
    parser.add_argument('kmer_db2', help='second kmer set')
    parser.add_argument('kmer_db_out', help='output kmer db')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kmers_binary_op, split_args=True)
    return parser

__commands__.append(('kmers_binary_op', parser_kmers_binary_op))

# =========================

def kmers_set_counts(kmer_db_in, value, kmer_db_out, threads=None):
    """Copy the kmer database, setting all kmer counts in the output to the given value."""

    tools.kmc.KmcTool().set_kmer_counts(kmer_db_in, value, kmer_db_out, threads=threads)

def parser_kmers_set_counts(parser=argparse.ArgumentParser()):
    parser.add_argument('kmer_db_in', help='input kmer db')
    parser.add_argument('value', type=int, help='all kmer counts in the output will be set to this value')
    parser.add_argument('kmer_db_out', help='output kmer db')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kmers_set_counts, split_args=True)
    return parser

__commands__.append(('kmers_set_counts', parser_kmers_set_counts))

# ========================

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
