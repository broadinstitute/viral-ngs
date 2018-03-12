#!/usr/bin/env python
"""Dummy CLI defining simple artificial commands for testing util/cmd.py and util/metadata.py"""

__commands__ = []

import os
import os.path
import shutil
import argparse
import time

import util.file
from util.argparse_arg_types import InFile, OutFile, InFiles, OutFiles, InFilesPrefix, OutFilesPrefix, InFile_OneOf

#############################################################################################

def get_file_size(in_fname, size_fname):
    """Writes size of input file to output file"""

    util.file.dump_file(size_fname, os.path.getsize(in_fname))

def parser_get_file_size(parser=argparse.ArgumentParser()):
    parser.add_argument('in_fname', type=InFile, help='File the size of which we want')
    parser.add_argument('size_fname', type=OutFile, help='File to which to write size of file `in_fname`')
    util.cmd.attach_main(parser, get_file_size, split_args=True)
    return parser

__commands__.append(('get_file_size', parser_get_file_size))

#############################################################################################

def get_file_info(args):
    """Gets info about various files.  Weird command with lots of options so we can exercise various aspects of metadata recording."""

    if args.fail: raise RuntimeError('Planned error')

    beg_time = time.time()
    
    all_fnames = args.in_fnames
    if args.in_fnames_pfx:
        all_fnames.extend([args.in_fnames_pfx+'.ex1', args.in_fnames_pfx+'.ex2'])
    if args.in_fnames_dir:
        all_fnames.extend(util.file.files_or_pipes_in_dir(args.in_fnames_dir))

    util.file.dump_file(args.info_fname, sum(map(os.path.getsize, all_fnames)) * args.factor)

    if args.copy_info_to:
        for ext in ('.out1', '.out2'):
            shutil.copyfile(args.info_fname, args.copy_info_to+ext)

    for f in (args.make_empty or ()):
        util.file.make_empty(f)

    end_time = time.time()

    return dict(__metadata__=dict(runtime=end_time-beg_time, nfiles=len(all_fnames)))

def parser_get_file_info(parser=argparse.ArgumentParser()):
    parser.add_argument('in_fnames', type=InFile, nargs='+', help='File(s) the size of which we want')
    parser.add_argument('info_fname', type=OutFile, help='File to which to write info about the files')
    parser.add_argument('--in-fnames-pfx', type=InFilesPrefix(suffixes=('.ex1', '.ex2')), help='Also get info for these files')
    parser.add_argument('--in-fnames-dir', type=InFiles(compute_fnames=util.file.files_or_pipes_in_dir),
                        help='Also get info for files in this dir')
    parser.add_argument('--copy-info-to', type=OutFilesPrefix(suffixes=('.out1', '.out2')), help='Copy results also to these files')
    parser.add_argument('--factor', type=int, default=1, help='Multiply result by this factor')
    parser.add_argument('--make-empty', type=OutFile, action='append', help='Create empty file here')
    parser.add_argument('--fail', action='store_true', help='Raise an exception')
    
    util.cmd.attach_main(parser, get_file_info, split_args=False)
    return parser

__commands__.append(('get_file_info', parser_get_file_info))

#############################################################################################

def concat_input_files(args):
    """Concat the input files"""
    in_fnames = None
    for fns in ([args.inp], [args.inp+'.ex1', args.inp+'.ex2']):
        if all(map(util.file.is_file_or_pipe, fns)):
            in_fnames = fns
    if os.path.isdir(args.inp):
        in_fnames = util.file.files_or_pipes_in_dir(args.inp)

    util.file.concat(in_fnames, args.out)

def parser_concat_input_files(parser=argparse.ArgumentParser()):
    parser.add_argument('inp', type=InFile_OneOf(InFile, InFilesPrefix(suffixes=('.ex1', '.ex2')),
                                                 InFiles(compute_fnames=util.file.files_or_pipes_in_dir)),
                        help='Names of input files')
    parser.add_argument('out', type=OutFile, help='write concatenation of input files to this file')

    util.cmd.attach_main(parser, concat_input_files, split_args=False)
    return parser

__commands__.append(('concat_input_files', parser_concat_input_files))

#############################################################################################

def concat_input_files_list(args):
    """Concat the input files"""
    util.file.concat(args.inp, args.out)

def parser_concat_input_files_list(parser=argparse.ArgumentParser()):
    parser.add_argument('inp', type=InFile, nargs='+',
                        help='Names of input files')
    parser.add_argument('out', type=OutFile, help='write concatenation of input files to this file')

    util.cmd.attach_main(parser, concat_input_files_list, split_args=False)
    return parser

__commands__.append(('concat_input_files_list', parser_concat_input_files_list))

#############################################################################################



def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
