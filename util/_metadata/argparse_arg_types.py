
import functools

import util.file
import util.misc

from util._metadata.recording import FileArg
from util._metadata.md_utils import _make_list

# ** InFile, OutFile etc

def InFile(val, compute_fnames=_make_list):
    """Argparse argument type for arguments that denote input files."""
    file_arg = FileArg(val, mode='r', compute_fnames=compute_fnames)
    util.file.check_paths(read=file_arg.fnames)
    return file_arg

def OutFile(val, compute_fnames=_make_list):
    """Argparse argument type for arguments that denote output files."""
    file_arg = FileArg(val, mode='w', compute_fnames=compute_fnames)
    util.file.check_paths(write=file_arg.fnames)
    return file_arg

def InFiles(compute_fnames):
    """Argparse argument type for a string from which names of input files can be computed"""
    return functools.partial(InFile, compute_fnames=compute_fnames)

def OutFiles(compute_fnames):
    """Argparse argument type for a string from which names of output files can be computed"""
    return functools.partial(OutFile, compute_fnames=compute_fnames)

def InFilesPrefix(suffixes):
    """Argparse argument type for a string that denotes the common prefix of a group of input files."""
    return InFiles(compute_fnames=functools.partial(util.misc.add_suffixes, suffixes=suffixes))

def OutFilesPrefix(suffixes):
    """Argparse argument type for a string that denotes the common prefix of a group of input files."""
    return OutFiles(compute_fnames=functools.partial(util.misc.add_suffixes, suffixes=suffixes))

