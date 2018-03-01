
import functools
import os

import util.file
import util.misc
import util.cmd_plugins

#from util._metadata.file_arg import FileArg
from util._metadata.md_utils import _make_list

class OptionalFile(str):
    """Marker for optional files returned by compute_fnames methods"""

    def __new__(cls, val):
        return str.__new__(cls, val)

#OptionalFile = FileArg.OptionalFile

# ** InFile, OutFile etc

def make_arg_handler(mode, compute_fnames):
    return 

def InFile(val, compute_fnames=_make_list):
    """Argparse argument type for arguments that denote input files."""
    #file_arg = FileArg(val, mode='r', compute_fnames=compute_fnames)
    
    #util.file.check_paths(read=file_arg.fnames)
    #return file_arg
    return util.cmd_plugins.cmd_plugin_mgr.hook.cmd_handle_file_arg(val=val, mode='r', compute_fnames=compute_fnames)

def OutFile(val, compute_fnames=_make_list):
    """Argparse argument type for arguments that denote output files."""
    #file_arg = FileArg(val, mode='w', compute_fnames=compute_fnames)
    #util.file.check_paths(write=file_arg.fnames)
    #return file_arg
    return util.cmd_plugins.cmd_plugin_mgr.hook.cmd_handle_file_arg(val=val, mode='w', compute_fnames=compute_fnames)
#    return make_arg_handler(mode='w', compute_fnames=compute_fnames)

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

def _InFile_OneOf(val, opts):
    for opt in opts:
        try:
            return opt(val)
        except EnvironmentError:
            pass

    raise EnvironmentError("Could not open: {}".format(val))

def InFile_OneOf(*opts):
    """Argparse argument type for an input arg that can be one of several types.  The first type for which all the input
    files exist, will be taken as the right one.
    """
    return functools.partial(_InFile_OneOf, opts=opts)
