
import functools
import os

import util.file
import util.misc
import util.cmd_plugins

class OptionalFile(str):
    """Marker for optional files returned by compute_fnames methods"""

    def __new__(cls, val):
        return str.__new__(cls, val)

# ** InFile, OutFile etc

def _call_arg_handler(val, mode, compute_fnames=util.misc.make_list):
    """Call any hook implementations defined to handle argparse args that denote input or output file(s).
    But first, check that input files are readable and output files are writable; throw an exception if not.
    Besides detecting problems early, before wasting resources to execute a command destined to fail,
    this check is used by InFile_OneOf() to determine what input files are actually denoted by `val`,
    when this depends on what input files exist.
    """
    required_fnames = [str(f) for f in compute_fnames(val) if not isinstance(f, OptionalFile)]
    rw = {'r':'read', 'w':'write'}[mode]
    util.file.check_paths(**{rw : required_fnames})

    return util.cmd_plugins.cmd_plugin_mgr.hook.cmd_handle_file_arg(val=val, mode=mode, compute_fnames=compute_fnames)

InFile = functools.partial(_call_arg_handler, mode='r')
OutFile = functools.partial(_call_arg_handler, mode='w')

def InFiles(compute_fnames):
    """Argparse argument type for a string from which names of input files can be computed"""
    return functools.partial(InFile, compute_fnames=compute_fnames)

def OutFiles(compute_fnames):
    """Argparse argument type for a string from which names of output files can be computed"""
    return functools.partial(OutFile, compute_fnames=compute_fnames)

def InFilesPrefix(suffixes):
    """Argparse argument type for a string that denotes the common prefix of a group of input files with a known list of suffixes."""
    return InFiles(compute_fnames=functools.partial(util.misc.add_suffixes, suffixes=suffixes))

def OutFilesPrefix(suffixes):
    """Argparse argument type for a string that denotes the common prefix of a group of input files with a known list of suffixes."""
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

def from_cmd_return_value(val, cmd_return_value=None):
    return cmd_return_value or []
