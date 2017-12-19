'''Automated recording of data provenance and metrics.  This lets you answer questions like: how was a given data file made?  Which workflow
versions tend to produce the best results, according to given metrics?

This module deals with data recording; ../data_utils.py deals with data querying.

The unit of recording in this module is one command (see cmd.py), as opposed to one Python function or one whole workflow.
To use this module, import InFile and OutFile from it; then, when defining argparse arguments for a command, in add_argument()
use "type=InFile" for input files or "type=OutFile" for output files. Until provenance tracking is enabled (see below how),
InFile and OutFile are defined to be simply str.  To enable provenance tracking, set the environment variable
VIRAL_NGS_METADATA to a writable directory.  Then, whenever a command is run, a file is written to that directory
recording the command, all its parameters, the input and output files, and details of the run environment.
Input and output files are identified by a hash of their contents, so that copied/moved/renamed files can be properly identified.
For each file, we also record metadata such as name, size, and modification date.


: for each data file, automatically keeping track of how it was created -- from which inputs, by which 
code version and with what parameters.


Note that only files actually used are listed in the metadata; optional files not specified in a given command invocation are not.


See data_utils.py for utils that create useful reports from provenance data.
'''

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "add_metadata_tracking", "is_metadata_tracking_enabled", "metadata_dir"]

# built-ins
import argparse
import logging
import os
import os.path
import sys
import stat
import platform
import shutil
import sys
import time
import uuid
import socket
import getpass
import pwd
import json
import traceback
import inspect
import collections
import re

# intra-module
import util.file
import util.misc
import util.version

# third-party
import concurrent.futures

_log = logging.getLogger(__name__)


VIRAL_NGS_METADATA_FORMAT='1.0.0'

def is_metadata_tracking_enabled():
    return 'VIRAL_NGS_METADATA' in os.environ

def metadata_dir():
    return os.environ['VIRAL_NGS_METADATA']

class FileArg(str):

    '''Argparse parameter type for input and output files, which provides error-checking and optional provenance tracking.'''
    
    def __new__(cls, *args, **kw):
        return str.__new__(cls, *args, **kw)

class InFile(FileArg):

    def __new__(cls, *args, **kw):

        # check that the file is readable, force an immediate error if it is not
        with open(args[0]):
            pass

        return FileArg.__new__(cls, *args, **kw)

class OutFile(FileArg):

    def __new__(cls, *args, **kw):
        fname = args[0]
        if not os.path.exists(fname):
            with open(fname, 'w'):
                pass
            os.unlink(fname)
        else:
            if not (os.path.isfile(fname) and os.access(fname, os.W_OK)):
                raise IOError('Output filel not writable: ' + fname)

        return FileArg.__new__(cls, *args, **kw)

if not is_metadata_tracking_enabled():
    InFile = str
    OutFile = str

class Hasher(object):
    """Manages computation of file hashes.

    It might also cache the actual contents of files, depending on size.
    """

    def __init__(self, hash_algorithm='sha1'):
        self.hash_algorithm = hash_algorithm

    def __call__(self, file):
        return self.hash_algorithm + '_' + util.file.hash_file(file, hash_algorithm=self.hash_algorithm)


def create_run_id(t=None):
    """Generate a unique ID for a run (set of steps run as part of one workflow)"""
    return util.file.string_to_file_name('__'.join(map(str, (time.strftime('%Y%m%d%H%M%S', time.localtime(t))[2:], getpass.getuser(),
                                                             os.path.basename(os.getcwd()), uuid.uuid4()))))[:210]

def set_run_id():
    """Generate and record in the environment a unique ID for a run (set of steps run as part of one workflow)."""
    os.environ['VIRAL_NGS_METADATA_RUN_ID'] = create_run_id()

def tag_code_version(tag, push_to=None):
    """Create a lightweight git tag for the current state of the project repository, even if the state is dirty.
    If the repository is dirty, use the 'git stash create' command to create a commit representing the current state,
    and tag that; else, tag the existing clean state.  If `push_to` is not None, push the tag to the specified remote.
    """
    def run_cmd(cmd):
        out = ''
        result = util.misc.run_and_print(cmd.strip().split(), silent=True)
        if result.returncode == 0:
            out = result.stdout
            if not isinstance(out, str):
                out = out.decode('utf-8')
            out = out.strip()
        return out

    with util.file.pushd_popd(util.version.get_project_path()):
        run_cmd('git tag ' + tag + ' ' + run_cmd('git stash create'))
        if push_to:
            run_cmd('git push ' + push_to + ' ' + tag)
    
def add_metadata_tracking(cmd_parser, cmd_main, cmd_main_orig):
    """Add provenance tracking to the given command.
    
    Args:
        cmd_parser: parser for a command defined in a script
        cmd_main: function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.
        cmd_main_orig: the original cmd_main function defined in the command script, if cmd_main is a wrapper around the original 
             function, else same as cmd_main

    Returns:
        a wrapper for cmd_main, which has the same signature but adds provenance tracking (if provenance tracking is configured)
    """
    if not is_metadata_tracking_enabled():
        return cmd_main

    cmd_module=os.path.splitext(os.path.basename(inspect.getfile(cmd_main_orig)))[0]
    cmd_name=cmd_main_orig.__name__

    cmd_parser.add_argument('--metadata', action='append', nargs=2, metavar=('ATTRIBUTE', 'VALUE'), help='attach metadata to step')
    cmd_parser.add_argument('--file-metadata', action='append', nargs=3, metavar=('FILE', 'ATTRIBUTE', 'VALUE'), 
                            help='attach metadata to input or output file')

    def _run_cmd_with_tracking(args):

        cmd_exception, cmd_exception_str, cmd_result = (None,)*3

        args_dict = vars(args)
        args_dict.pop('func_main', None)

        cmd_metadata = dict(args_dict.get('metadata', []))
        if hasattr(args, 'metadata'):
            delattr(args, 'metadata')

        file_metadata = args_dict.get('file_metadata', [])
        if hasattr(args, 'file_metadata'):
            delattr(args, 'file_metadata')

        try:
            beg_time = time.time()

            # *** Run the actual command ***
            cmd_result = cmd_main(args)
        except Exception as e:
            cmd_exception = e
            cmd_exception_str = traceback.format_exc()
        finally:
            try:
                if is_metadata_tracking_enabled():
                    end_time = time.time()

                    run_id = os.environ.get('VIRAL_NGS_METADATA_RUN_ID', create_run_id(beg_time))
                    step_id = '__'.join(map(str, (create_run_id(beg_time), cmd_module, cmd_name)))

                    code_repo = os.path.join(metadata_dir(), 'code_repo')
                    tag_code_version('cmd_' + step_id, push_to=code_repo if os.path.isdir(code_repo) else None)

                    pgraph = dict(format=VIRAL_NGS_METADATA_FORMAT, in_files={}, out_files={})
                    hasher = Hasher()


                    # Sometimes file names accessed by a command are not passed as command args, but computed within the command.
                    # The command can tell us about such files by returning a dict with the key 'files' mapped to a dict
                    # containing additional (arg name, file name(s)) mappings.
                    info_from_cmd = cmd_result if isinstance(cmd_result, collections.Mapping) else {}
                    args_dict.update(info_from_cmd.get('files', {}))

                    pgraph['step']=dict(step_id=step_id, run_id=run_id,
                                        metadata_dir=metadata_dir(),
                                        cmd_module=cmd_module, cmd_name=cmd_name,
                                        beg_time=beg_time, end_time=end_time, duration=end_time-beg_time,
                                        exception=cmd_exception_str,
                                        viral_ngs_version=util.version.get_version(),
                                        viral_ngs_path=util.version.get_project_path(),
                                        viral_ngs_path_real=os.path.realpath(util.version.get_project_path()),
                                        platform=platform.platform(), cpus=util.misc.available_cpu_count(), host=socket.getfqdn(),
                                        user=getpass.getuser(),
                                        cwd=os.getcwd(),
                                        argv=tuple(sys.argv),
                                        args=args_dict,
                                        **cmd_metadata)

                    # The command can, through its return value, pass us metadata to attach either to input/output files or to the
                    # step itself (the latter represented by the key of None)
                    file2metadata = info_from_cmd.get('metadata', {})
                    pgraph['step'].update(file2metadata.get(None, {}))
                    for file, attr, val in file_metadata:
                        file_metadata.setdefault(file, {})[attr] = val

                    for arg, val in args_dict.items():
                        for i, v in enumerate(util.misc.make_seq(val)):
                            if v and isinstance(v, FileArg) and (isinstance(v, InFile) or cmd_exception is None) and \
                               os.path.isfile(v) and not stat.S_ISFIFO(os.stat(v).st_mode):

                                # Metadata for a file is considered to reside not on the file node, but 
                                # on the edge between the file node and the step node.
                                # That way, if multiple steps read the same input or produce the same output, metadata from the
                                # different steps will remain separate ("same" here refers to file contents rather than identity).

                                file_info = dict(fname=v, realpath=os.path.realpath(v), arg=arg, arg_order=i, **file2metadata.get(v, {}))
                                try:
                                    file_stat = os.stat(v)
                                    file_info.update(size=file_stat[stat.ST_SIZE])
                                    file_info.update(owner=pwd.getpwuid(file_stat[stat.ST_UID]).pw_name)
                                except Exception:
                                    _log.warning('Error getting file info ({})'.format(traceback.format_exc()))
                                                     
                                pgraph['in_files' if isinstance(v, InFile) else 'out_files'][hasher(v)] = file_info

                    util.file.dump_file(os.path.join(metadata_dir(), step_id+'.json'),
                                        json.dumps(pgraph, sort_keys=True, indent=4))

            except Exception:
                _log.warning('Error recording metadata ({})'.format(traceback.format_exc()))
                
        if cmd_exception:
            raise cmd_exception

    return _run_cmd_with_tracking

if __name__ == '__main__':
    print(os.path.realpath(util.version.get_project_path()))
