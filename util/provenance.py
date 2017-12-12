'''Provenance tracking: for each file, automatically keeping track of how it was created -- from which other files and by which 
code version.  Also facilitates uniform representation of operation parameters and result metrics, enabling various queries.

See data_utils.py for utils that create useful reports from provenance data.
'''

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "add_provenance_tracking", "is_provenance_tracking_enabled", "VIRAL_NGS_PROVENANCE_LOCATION"]

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
import json
import traceback
import inspect
import collections

# intra-module
import util.file
import util.misc
import util.version

# third-party
import concurrent.futures

_log = logging.getLogger(__name__)


VIRAL_NGS_PROVENANCE_DATA_FORMAT='1.0.0'
VIRAL_NGS_PROVENANCE_LOCATION='VIRAL_NGS_PROVENANCE_LOCATION'

def is_provenance_tracking_enabled():
    return VIRAL_NGS_PROVENANCE_LOCATION in os.environ

class FileArg(str):

    '''Argparse parameter type for input and output files, which provides error-checking and optional provenance tracking.'''
    
    def __new__(cls, *args, **kw):
        return str.__new__(cls, *args, **kw)

class InFile(FileArg):

    def __new__(cls, *args, **kw):
        return FileArg.__new__(cls, *args, **kw)

class OutFile(FileArg):

    def __new__(cls, *args, **kw):
        return FileArg.__new__(cls, *args, **kw)

if not is_provenance_tracking_enabled():
    InFile = str
    OutFile = str

class Hasher(object):
    """Manages computation of file hashes.

    It might also cache the actual contents of files, depending on size.
    """

    def __init__(self):
        pass

    def __call__(self, file):
        return util.file.hash_file(file, hash_algorithm='sha1')

def add_provenance_tracking(cmd_parser, cmd_func, cmd_module, cmd_name):
    """Add provenance tracking to the given command.
    
    Args:
        cmd_parser: parser for a command defined in a script
        cmd_func: function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.

    Returns:
        a wrapper for cmd_func, which has the same signature but adds provenance tracking (if provenance tracking is configured)
    """

    def _run_cmd_with_tracking(args):

        cmd_exception = None
        cmd_result = None

        try:
            beg_time = time.time()
            cmd_result = cmd_func(args)
        except Exception as e:
            cmd_exception = e
        finally:
            try:
                if is_provenance_tracking_enabled():
                    end_time = time.time()

                    step_id = '__'.join(map(str, (time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime(beg_time)), 
                                                  cmd_module, cmd_name, uuid.uuid4()))).replace('-','_')

                    pgraph = dict(format=VIRAL_NGS_PROVENANCE_DATA_FORMAT, in_files={}, out_files={})
                    hasher = Hasher()

                    args_dict = vars(args)
                    args_dict.pop('func_main', None)
                    pgraph['step']=dict(step_id=step_id, run_id=os.get('VIRAL_NGS_RUN_ID', ''), 
                                        cmd_module=cmd_module, cmd_name=cmd_name,
                                        beg_time=beg_time, end_time=end_time, duration=end_time-beg_time,
                                        exception=str(cmd_exception),
                                        viral_ngs_version=util.version.get_version(),
                                        platform=platform.platform(), cpus=util.misc.available_cpu_count(), host=socket.getfqdn(), 
                                        user=getpass.getuser(),
                                        cwd=os.getcwd(),
                                        argv=tuple(sys.argv),
                                        args=args_dict, **(cmd_result if isinstance(cmd_result, collections.Mapping) else {}))
                    
                    for arg, val in args_dict.items():
                        for i, v in enumerate(util.misc.make_seq(val)):
                            if v and isinstance(v, FileArg) and (isinstance(v, InFile) or cmd_exception is None) and \
                               os.path.isfile(v) and not stat.S_ISFIFO(os.stat(v).st_mode):
                                file_hash = hasher(v)
                                pgraph['in_files' if isinstance(v, InFile) else 'out_files'][file_hash] = \
                                    dict(fname=v, realpath=os.path.realpath(v), arg=arg, arg_order=i)

                    util.file.dump_file(os.path.join(os.environ['VIRAL_NGS_PROV'], step_id+'.json'),
                                        json.dumps(pgraph, sort_keys=True, indent=4))

            except Exception:
                _log.warning('Error recording provenance ({})'.format(traceback.format_exc()))
                
        if cmd_exception:
            raise cmd_exception

    return _run_cmd_with_tracking
