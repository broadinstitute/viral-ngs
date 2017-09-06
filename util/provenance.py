'''Provenance tracking: for each file, automatically keeping track of how it was created -- from which other files and by which 
code version.'''

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "add_provenance_tracking"]

import hashlib
import inspect
import json
import logging
import os
import os.path
import shutil
import subprocess
import sys
import time
import uuid

import util.file

class FileArg(object):

    '''Argparse parameter type for input and output files, which provides error-checking and optional provenance tracking.'''
    
    def __init__(self, fname, mode):
        self.fname = fname
        self.mode = mode
        assert mode in ('r', 'w')

class InFile(FileArg):

    def __init__(self, fname):
        super(FileArg, self).__init(fname, 'r')

class OutFile(FileArg):

    def __init__(self, fname):
        super(FileArg, self).__init(fname, 'w')

def add_provenance_tracking(cmd_parser, cmd_func):
    '''Add provenance tracking to the given command.
    
    Args:
        cmd_parser - parser for a command defined in a script
        cmd_func - function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.
    '''

    if 'VIRAL_NGS_PROV' not in os.environ: return cmd_func

    prov_dir = os.environ['VIRAL_NGS_PROV']

    join = os.path.join
    dump_file = util.file.dump_file

    time_str = time.strftime('%Y-%m-%d-%H-%M-%S')
    step_id = '-'.join(map(str, (time_str, cmd_func.__module__, cmd_func.__name__, uuid.uuid4())))
    step_dir = join(prov_dir, step_id)

    def _run_cmd_with_tracking(args):

        os.mkdir(step_dir)
        dump_file(join(step_dir, 'env.json'), json.dumps(dict(args=args,
                                                              env=os.environ, argv=sys.argv, start=time_str,
                                                              func_name=cmd_func.__name__,
                                                              func_module=cmd_func.__module__,
                                                              func_file=inspect.getfile(cmd_func))))



        args2 = dict()
        hash_at_end = dict()
        args = vars(args)
        for arg, val in args.items():
            if not isinstance(val, FileArg):
                args2[k] = val
            else:
                args2[k] = val.fname
                if 'r' in val.mode:
                    dump_file(join(step_dir, arg+'.hash'), util.file.hash_file(val.fname, hash_algorithm='sha1'))
                else:
                    hash_at_end[arg] = val.fname

        exc = None
        try: 
            cmd_impl(args2)
        except Exception as e:
            exc = e
        finally:
            dump_file(join(step_dir, 'done.txt'), time.strftime('%Y-%m-%d-%H-%M-%S'))
            if exc is not None:
                dump_file(join(step_dir, 'exception.txt'), str(exc))
                raise exc
            else:
                for arg, fname in hash_at_end.items():
                    dump_file(join(step_dir, arg+'.hash'), util.file.hash_file(fname, hash_algorithm='sha1'))

    return _run_cmd_with_tracking
