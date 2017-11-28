'''Provenance tracking: for each file, automatically keeping track of how it was created -- from which other files and by which 
code version.  Also facilitates uniform tracking of '''

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "add_provenance_tracking"]

# built-ins
import argparse
import collections
import hashlib
import functools
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

# intra-module
import util.file
import util.misc
import util.version

# third-party
import networkx

class FileArg(object):

    '''Argparse parameter type for input and output files, which provides error-checking and optional provenance tracking.'''
    
    def __init__(self, fname, mode):
        self.fname = fname
        self.mode = mode
        assert mode in ('r', 'w')

class InFile(FileArg):

    def __init__(self, fname):
        super(InFile, self).__init__(fname, 'r')

class OutFile(FileArg):

    def __init__(self, fname):
        super(OutFile, self).__init__(fname, 'w')



def add_provenance_tracking(cmd_parser, cmd_func):
    '''Add provenance tracking to the given command.
    
    Args:
        cmd_parser - parser for a command defined in a script
        cmd_func - function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.
    '''

    def dict2namespace(d):
        ns = argparse.Namespace()
        for k, v in d.items():
            setattr(ns, k, v)
        return ns

    @functools.wraps(cmd_func)
    def _unwrap_file_args(args):
        cmd_func(dict2namespace(dict((arg, val.fname if isinstance(val, FileArg) else val) for arg, val in vars(args).items())))

    if 'VIRAL_NGS_PROV' not in os.environ: return _unwrap_file_args

    prov_dir = os.environ['VIRAL_NGS_PROV']

    join = os.path.join
    dump_file = util.file.dump_file

    def _run_cmd_with_tracking(args):

        time_str = time.strftime('%Y-%m-%d-%H-%M-%S')
        step_id = '-'.join(map(str, (time_str, cmd_func.__module__, cmd_func.__name__, uuid.uuid4())))
        step_dir = join(prov_dir, step_id)

        os.mkdir(step_dir)
        os.mkdir(join(step_dir, 'inputs'))
        os.mkdir(join(step_dir, 'outputs'))
        os.mkdir(join(step_dir, 'metrics'))

        args2 = dict()
        hash_at_end = dict()
        args = vars(args)
        for arg, val in args.items():
            if not isinstance(val, FileArg):
                args2[arg] = val
            else:
                args2[arg] = val.fname
                if 'r' in val.mode:
                    dump_file(join(step_dir, 'inputs', arg+'.hash'), util.file.hash_file(val.fname, hash_algorithm='sha1'))
                else:
                    hash_at_end[arg] = val.fname

        dump_file(join(step_dir, 'args.json'), json.dumps(dict(args=str(args2))))
        dump_file(join(step_dir, 'version.json'), json.dumps(dict(viral_ngs_version=util.version.get_version())))
        dump_file(join(step_dir, 'env.json'), json.dumps(dict(env=dict(os.environ))))
        dump_file(join(step_dir, 'cpus.json'), json.dumps(dict(cpus=util.misc.available_cpu_count())))
        dump_file(join(step_dir, 'argv.json'), json.dumps(dict(argv=tuple(sys.argv))))
        dump_file(join(step_dir, 'start.json'), json.dumps(dict(start=time_str)))

        exc = None
        start_time = time.time()
        try: 
            cmd_func(dict2namespace(args2))
        except Exception as e:
            exc = e
        finally:
            duration = time.time() - start_time
            end_time_str = time.strftime('%Y-%m-%d-%H-%M-%S')
            dump_file(join(step_dir, 'done.txt'), end_time_str)
            dump_file(join(step_dir, 'metrics', '_runtime.dat'), str(duration))
            if exc is not None:
                dump_file(join(step_dir, 'exception.txt'), str(exc))
                raise exc
            else:
                for arg, fname in hash_at_end.items():
                    dump_file(join(step_dir, 'outputs', arg+'.hash'), util.file.hash_file(fname, hash_algorithm='sha1'))
                    if fname.endswith('.metrics.tsv'):
                        shutil.copyfile(fname, join(step_dir, 'metrics', os.path.basename(fname)))
                dump_file(join(step_dir, 'success.txt'), end_time_str)

    return _run_cmd_with_tracking

class ProvenanceGraph(object):
    '''A graph representing a network of operations (steps), and their inputs and outputs.'''

    _instance = None

    @classmethod
    def instance(cls):
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance
    
    def __init__(self):
        self.graph = networkx.DiGraph()

    def load_provenance_data(path):
        '''Load raw provenance data into this graph.'''
        
        # The provenance data recorded by `add_provenance_tracking` can be in ad-hoc, rough or incomplete form.
        # Here we load all these pieces of data (the part 
        
