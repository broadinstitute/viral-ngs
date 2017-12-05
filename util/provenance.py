'''Provenance tracking: for each file, automatically keeping track of how it was created -- from which other files and by which 
code version.  Also facilitates uniform representation of operation parameters and result metrics, enabling various queries.'''

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "OutMetricsFile", "add_provenance_tracking"]

# built-ins
import argparse
import logging
import os
import os.path
import platform
import shutil
import sys
import time
import uuid
import socket
import getpass
import json
import traceback

# intra-module
import util.file
import util.misc
import util.version

# third-party
import networkx
import networkx.readwrite.json_graph
import networkx.drawing.nx_pydot
import concurrent.futures

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

class OutMetricsFile(OutFile):

    def __new__(cls, *args, **kw):
        return OutFile.__new__(cls, *args, **kw)
    
def is_provenance_tracking_enabled():
    return 'VIRAL_NGS_PROV' in os.environ


class Hasher(object):
    """Manages computation of file hashes.

    It might also cache the actual contents of files, depending on size.
    """

    def __init__(self):
        pass

    def __call__(self, file):
        return util.file.hash_file(file, hash_algorithm='sha1')

def add_provenance_tracking(cmd_parser, cmd_func):
    """Add provenance tracking to the given command.
    
    Args:
        cmd_parser: parser for a command defined in a script
        cmd_func: function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.

    Returns:
        a wrapper for cmd_func, which has the same signature but adds provenance tracking (if provenance tracking is configured)
    """

    def _run_cmd_with_tracking(args):

        exception = None

        try:
            beg_time = time.time()
            print('calling', cmd_func, 'at time', beg_time)
            cmd_func(args)
            print('finished', cmd_func, 'at time', beg_time)
        except Exception as e:
            print('saw exception', str(e))
            traceback.print_tb(e)
            exception = e
        finally:
            if is_provenance_tracking_enabled():
                end_time = time.time()

                step_id = '-'.join(map(str, (time.strftime('%Y-%m-%d-%H-%M-%S', time.localtime(beg_time)), 
                                             cmd_func.__module__, cmd_func.__name__, uuid.uuid4())))

                pgraph = networkx.DiGraph()
                hasher = Hasher()

                pgraph.add_node(step_id, beg_time=beg_time, end_time=end_time, duration=end_time-beg_time,
                                exception=str(exception),
                                viral_ngs_version=util.version.get_version(),
                                platform=platform.platform(), cpus=util.misc.available_cpu_count(), host=socket.getfqdn(), 
                                user=getpass.getuser(),
                                argv=tuple(sys.argv))

                for arg, val in vars(args).items():
                    for i, v in enumerate(util.misc.make_seq(val)):
                        if v and isinstance(v, FileArg):
                            file_hash = hasher(v)
                            edge_attrs = dict(arg=arg)
                            if len(val) > 1:
                                edge_attrs.update(arg_order=i)
                            if isinstance(v, InFile):
                                pgraph.add_edge(file_hash, step_id, **edge_attrs)
                            elif exception is None:
                                assert isinstance(v, OutFile)
                                pgraph.add_edge(step_id, file_hash, **edge_attrs)

                util.file.dump_file(os.path.join(os.environ['VIRAL_NGS_PROV'], step_id+'.graphml'),
                                    networkx.readwrite.json_graph.jit_data(pgraph))
                networkx.drawing.nx_pydot.write_dot(pgraph, os.path.join(os.environ['VIRAL_NGS_PROV'], step_id+'.dot'))

    return _run_cmd_with_tracking

if __name__ == '__main__':
    print(bool(OutFile(None)))

