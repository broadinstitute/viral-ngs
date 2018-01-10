'''Automated recording of data provenance and metrics.  This lets you answer questions like: how was a given data file made?  Which workflow
versions tend to produce the best results, according to given metrics?

This module deals with data recording; ../data_utils.py deals with data querying.

The unit of recording in this module is one command (see cmd.py), as opposed to one Python function or one whole workflow.
To use this module, import InFile and OutFile from it; then, when defining argparse arguments for a command, in add_argument()
use "type=InFile" for input files or "type=OutFile" for output files. Until provenance tracking is enabled (see below how),
InFile and OutFile are defined to be simply str.  To enable provenance tracking, set the environment variable
VIRAL_NGS_METADATA_PATH to a writable directory.  Then, whenever a command is run, a file is written to that directory
recording the command, all its parameters, the input and output files, and details of the run environment.
Input and output files are identified by a hash of their contents, so that copied/moved/renamed files can be properly identified.
For each file, we also record metadata such as name, size, and modification date.


: for each data file, automatically keeping track of how it was created -- from which inputs, by which 
code version and with what parameters.


Note that only files actually used are listed in the metadata; optional files not specified in a given command invocation are not.


See data_utils.py for utils that create useful reports from provenance data.

Environment variables used:

   VIRAL_NGS_METADATA_PATH: location to which metadata should be recorded.  If this environment variable is not set, metadata recording
       is disabled, and this module has no effect.

Implementation notes:

Metadata recording is done on a best-effort basis.  If metadata recording fails for any reason, a warning is logged, but no error is raised.
'''

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "FilePrefix", "add_metadata_tracking", "is_metadata_tracking_enabled", "metadata_dir", 
           "load_provenance_graph"]

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
import functools
import itertools
import operator

# intra-module
import util.file
import util.misc
import util.version


_log = logging.getLogger(__name__)

#########################
# Recording of metadata #
#########################

VIRAL_NGS_METADATA_FORMAT='1.0.0'

def metadata_dir():
    """Returns the directory to which metadata was recorded, as specified by the environment variable VIRAL_NGS_METADATA_PATH.
    Raises an error if the environment variable is not defined.
    """
    return os.environ['VIRAL_NGS_METADATA_PATH']

def is_metadata_tracking_enabled():
    return 'VIRAL_NGS_METADATA_PATH' in os.environ

    # check also that the only VIRAL_NGS_METADATA env vars are known ones


def _make_list(x):
    return [x]

class FileArg(str):

    '''The value of an argparse parser argument representing input or output file(s).  In addition to the string representing the
    argument value, keeps track of any filenames derived from the argument, and has methods for capturing metadata about the
    files to which they point.'''
    
    def __new__(cls, val, mode, val2fnames=_make_list):
        """Construct a FileArg.

        Args:
           val: the value of the command-line argument
           mode: 'r' if `val` points to input file(s), 'w' if to output files
           val2fnames: function that will compute, from `val, the list of actual filenames of the file(s) denoted by this 
             command-line argument.  By default, this is just one file and `val` contains its full name.  But `val` can be a 
             common prefix for a set of files with a given list of suffixes, or `val` can be a directory denoting all the files
             in the directory or just those matching a wildcard; and in those cases, val2fnames will compute the actual file names.
        """
        s = str.__new__(cls, val)
        s.val, s.mode, s.val2fnames = val, mode, val2fnames
        return s

    def get_fnames(self):
        """Return the list of filename(s) specified by this argument."""
        return self.val2fnames(self.val)

    def to_dict(self, hasher, out_files_exist):
        """Return a dict representing metadata about the file(s) denoted by this argument.

        Args:
            hasher: callable for computing the hash value of a file
            out_files_exist: if False, don't expect output files to exist (because the command raised an exception)
        """

        def file2dict(file_arg):
            """Compute a dictionary of info about one file"""
            file_info = dict(fname=file_arg, realpath=os.path.realpath(file_arg), abspath=os.path.abspath(file_arg))
            if self.mode=='r' or out_files_exist:
                file_info.update(hash=hasher(file_arg))

                try:
                    file_stat = os.stat(file_arg)
                    file_info.update(size=file_stat[stat.ST_SIZE],
                                     mtime=file_stat[stat.ST_MTIME], ctime=file_stat[stat.ST_CTIME])
                    file_info.update(owner=pwd.getpwuid(file_stat[stat.ST_UID]).pw_name)
                except Exception:
                    _log.warning('Error getting file info for {} ({})'.format(file_arg, traceback.format_exc()))
            return file_info

        return dict(__FileArg__=True, val=self.val, mode=self.mode, files=map(file2dict, self.get_fnames()))

def InFile(val, *args, **kwargs):
    """Argparse argument type for arguments that denote input files."""
    file_arg = FileArg(val, mode='r', *args, **kwargs)
    util.file.check_paths(read=file_arg.get_fnames())
    return file_arg

def OutFile(val, *args, **kwargs):
    """Argparse argument type for arguments that denote output files."""
    file_arg = FileArg(val, mode='w', *args, **kwargs)
    util.file.check_paths(write=file_arg.get_fnames())
    return file_arg

def _add_suffixes(val, suffixes):
    return [val+sfx for sfx in suffixes]

def FilePrefix(InFile_or_OutFile, suffixes):
    """Argparse argument type for arguments that denote a prefix for a set of input or output files with known extensions.

    Usage examples: 

        # specify a base name for BLAST database files
        parser.add_argument(--blastDbs, type=FilePrefix(InFile, suffixes=('.nhr', '.nin', '.nsq')))
        # specify a FASTA file for which we also create an index
        parser.add_argument(outFasta, type=FilePrefix(OutFile, suffixes=('', '.fai')))
    """
    return functools.partial(InFile_or_OutFile, val2fnames=functools.partial(_add_suffixes, suffixes=suffixes))

class Hasher(object):
    """Manages computation of file hashes.

    It might also cache the actual contents of files, depending on size.
    """

    def __init__(self, hash_algorithm='sha1'):
        self.hash_algorithm = hash_algorithm

    def __call__(self, file):
        file_hash = ''
        try:
            if os.path.isfile(file) and not stat.S_ISFIFO(os.stat(file).st_mode):
                file_hash = self.hash_algorithm + '_' + util.file.hash_file(file, hash_algorithm=self.hash_algorithm)
        except Exception:
            _log.warning('Cannot compute hash for {}: {}'.format(file, traceback.format_exc()))
        return file_hash


def create_run_id(t=None):
    """Generate a unique ID for a run (set of steps run as part of one workflow)."""
    return util.file.string_to_file_name('__'.join(map(str, (time.strftime('%Y%m%d%H%M%S', time.localtime(t))[2:], getpass.getuser(),
                                                             os.path.basename(os.getcwd()), uuid.uuid4()))))[:210]

def set_run_id():
    """Generate and record in the environment a unique ID for a run (set of steps run as part of one workflow)."""
    os.environ['VIRAL_NGS_METADATA_RUN_ID'] = create_run_id()

def tag_code_version(tag, push_to=None):
    """Create a lightweight git tag for the current state of the project repository, even if the state is dirty.
    If the repository is dirty, use the 'git stash create' command to create a commit representing the current state,
    and tag that; else, tag the existing clean state.  If `push_to` is not None, push the tag to the specified git remote.
    Return the git hash for the created git tag.  In case of any error, print a warning and return an empty string.
    """

    def run_cmd(cmd):
        """Run a command and return its output"""
        out = ''
        result = util.misc.run_and_print(cmd.strip().split(), silent=True)
        if result.returncode == 0:
            out = result.stdout
            if not isinstance(out, str):
                out = out.decode('utf-8')
            out = out.strip()
        return out

    code_hash = ''

    try:
        with util.file.pushd_popd(util.version.get_project_path()):
            code_hash = run_cmd('git stash create') or run_cmd('git log -1 --format="%H"')
            run_cmd('git tag ' + tag + ' ' + code_hash)
            if push_to:
                run_cmd('git push ' + push_to + ' ' + tag)
    except Exception:
        _log.warning('Could not create git tag: {}'.format(traceback.format_exc()))

    return code_hash

def add_metadata_tracking(cmd_parser, cmd_main, cmd_main_orig):
    """Add provenance tracking to the given command.  

    Called from util.cmd.attach_main().
    
    Args:
        cmd_parser: parser for a command defined in a script
        cmd_main: function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.
        cmd_main_orig: the original cmd_main function defined in the command script, if cmd_main is a wrapper around the original 
             function, else same as cmd_main.

    Returns:
        a wrapper for cmd_main, which has the same signature but adds metadata recording if enabled.
    """

    cmd_parser.add_argument('--metadata', action='append', nargs=2, metavar=('ATTRIBUTE', 'VALUE'),
                            help='attach metadata to this step (step=this specific execution of this command)')

    cmd_module=os.path.splitext(os.path.basename(inspect.getfile(cmd_main_orig)))[0]
    cmd_module2=os.path.splitext(os.path.basename(sys.argv[0]))[0]
    assert cmd_module == cmd_module2

    def _run_cmd_with_tracking(args):

        cmd_name = args.command
        args_dict = vars(args).copy()
        delattr(args, 'metadata')

        cmd_exception, cmd_exception_str, cmd_result = (None,)*3

        try:
            _log.info('calling command {}.{}; metadata tracking: {}'.format(cmd_module, cmd_name, is_metadata_tracking_enabled()))
            beg_time = time.time()

            # *** Run the actual command ***
            cmd_result = cmd_main(args)
        except Exception as e:
            cmd_exception = e
            cmd_exception_str = traceback.format_exc()
        finally:
            try:  # if any errors happen during metadata recording just issue a warning
                if is_metadata_tracking_enabled():
                    end_time = time.time()
                    _log.info('command {}.{} finished in {}s; recording metadata to {}'.format(cmd_module, cmd_name, end_time-beg_time,
                                                                                               metadata_dir()))

                    #
                    # Record data pertaining to the whole step
                    # 

                    # run_id is the same for all steps run as part of a single workflow.
                    # if not given in the environment, create a run_id for a one-step workflow consisting of just this step.
                    run_id = os.environ.get('VIRAL_NGS_METADATA_RUN_ID', create_run_id(beg_time))
                    step_id = '__'.join(map(str, (create_run_id(beg_time), cmd_module, cmd_name)))

                    # record the code version used to run this step
                    code_repo = os.path.join(metadata_dir(), 'code_repo')
                    code_hash = tag_code_version('cmd_' + step_id, push_to=code_repo if os.path.isdir(code_repo) else None)

                    step_data = dict(format=VIRAL_NGS_METADATA_FORMAT)

                    metadata_from_cmd_line = args_dict.pop('metadata', {}) or {}
                    metadata_from_cmd_return = cmd_result if isinstance(cmd_result, collections.Mapping) and '__metadata__' in cmd_result \
                                               else {}

                    args_dict.pop('func_main', '')

                    hasher = Hasher()

                    def transform_val(val):
                        if isinstance(val, FileArg): return val.to_dict(hasher, cmd_exception is None)
                        if isinstance(val, list): return map(transform_val, val)
                        return val

                    args_dict = dict((k, transform_val(v)) for k, v in args_dict.items())

                    step_data['step'] = dict(step_id=step_id, run_id=run_id,
                                             cmd_module=cmd_module, cmd_name=cmd_name,
                                             version_info=dict(viral_ngs_version=util.version.get_version(),
                                                               viral_ngs_path=util.version.get_project_path(),
                                                               viral_ngs_path_real=os.path.realpath(util.version.get_project_path()),
                                                               code_hash=code_hash),
                                             run_env=dict(metadata_dir=metadata_dir(),
                                                          platform=platform.platform(), 
                                                          cpus=util.misc.available_cpu_count(), host=socket.getfqdn(),
                                                          user=getpass.getuser(),
                                                          cwd=os.getcwd()),
                                             run_info=dict(beg_time=beg_time, end_time=end_time, duration=end_time-beg_time,
                                                           exception=cmd_exception_str,
                                                           argv=tuple(sys.argv)),
                                             args={k: transform_val(v) for k, v in args_dict.items()},
                                             metadata_from_cmd_line=metadata_from_cmd_line,
                                             metadata_from_cmd_return=metadata_from_cmd_return)

                    util.file.dump_file(os.path.join(metadata_dir(), step_id+'.json'),
                                        json.dumps(step_data, sort_keys=True, indent=4, default=str))

            except Exception:
                # metadata recording is not an essential operation, so if anything goes wrong we just print a warning
                _log.warning('Error recording metadata ({})'.format(traceback.format_exc()))

        if cmd_exception:
            raise cmd_exception

    return _run_cmd_with_tracking

########################
# Analysis of metadata #
########################

class ProvenanceGraph(object):

    '''Provenance graph representing a set of computations.  It has two types of nodes: steps and files.  A step node represents
    a computation step (a specific execution of a specific command); a file node represents a data file.  Edges go from input files
    of a step to the step, and from a step to its output files.  Various information about the parameters and execution of a step
    is represented as attributes of the step node.'''

    def __init__(self):
        self.pgraph = networkx.DiGraph()

    def load(path=None):
        """Read the provenance graph."""
        
        if path is None: path = metadata_dir()

        for step_record_fname in os.listdir(path):
            _log.info('loading step {}'.format(step_record_fname))
            if step_record_fname.endswith('.json'):
                step_record = json.loads(util.file.slurp_file(os.path.join(path, step_record_fname)))
                step_id = step_record['step']['step_id']
                pgraph.add_node(step_id, node_kind='step', **step_record['step'])
                for file_kind in ('in_files', 'out_files'):
                    for port_name, file_info in step_record[file_kind].items():
                        port_id = step_id + '__' + port_name
                        pgraph.add_node(port_id, **file_info)
                        pgraph.add_node(port_id, node_kind=('input' if file_kind=='in_files' else 'output'))
                        pgraph.add_edge(*((port_id, step_id) if file_kind=='in_files' else (step_id, port_id)))
                        if 'hash' in file_info:
                            file_id = file_info['hash']
                            if file_id:
                                pgraph.add_node(file_id, node_kind='file')
                                pgraph.add_edge(*((file_id, port_id) if file_kind=='in_files' else (port_id, file_id)))

        return pgraph


def compute_paths(paths_file):
    """Compute computation paths, and their attributes."""

    g = load_provenance_graph()

    #start_nodes = [ n for n in pgraph if dict().viewitems() <= g.nodes[n].viewitems() ]


if __name__ == '__main__':
    z=OutFile('hi')
    print(z,type(z), z.suffixes, isinstance(z,str), str(z))
