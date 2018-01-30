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

# * Prelims

__author__ = "ilya@broadinstitute.org"
__all__ = ["InFile", "OutFile", "InFiles", "OutFiles", "add_metadata_tracking", "is_metadata_tracking_enabled", "metadata_dir", 
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
import copy
import traceback
import collections
import functools
import itertools
import operator
import binascii
import fnmatch

# intra-module
import util.file
import util.misc
import util.version

# third-party
import networkx
import networkx.algorithms.dag

_log = logging.getLogger(__name__)

# * Recording of metadata

VIRAL_NGS_METADATA_FORMAT='1.0.0'

def metadata_dir():
    """Returns the directory to which metadata was recorded, as specified by the environment variable VIRAL_NGS_METADATA_PATH.
    Raises an error if the environment variable is not defined.
    """
    return os.environ['VIRAL_NGS_METADATA_PATH']

def is_metadata_tracking_enabled():
    return 'VIRAL_NGS_METADATA_PATH' in os.environ

    # check also that the only VIRAL_NGS_METADATA env vars are known ones


def _make_list(*x): return x

# ** class FileArg
class FileArg(object):

    '''The value of an argparse parser argument representing input or output file(s).  In addition to the string representing the
    argument value, keeps track of any filenames derived from the argument, and has methods for capturing metadata about the
    files to which they point.'''
    
    def __init__(self, val, mode, compute_fnames=_make_list):
        """Construct a FileArg.

        Args:
           val: the value of the command-line argument
           mode: 'r' if `val` points to input file(s), 'w' if to output files
           compute_fnames: function that will compute, from `val`, the list of actual filenames of the file(s) denoted by this 
             command-line argument.  By default, this is just one file and `val` contains its full name.  But `val` can be a 
             common prefix for a set of files with a given list of suffixes, or `val` can be a directory denoting all the files
             in the directory or just those matching a wildcard; and in those cases, compute_fnames will compute the actual file names
             by some non-trivial operation.
        """
        self.val, self.mode, self.compute_fnames = val, mode, compute_fnames

    def get_fnames(self):
        """Return the list of filename(s) specified by this command-line argument."""
        return self.compute_fnames(self.val)

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
                    file_info.update(inode=file_stat[stat.ST_INO], device=file_stat[stat.ST_DEV])
                except Exception:
                    _log.warning('Error getting file info for {} ({})'.format(file_arg, traceback.format_exc()))
            return file_info

        return dict(__FileArg__=True, val=self.val, mode=self.mode, files=list(map(file2dict, self.get_fnames())))

    @staticmethod
    def is_from_dict(val):
        return isinstance(val, collections.Mapping) and '__FileArg__' in val and isinstance(val['files'], list)

    def __str__(self):
        return '{}({})'.format('InFile' if self.mode=='r' else 'OutFile', self.val)

    def __repr__(self): return str(self)
        
# ** InFile, OutFile etc

def InFile(val, compute_fnames=_make_list):
    """Argparse argument type for arguments that denote input files."""
    file_arg = FileArg(val, mode='r', compute_fnames=compute_fnames)
    util.file.check_paths(read=file_arg.get_fnames())
    return file_arg

def OutFile(val, compute_fnames=_make_list):
    """Argparse argument type for arguments that denote output files."""
    file_arg = FileArg(val, mode='w', compute_fnames=compute_fnames)
    util.file.check_paths(write=file_arg.get_fnames())
    return file_arg

def InFiles(compute_fnames):
    """Argparse argument type for a string from which names of input files can be computed"""
    return functools.partial(InFile, compute_fnames=compute_fnames)

def OutFiles(compute_fnames):
    """Argparse argument type for a string from which names of output files can be computed"""
    return functools.partial(OutFile, compute_fnames=compute_fnames)


# ** Hashing of files

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

# ** Getting the environment
    
def create_run_id(t=None):
    """Generate a unique ID for a run (set of steps run as part of one workflow)."""
    return util.file.string_to_file_name('__'.join(map(str, (time.strftime('%Y%m%d%H%M%S', time.localtime(t))[2:], getpass.getuser(),
                                                             os.path.basename(os.getcwd()), uuid.uuid4()))))[:210]

def set_run_id():
    """Generate and record in the environment a unique ID for a run (set of steps run as part of one workflow)."""
    os.environ['VIRAL_NGS_METADATA_RUN_ID'] = create_run_id()

def _shell_cmd(cmd):
    """Run a command and return its output; if command fails, return empty string."""
    out = ''
    result = util.misc.run_and_print(cmd.strip().split(), silent=True)
    if result.returncode == 0:
        out = result.stdout
        if not isinstance(out, str):
            out = out.decode('utf-8')
        out = out.strip()
    return out
    
def tag_code_version(tag, push_to=None):
    """Create a lightweight git tag for the current state of the project repository, even if the state is dirty.
    If the repository is dirty, use the 'git stash create' command to create a commit representing the current state,
    and tag that; else, tag the existing clean state.  If `push_to` is not None, push the tag to the specified git remote.
    Return the git hash for the created git tag.  In case of any error, print a warning and return an empty string.
    """

    code_hash = ''

    try:
        with util.file.pushd_popd(util.version.get_project_path()):
            code_hash = _shell_cmd('git stash create') or _shell_cmd('git log -1 --format=%H')
            _shell_cmd('git tag ' + tag + ' ' + code_hash)
            if push_to:
                _shell_cmd('git push ' + push_to + ' ' + tag)
    except Exception:
        _log.warning('Could not create git tag: {}'.format(traceback.format_exc()))

    return code_hash

def get_conda_env():
    """Return the active conda environment"""
    return _shell_cmd('conda env export')

# ** add_metadata_tracking
def add_metadata_tracking(cmd_parser, cmd_main):
    """Add provenance tracking to the given command.  

    Called from util.cmd.attach_main().
    
    Args:
        cmd_parser: parser for a command defined in a script
        cmd_main: function implementing the command. Function takes one parameter, an argparse.Namespace, giving the values of the command's
             arguments.

    Returns:
        a wrapper for cmd_main, which has the same signature but adds metadata recording if enabled.
    """

    cmd_parser.add_argument('--metadata', action='append', nargs=2, metavar=('ATTRIBUTE', 'VALUE'),
                            help='attach metadata to this step (step=this specific execution of this command)')

    @functools.wraps(cmd_main)
    def _run_cmd_with_tracking(args):
        """Call the command implementation `cmd_main` with the arguments `args` parsed by `cmd_parser`, and record various
        metadata about the invocation."""

# *** Before calling cmd impl
        args_dict = vars(args).copy()

        # save any metadata specified on the command line.  then drop the 'metadata' argument from the args dict, since
        # the original command implementation `cmd_main` does not recognize this arg.
        metadata_from_cmd_line = { k[len('VIRAL_NGS_METADATA_VALUE_'):] : v
                                   for k, v in os.environ.items() if k.startswith('VIRAL_NGS_METADATA_VALUE_') }
        metadata_from_cmd_line.update(dict(args_dict.pop('metadata', {}) or {}))

        # for args denoting input or output files, for which 'type=InFile' or 'type=OutFile' was used when adding the args to
        # the parser, the corresponding values will be of type FileArg, rather than strings.  We must convert these values
        # to str before calling the original command implementation `cmd_main`.
        def replace_file_args(val):
            return (isinstance(val, FileArg) and val.val) or (isinstance(val, (list, tuple)) and list(map(replace_file_args, val))) or val

        args_new = argparse.Namespace(**{arg: replace_file_args(val) for arg, val in args_dict.items()})

        cmd_module=os.path.splitext(os.path.basename(sys.argv[0]))[0]
        cmd_name = args_dict.get('command', cmd_main.__name__)
        
        # Determine the run id and the step id for this step.  A step is a particular invocation of a command; a run is a set
        # of steps invoked as part of one workflow, such as one Snakemake invocation.
        # run_id is the same for all steps run as part of a single workflow.
        # if not given in the environment, create a run_id for a one-step workflow consisting of just this step.
        run_id = os.environ.get('VIRAL_NGS_METADATA_RUN_ID', create_run_id(beg_time))
        step_id = '__'.join(map(str, (create_run_id(beg_time), cmd_module, cmd_name)))

        # Sometimes, the implementation of a command will invoke another command as a subcommand.
        # We keep, in an environment variable, the list of any steps already running, and record this info as part of step metadata.
        save_steps_running = os.environ.get('VIRAL_NGS_METADATA_STEPS_RUNNING', '')
        os.environ['VIRAL_NGS_METADATA_STEPS_RUNNING'] = ((save_steps_running+':') if save_steps_running else '') + step_id

        cmd_exception, cmd_exception_str, cmd_result = None, None, None

        try:
            beg_time = time.time()

            # *** Run the actual command ***
            cmd_result = cmd_main(args_new)
        except Exception as e:
            cmd_exception = e
            cmd_exception_str = traceback.format_exc()
        finally:
            os.environ['VIRAL_NGS_METADATA_STEPS_RUNNING'] = save_steps_running
            try:  # if any errors happen during metadata recording just issue a warning
                # If command was cancelled by the user by Ctrl-C, skip the metadata recording; but if it failed with an exception,
                # still record that.
                if not isinstance(cmd_exception, KeyboardInterrupt):
# *** Record metadata after cmd impl returns
                    end_time = time.time()

                    _log.info('command {}.{} finished in {}s; exception={}'.format(cmd_module, cmd_name, end_time-beg_time, 
                                                                                   cmd_exception_str))
                    _log.info('recording metadata to {}'.format(metadata_dir()))

                    # record the code version used to run this step
                    code_repo = os.path.join(metadata_dir(), 'code_repo')
                    code_hash = tag_code_version('cmd_' + step_id, push_to=code_repo if os.path.isdir(code_repo) else None)

                    step_data = dict(__viral_ngs_metadata__=True, format=VIRAL_NGS_METADATA_FORMAT)

                    metadata_from_cmd_return = cmd_result if isinstance(cmd_result, collections.Mapping) and '__metadata__' in cmd_result \
                                               else {}

                    args_dict.pop('func_main', '')

                    hasher = Hasher()

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
                                                          cwd=os.getcwd(), conda_env=get_conda_env()),
                                             run_info=dict(beg_time=beg_time, end_time=end_time, duration=end_time-beg_time,
                                                           exception=cmd_exception_str,
                                                           argv=tuple(sys.argv)),
                                             args=args_dict,
                                             metadata_from_cmd_line=metadata_from_cmd_line,
                                             metadata_from_cmd_return=metadata_from_cmd_return,
                                             enclosing_steps=save_steps_running)

                    def write_obj(x):
                        return (isinstance(x, FileArg) and x.to_dict(hasher, cmd_exception is None)) or str(x)
                    
                    json_str = json.dumps(step_data, sort_keys=True, indent=4, default=write_obj)
                    crc32 = format(binascii.crc32(json_str.encode()) & 0xffffffff, '08x')

                    util.file.dump_file(os.path.join(metadata_dir(), '{}.crc32_{}.json'.format(step_id, crc32)), json_str)
                    _log.info('metadata recording took {}s'.format(time.time() - end_time))

            except Exception:
                # metadata recording is not an essential operation, so if anything goes wrong we just print a warning
                _log.warning('Error recording metadata ({})'.format(traceback.format_exc()))

        if cmd_exception:
            _log.warning('Command failed with exception: {}'.format(cmd_exception_str))
            raise cmd_exception

    return _run_cmd_with_tracking

def dict_has_keys(d, keys_str):
    """Test whether a `d` is a dict containing all the given keys (given as tokens of `keys_str`)"""
    return isinstance(d, collections.Mapping) and d.keys() >= set(keys_str.split())

def is_valid_step_record(d):
    """Test whether `d` is a dictionary containing all the expected elements of a step, as recorded by the code above"""
    return dict_has_keys(d, 'format step') and \
        dict_has_keys(d['step'], 'args step_id cmd_module')
#    return dict_has_keys(d, '__viral_ngs_metadata__ format step') and \
#        dict_has_keys(d['step'], 'args step_id cmd_module version_info run_env run_info')

##########################

# ** disabling metadata tracking    
###
#
# Section: disabling metadata tracking
#
# If metadata tracking is not enabled, disable the code in this module, in a simple way that ensures it won't affect normal operations.
# operations.  We do this by replacing externally called routines with simple stubs.
#
###

def _return_str(*args, **kw):
    return str

def _add_metadata_tracking_dummy(cmd_parser, cmd_main):
    """Add the --metadata command-line argument to `cmd_parser`, but make it a no-op; return `cmd_main` unchanged."""
    cmd_parser.add_argument('--metadata', nargs=2,
                            metavar=('ATTRIBUTE', 'VALUE'),
                            help='(DISABLED because metadata tracking is disabled, set environment variable '
                            'VIRAL_NGS_METADATA_PATH to enable) attach metadata to this step (step=this specific execution of this command)')
    @functools.wraps(cmd_main)
    def _run_cmd_with_tracking(args):
        delattr(args, 'metadata')
        return cmd_main(args)
    
    return _run_cmd_with_tracking

if not is_metadata_tracking_enabled():
    InFile, OutFile, InFiles, OutFiles = str, str, _return_str, _return_str
    add_metadata_tracking = _add_metadata_tracking_dummy

# * Analysis of metadata 
##########################

# ** class ProvenanceGraph
class ProvenanceGraph(networkx.DiGraph):

    '''Provenance graph representing a set of computations.  It has two types of nodes: steps and files.  A step node
    represents a computation step (a specific execution of a specific command); a file node represents a particular file
    that was read/written by one or more steps.

    Edges go from input files of a step to the step, and from a step to its output files, so the graph is bipartite.
    Various information about the parameters and execution of a step is represented as attributes of the step node.

    Fields:
       hash2files: map from file hash to file nodes with that hash
    '''

    FileNode = collections.namedtuple('FileNode', 'realpath file_hash mtime')

    def __init__(self):
        """Initialize an empty provenance graph."""
        super(ProvenanceGraph, self).__init__()

    @property
    def file_nodes(self): 
        """The nodes representing files"""
        return [n for n in self if self.nodes[n]['node_kind'] == 'file']

    @property
    def step_nodes(self): 
        """The nodes representing steps"""
        return [n for n in self if self.nodes[n]['node_kind'] == 'step']

# *** Loading pgraph
    def load(self, metadata_path=None, max_age_days=10000):
        """Read the provenance graph from the metadata directory.

        Args:
           path: path to the metadata directory; if None (default), use the path specified by the environment var VIRAL_NGS_METADATA_PATH.
           max_age_days: ignore steps older than this many days
        """

        metadata_path = metadata_path or metadata_dir()

# **** Load steps and files
        for step_record_fname in os.listdir(metadata_path or metadata_dir()):
            if not step_record_fname.endswith('.json'): continue

            _log.info('loading step {}'.format(step_record_fname))
            json_str = util.file.slurp_file(os.path.join(metadata_path, step_record_fname))
            #crc32 = format(binascii.crc32(json_str.encode()) & 0xffffffff, '08x')
            #assert step_record_fname.endswith('.{}.json'.format(crc32))

            step_record = json.loads(json_str)

            # check that all expected elements are in the json file -- some older-format files might be missing data
            if not is_valid_step_record(step_record): 
                _log.warning('not valid step: {}'.format(step_record_fname))
                continue
            if step_record['step']['run_info']['exception']: continue  # for now, ignore steps that failed
            if step_record['step'].get('enclosing_steps', ''): continue  # for now, skip steps that are sub-steps of other steps
            if ((time.time() - step_record['step']['run_info']['beg_time']) / (60*60*24)) > max_age_days: continue

            step_id = step_record['step']['step_id']
            self.add_node(step_id, node_kind='step', **step_record['step'])

            # Add input and output files of this step as data nodes.
            for arg, val in step_record['step']['args'].items():
                def gather_files(val):
                    """Return a flattened list of files denoted by this command argument (or empty list if it does not denote files)."""
                    return (FileArg.is_from_dict(val) and [val]) \
                        or (isinstance(val, (list, tuple)) and functools.reduce(operator.concat, list(map(gather_files, val)), [])) or []

                for files in gather_files(val):
                    assert FileArg.is_from_dict(files)
                    for f in files['files']:
                        if dict_has_keys(f, 'hash fname realpath size'):
                            file_node = self.FileNode(f['realpath'], f['hash'], f['mtime'])
                            self.add_node(file_node, node_kind='file', **f)
                            e = (file_node, step_id) if files['mode'] == 'r' else (step_id, file_node)
                            self.add_edge(*e, arg=arg)

                            # gather any per-file metadata specified on the command line.  currently, such metadata can be specified
                            # for a given command arg, and applied to all files denoted by the arg.
                            for metadata_attr, metadata_val in dict(step_record['step']['metadata_from_cmd_line']).items():
                                pfx = 'file.{}.'.format(arg)
                                if metadata_attr.startswith(pfx):
                                    self[e[0]][e[1]][metadata_attr[len(pfx):]] = metadata_val
                    # end: for f in files['files']
                # end: for files in gather_files(val)
            # end: for arg, val in step_record['step']['args'].items()
        # end: for step_record_fname in os.listdir(path)

# **** Check for anomalies

        for f in self.file_nodes:
            if self.in_degree[f] > 1:
                print('ANOMALY: file with indegree {}: {}'.format(self.in_degree[f], f))
                for e in self.in_edges(nbunch=f, data=True):
                    print(e[0], e[1], json.dumps(e[2], indent=4))

# **** Reconstruct missing connections between steps and files

        hash2files = collections.defaultdict(set)
        realpath2files = collections.defaultdict(set)
        for f in self.file_nodes:
            hash2files[f.file_hash].add(f)
            realpath2files[f.realpath].add(f)

        for f in self.file_nodes:
            if not self.pred[f]:
                # file f is read by some step, but we don't have a record of a step that wrote this exact file.
                # do we have 
                f2s = [ f2 for f2 in (hash2files[f.file_hash] & realpath2files[f.realpath]) if f2 != f and f2.mtime < f.mtime ]
                if f2s:
                    f2s = sorted(f2s, key=operator.attrgetter('mtime'))
                    for s in list(self.succ[f]):
                        f2s_bef = [ f2 for f2 in f2s if f2.mtime < self.nodes[s]['run_info']['beg_time'] ]
                        if f2s_bef:
                            f2 = f2s_bef[-1]
                            assert f2 != f and f2.realpath == f.realpath and f2.file_hash == f.file_hash and f2.mtime < f.mtime
                            self.add_edge(f2, s, **self.edges[f, s])
                            self.remove_edge(f, s)
                            _log.info('reconnected: {}->{}'.format(f2.realpath, s))

        # process hardlinks: if a file doesn't have a creation record, but is the same as a file that does, then
        # reroute input links to the file that does have a creation record.

        if False:

            new_edges = []
            del_nodes = []
            for file_node in self.file_nodes:
                #assert self.in_degree[file_node] <= 1
                if self.in_degree[file_node] == 0:
                    if not os.path.isfile(file_node.realpath): continue
                    # loop over files which have contents identical to file_node, but which have a creation record
                    # (i.e. are the output of a known step)
                    for file_node_2 in self.hash2files[file_node.file_hash]:
                        if file_node_2 != file_node and file_node_2.mtime == file_node.mtime and self.in_degree[file_node_2] > 0 and \
                           os.path.isfile(file_node_2.realpath) and os.path.samefile(file_node.realpath, file_node_2.realpath):
                            # merge file_node into file_node_2: for any steps that read file_node, reroute them to read file_node_2 instead.
                            for n, step_reading_n, edge_attrs in self.out_edges(file_node, data=True):
                                assert n == file_node
                                new_edges.append((file_node_2, step_reading_n, copy.deepcopy(edge_attrs)))
                                _log.info('rerouted input of {} from {} to {}'.format(step_reading_n, n, file_node_2))
                            del_nodes.append(n)
                            break
                    # end: for file_node_2 in self.hash2files[file_node.file_hash]:
                # end: if self.in_degree[file_node] == 0
            # end: for file_node in file_nodes
            self.add_edges_from(new_edges)
            self.remove_nodes_from(del_nodes)

        assert networkx.algorithms.dag.is_directed_acyclic_graph(self)

        #
        # Print nodes with unknown origin
        # 

        unknown_origin_files = []
        for file_node in self.file_nodes:
            #assert self.in_degree[file_node] <= 1
            if self.in_degree[file_node] > 1:
                _log.warning('indegree of {} is {}'.format(file_node, self.in_degree[file_node]))
            if self.in_degree[file_node] == 0 and os.path.isfile(file_node.realpath):
                unknown_origin_files.append(file_node.realpath)
                #_log.warning('UNKNOWN ORIGIN: {}'.format(file_node))

        print('UNKNOWN_ORIGIN:\n', '\n'.join(unknown_origin_files))

    # end: def load(self, path=None)

# *** write_dot    
    def write_dot(self, dotfile, nodes=None, ignore_cmds=(), ignore_exts=()):
        ignore_exts = ()

        def get_val(d, keys):
            """Fetch a value from a nested dict using a sequence of keys"""
            if not keys: return d
            keys = util.misc.make_seq(keys)
            if isinstance(d, collections.Mapping) and keys[0] in d:
                return get_val(d[keys[0]], keys[1:])
            return None

        node2id = {}

        def fix_name(s):
            """Return a string based on `s` but valid as a GraphViz node name"""
            return 'n' + str(node2id.setdefault(s, len(node2id)))
        # util.file.string_to_file_name(str(s)).replace('-', '_').replace('.', '_')

        ignored = set()

        with open(dotfile, 'wt') as out:
            out.write('digraph G {\n')
            for n in self:
                if not nodes or n in nodes:
                    n_attrs = self.nodes[n]
                    if self.nodes[n]['node_kind'] == 'step':
                        label = get_val(n_attrs, 'metadata_from_cmd_line step_name'.split()) or get_val(n_attrs, 'cmd_name') or 'unknown_cmd'
                        if label in ignore_cmds: 
                            ignored.add(n)
                            continue
                        shape = 'invhouse'
                    else:
                        label = get_val(n_attrs, 'fname') or 'noname'
                        if any(label.endswith(e) for e in ignore_exts): 
                            ignored.add(n)
                            continue
                        #if ',' in label:
                        #    label = ','.join(filter(lambda f: not f.startswith('/tmp'), label.split(',')))
                        shape = 'oval'
                        
                    out.write('{} [label="{}", shape={}];\n'.format(fix_name(n), label, shape))

            for u, v, arg in self.edges(data='arg'):
                if nodes and (u not in nodes or v not in nodes): continue
                if u in ignored or v in ignored: continue
                out.write('{} -> {} [label="{}"];\n'.format(fix_name(u), fix_name(v), arg))
            out.write('}\n')
    # end: def write_dot(self, dotfile, nodes=None, ignore_cmds=(), ignore_exts=()):
        
# end: class ProvenanceGraph(object)

# ** compute_paths
def compute_paths():
    """Compute computation paths, and their attributes."""

    G = ProvenanceGraph()
    G.load()

    beg_nodes = [f for f in G.file_nodes if fnmatch.fnmatch(f.realpath, '*/data/00_raw/*.bam')]
    end_nodes = [f for f in G.file_nodes if fnmatch.fnmatch(f.realpath, '*/data/02_assembly/*.fasta')]
    for i, e in enumerate(end_nodes):
        ancs = networkx.algorithms.dag.ancestors(G, e)
        ancs.add(e)
        if True: #set(beg_nodes) & ancs:

            print('===============================beg ', i)
            for a in ancs:
                if G.nodes[a]['node_kind'] == 'step':
                    print('-----------------------', json.dumps(G.nodes[a], indent=4))

            print('===============================end ', i)

            dot_fname = 'pgraph{:03}.dot'.format(i)
            svg_fname = 'pgraph{:03}.svg'.format(i)
            G.write_dot(dot_fname, nodes=ancs, ignore_cmds=['main_vcf_to_fasta'], 
                        ignore_exts=['.fai', '.dict', '.nix']+['.bitmask', '.nhr', '.nin', '.nsq']+
                        ['.srprism.'+ext for ext in 'amp idx imp map pmp rmp ss ssa ssd'.split()])
            _shell_cmd('dot -Tsvg -o {} {}'.format(svg_fname, dot_fname))
            _log.info('created {}'.format(svg_fname))
        

    # beg_edges = [e[:2] for e in g.out_edges(nbunch=data_nodes, data=True) if e[2]['ext']=='.bam' and 'data/00_raw' in e[2]['fname']]
    # end_edges = [e[:2] for e in g.in_edges(nbunch=data_nodes, data=True) if e[2].get('role', None) == 'final_assembly']

    # beg_nodes = set(map(operator.itemgetter(0), beg_edges))
    # end_nodes = set(map(operator.itemgetter(1), end_edges))

    # for end_edge in end_edges:
    #     #with networkx.utils.contextmanagers.reversed(g):
    #     #    pred = networkx.predecessor(g, end_node)
    #     end_node = end_edge[0]
    #     print('end_node=', end_node)
    #     prereqs = pgraph.find_prereq_steps(end_node)
    #     reachable_beg_nodes = prereqs# & beg_nodes
    #     print('end_node=', end_node, 'reachable_beg_nodes=({})'.format(len(reachable_beg_nodes)), 
    #           '\n'.join(map(str,['\nn={} \ninp={} \nout={}'.format(n, '\n'.join(g.nodes[n].get('inputs',[])), '\n'.join(g.nodes[n].get('outputs',[]))) for n in reachable_beg_nodes])))

    #     # for b in reachable_beg_nodes:
    #     #     print('===========path from', type(b), g.nodes[b])
    #     #     n = b
    #     #     while n != end_node:
    #     #         p = pred[n][0]
    #     #         print(p)
    #     #         print('     in=', g.nodes[p].get('inputs', []))
    #     #         print('    out=', g.nodes[p].get('outputs', []))
    #     #         print(' ')
    #     #         n = p
    

    #start_nodes = [ n for n in pgraph if dict().viewitems() <= g.nodes[n].viewitems() ]

################################################################    

def _setup_logger(log_level):
    loglevel = getattr(logging, log_level.upper(), None)
    assert loglevel, "unrecognized log level: %s" % log_level
    _log.setLevel(loglevel)
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s"))
    _log.addHandler(h)


if __name__ == '__main__':
    #import assembly
    _setup_logger('INFO')
    compute_paths()

    if False:
        pgraph = ProvenanceGraph()
        pgraph.load(max_age_days=100)
        dot_fname = 'pgraph.dot'
        svg_fname = 'pgraph.svg'
        pgraph.write_dot(dot_fname, ignore_cmds=['main_vcf_to_fasta'], 
                         ignore_exts=['.fai', '.dict', '.nix']+['.bitmask', '.nhr', '.nin', '.nsq']+
                         ['.srprism.'+ext for ext in 'amp idx imp map pmp rmp ss ssa ssd'.split()])
        _shell_cmd('dot -Tsvg -o {} {}'.format(svg_fname, dot_fname))


    if False:
        for i, C in enumerate(networkx.connected_components(pgraph.pgraph.to_undirected())):
            dot_fname = 'pgraph{:03}.dot'.format(i)
            svg_fname = 'pgraph{:03}.svg'.format(i)
            pgraph.write_dot(dot_fname, nodes=C, ignore_cmds=['main_vcf_to_fasta'], 
                             ignore_exts=['.fai', '.dict', '.nix']+['.bitmask', '.nhr', '.nin', '.nsq']+
                             ['.srprism.'+ext for ext in 'amp idx imp map pmp rmp ss ssa ssd'.split()])
            _shell_cmd('dot -Tsvg -o {} {}'.format(svg_fname, dot_fname))
            _log.info('created {}'.format(svg_fname))


