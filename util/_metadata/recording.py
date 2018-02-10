"""Recording of metadata"""

# built-ins
import os
import os.path
import stat
import sys
import argparse
import getpass
import socket
import collections
import traceback
import functools
import time
import uuid
import binascii
import json
import platform

# intra-module
import util.file
import util.misc
import util.version

from .file_arg import FileArg
from .hashing import Hasher
from util._metadata import caching
from .md_utils import _make_list, _shell_cmd, _mask_secret_info, dict_has_keys
from . import _log, metadata_dir, is_metadata_tracking_enabled

# third-party
import fs
import fs.path
import fs_s3fs

VIRAL_NGS_METADATA_FORMAT='1.0.0'

# ** run_id management

def create_run_id(t=None):
    """Generate a unique ID for a run (set of steps run as part of one workflow)."""
    return util.file.string_to_file_name('__'.join(map(str, (time.strftime('%Y%m%d%H%M%S', time.localtime(t))[2:], getpass.getuser(),
                                                             os.path.basename(os.getcwd()), uuid.uuid4()))))[:210]

def set_run_id():
    """Generate and record in the environment a unique ID for a run (set of steps run as part of one workflow)."""
    os.environ['VIRAL_NGS_METADATA_RUN_ID'] = create_run_id()

# ** Getting the execution environment

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

def add_metadata_arg(cmd_parser, help_extra=''):
    """Add --metadata arg to `cmd_parser`"""
    if not getattr(cmd_parser, 'metadata_arg_added', False):
        cmd_parser.add_argument('--metadata', nargs=2, metavar=('ATTRIBUTE', 'VALUE'), action='append',
                                help='attach metadata to this step (step=this specific execution of this command)' + help_extra)
        setattr(cmd_parser, 'metadata_arg_added', True)

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
    add_metadata_arg(cmd_parser)

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
            if isinstance(val, FileArg): return val.val
            if isinstance(val, (list, tuple)): return list(map(replace_file_args, val))
            return val

        args_new = argparse.Namespace(**{arg: replace_file_args(val) for arg, val in args_dict.items()})

        cmd_module=os.path.splitext(os.path.basename(sys.argv[0]))[0]
        cmd_name = args_dict.get('command', cmd_main.__name__)
        
        # Determine the run id and the step id for this step.  A step is a particular invocation of a command; a run is a set
        # of steps invoked as part of one workflow, such as one Snakemake invocation.
        # run_id is the same for all steps run as part of a single workflow.
        # if not given in the environment, create a run_id for a one-step workflow consisting of just this step.
        beg_time = time.time()
        run_id = os.environ.get('VIRAL_NGS_METADATA_RUN_ID', create_run_id(beg_time))
        step_id = '__'.join(map(str, (create_run_id(beg_time), cmd_module, cmd_name)))

        # Sometimes, the implementation of a command will invoke another command as a subcommand.
        # We keep, in an environment variable, the list of any steps already running, and record this info as part of step metadata.
        save_steps_running = os.environ.get('VIRAL_NGS_METADATA_STEPS_RUNNING', '')
        os.environ['VIRAL_NGS_METADATA_STEPS_RUNNING'] = ((save_steps_running+':') if save_steps_running else '') + step_id

        caching.reuse_cached_step(cmd_module, cmd_name, args_dict)

        cmd_exception, cmd_exception_str, cmd_result = None, None, None

        if 'VIRAL_NGS_TESTMON' in os.environ:
            testmon_data = testmon_core.TestmonData(os.path.realpath(util.version.get_project_path()))
            testmon_data.read_data()
            testmon_data.read_source()
            testmon = testmon_core.Testmon([util.version.get_project_path()], set(['singleprocess']))
            testmon.start()

        try:
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
                if is_metadata_tracking_enabled() and not isinstance(cmd_exception, KeyboardInterrupt):

                    if 'VIRAL_NGS_TESTMON' in os.environ:
                        if cmd_exception: testmon.stop()
                        else:
                            print('saving testmon')
                            testmon.stop_and_save(testmon_data, util.version.get_project_path(), cmd_module+'::'+cmd_name,
                                                  [dict(outcome='ok')])
                            print('writing data')
                            testmon_data.write_data()
                            testmon.close()
                            testmon_data.close_connection()


# *** Record metadata after cmd impl returns
                    end_time = time.time()

                    _log.info('command {}.{} finished in {}s; exception={}'.format(cmd_module, cmd_name, end_time-beg_time, 
                                                                                   cmd_exception_str))
                    _log.info('recording metadata to {}'.format(_mask_secret_info(metadata_dir())))

                    # record the code version used to run this step
                    code_repo = os.path.join(metadata_dir(), 'code_repo')
                    code_hash = tag_code_version('cmd_' + step_id, push_to=code_repo if os.path.isdir(code_repo) else None)

                    # The function that implements the command can pass us some metadata to be included in the step record,
                    # by returning a mapping with '__metadata__' as one key.  The remaining key-value pairs of the mapping are thenn
                    # treated as metadata.
                    metadata_from_cmd_return = cmd_result if isinstance(cmd_result, collections.Mapping) and '__metadata__' in cmd_result \
                                               else {}

                    args_dict.pop('func_main', '')  # 

                    step_data = dict(__viral_ngs_metadata__=True, format=VIRAL_NGS_METADATA_FORMAT)
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
                    #
                    # Serialize the record of this step to json.  In the process, for any FileArg args of the command,
                    # gather hashsums and other file info for the denoted file(s).
                    #

                    hasher = Hasher()

                    def write_obj(x):
                        """If `x` is a FileArg, return a dict representing it, else return a string representation of `x`.
                        Used for json serialization below."""
                        if not isinstance(x, FileArg): return str(x)
                        file_info = x.gather_file_info(hasher, out_files_exist=cmd_exception is None)
                        caching.cache_results(file_info)
                        return file_info
                    
                    json_str = json.dumps(step_data, sort_keys=True, indent=4, default=write_obj)

                    # as a sanity check, we compute the CRC of the json file contents, and make that part of the filename.
                    crc32 = format(binascii.crc32(json_str.encode()) & 0xffffffff, '08x')
                    json_fname = '{}.crc32_{}.json'.format(step_id, crc32)
                    
                    with fs.open_fs(metadata_dir()) as metadata_fs:
                        metadata_fs.settext(json_fname, json_str)

                    _log.info('metadata recording took {}s'.format(time.time() - end_time))

            except Exception:
                # metadata recording is not an essential operation, so if anything goes wrong we just print a warning
                _log.warning('Error recording metadata ({})'.format(traceback.format_exc()))

        if cmd_exception:
            _log.warning('Command failed with exception: {}'.format(cmd_exception_str))
            raise cmd_exception

    return _run_cmd_with_tracking

def is_valid_step_record(d):
    """Test whether `d` is a dictionary containing all the expected elements of a step, as recorded by the code above"""
    return dict_has_keys(d, 'format step') and \
        dict_has_keys(d['step'], 'args step_id cmd_module')
