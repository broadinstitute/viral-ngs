"""Caching and reuse of computation results"""

import os
import os.path
import argparse
import functools
import collections
import shutil
import json
import contextlib

import util.file

from . import metadata_dir
from .file_arg import FileArg
from util._metadata import recording
from .hashing import Hasher
from .testmon import testmon_core

import fs

# * Recording which tests run which commands

test2cmds = collections.defaultdict(set)
test_running = None

def instrument_module(module_globals):
    if 'VIRAL_NGS_GATHER_CMDS_IN_TESTS' not in os.environ: return
    for cmd_name, cmd_parser in module_globals['__commands__']:
        p = argparse.ArgumentParser()
        cmd_parser(p)
        cmd_main = p.get_default('func_main')
        if hasattr(cmd_main, '_cmd_main_orig'):
            cmd_main = cmd_main._cmd_main_orig
        print('cmd_name=', cmd_name, 'cmd_main=', cmd_main, 'module=', cmd_main.__module__, 'name', cmd_main.__name__)

        def make_caller(cmd_name, cmd_main):

            @functools.wraps(cmd_main)
            def my_cmd_main(*args, **kw):
                print('\n***FUNC IMPLEMENTING CMD ', module_globals['__name__'], cmd_name)
                test2cmds[test_running].add((module_globals['__name__'], cmd_name))
                return cmd_main(*args, **kw)

            return my_cmd_main

        my_cmd_main = make_caller(cmd_name, cmd_main)

        assert module_globals[cmd_main.__name__] == cmd_main
        module_globals[cmd_main.__name__] = my_cmd_main

def record_test_start(nodeid):
    global test_running
    test_running = nodeid

def tests_ended():
    #print('************GATHERED:\n{}'.format('\n'.join(map(str, test2cmds.items()))))
    pass

# * Caching of results

class FileCache(object):
    """Manages a cache of data files produced by commands."""

    def __init__(self, cache_dir):
        self.cache_dir = cache_dir

    def save_file(self, fname, file_hash):
        """Save the given file to the cache."""
        if not os.path.isfile(os.path.join(self.cache_dir, file_hash)):
            # FIXME: race condition
            shutil.copyfile(fname, os.path.join(self.cache_dir, file_hash))

    def has_file_with_hash(self, file_hash):
        return os.path.isfile(os.path.join(self.cache_dir, file_hash))

    def fetch_file(self, file_hash, dest_fname):
        shutil.copyfile(os.path.join(self.cache_dir, file_hash), dest_fname)

def cache_results(file_info):
    """If caching of results is on, and file_info contains output file(s), save them in the cache."""
    if 'VIRAL_NGS_DATA_CACHE' not in os.environ: return
    if file_info['mode'] == 'r': return
    cache = FileCache(os.environ['VIRAL_NGS_DATA_CACHE'])
    for f in file_info['files']:
        if 'hash' in f:
            cache.save_file(f['abspath'], f['hash'])

def reuse_cached_step(cmd_module, cmd_name, args):
    """If this step has been run with the same args before, and we have saved the results, reuse the results."""
    if 'VIRAL_NGS_DATA_CACHE' not in os.environ: return

    def replace_file_args(val):
        if isinstance(val, FileArg):
            if val.mode == 'w': return '_out_file_arg'
            file_info = val.gather_file_info(hasher=Hasher(), out_files_exist=False)
            return [f['hash'] for f in file_info['files']]
        if FileArg.is_from_dict(val):
            if val['mode'] == 'w': return '_out_file_arg'
            return [f['hash'] for f in val['files']]

        if isinstance(val, (list, tuple)): return list(map(replace_file_args, val))
        return val

    cur_args = {arg: replace_file_args(val) for arg, val in args.items() if arg != 'func_main'}

    def flatten_file_args(val):
        if isinstance(val, FileArg):
            if val.mode == 'r': return []
            return list(val.fnames)
        if FileArg.is_from_dict(val):
            if val['mode'] == 'r': return []
            return [f['hash'] for f in val['files']]

        if isinstance(val, (list, tuple)): return functools.reduce(operator.concat, list(map(flatten_file_args, val)), [])
        return []

    flat_cur_args = {arg: flatten_file_args(val) for arg, val in args.items() if arg != 'func_main'}

    with fs.open_fs(metadata_dir()) as metadata_fs:
        for step_record_fname in metadata_fs.listdir('/'):
            if step_record_fname.endswith('.json') and cmd_name in step_record_fname:
                json_str = util.file.slurp_file(os.path.join(metadata_dir(), step_record_fname))
                step_record = json.loads(json_str)
                if not recording.is_valid_step_record(step_record): continue
                if step_record['step']['run_info']['exception']: continue  # for now, ignore steps that failed
                if step_record['step'].get('enclosing_steps', ''): continue  # for now, skip steps that are sub-steps of other steps
                if step_record['step']['cmd_module'] != cmd_module or step_record['step']['cmd_name'] != cmd_name: continue
                cached_args = {arg: replace_file_args(val) for arg, val in step_record['step']['args'].items()}
                if cached_args == cur_args:
                    if os.path.isfile(os.path.join(metadata_dir(), step_record_fname+'.testmondata')):
                        print('testing for reuse: {}'.format(os.path.join(metadata_dir(), step_record_fname+'.testmondata')))

                        with contextlib.closing(testmon_core.TestmonData(os.path.realpath(util.version.get_project_path()), 
                                                                         datafile=os.path.join(metadata_dir(), 
                                                                                               step_record_fname+'.testmondata'))) as testmon_data:

                            testmon_data.read_data()
                            testmon_data.read_source()
                            if not testmon_data.test_should_run(step_record['step']['cmd_module']+'::'+step_record['step']['cmd_name']):
                                # test that the cache has all the output files

                                print('CAN REUSE! {}'.format(step_record_fname))
                                flat_cached_args = {arg: flatten_file_args(val) for arg, val in step_record['step']['args'].items()}
                                all_found=True
                                cache = FileCache(os.environ['VIRAL_NGS_DATA_CACHE'])

                                for arg in flat_cur_args:
                                    for f_hash, f in zip(flat_cached_args[arg], flat_cur_args[arg]):
                                        if not cache.has_file_with_hash(f_hash):
                                            print('not cached: {}'.format(f))
                                            all_found = False
                                        else:
                                            print('cached! {}'.format(f))

                                if all_found:

                                    for arg in flat_cur_args:
                                        for f_hash, f in zip(flat_cached_args[arg], flat_cur_args[arg]):
                                            cache.fetch_file(f_hash, f)
                                            print('fetched from cache: {}'.format(f))

                                    return step_record_fname

# end: def reuse_cached_step(cmd_module, cmd_name, args):
