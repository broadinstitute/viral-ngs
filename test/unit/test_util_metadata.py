"""Unit tests for util.metadata.

Integration tests are implemented as a pytest plugin in tests/conftest.py , which causes the entire test suite to be run under metadata
recording.
"""

__author__ = "ilya@broadinstitute.org"

import os
import os.path
import collections
import argparse
import json
import shutil
import warnings
import subprocess
import contextlib
import functools

import util.cmd
import util.file
import util.misc
from util._metadata import metadata_db, file_arg, md_utils
from test import tst_inp, TestCaseWithTmp

from test.unit import tst_cmds

import pytest

class TestMdUtils(object):
    """Test md_utils"""

    def test_misc(self):
        """Misc tests"""

        with pytest.warns(UserWarning):
            with md_utils.errors_as_warnings():
                raise RuntimeError('test error')

        with pytest.raises(KeyboardInterrupt):
            with md_utils.errors_as_warnings():
                raise KeyboardInterrupt('Someone pressed Ctrl-C')

        with pytest.raises(RuntimeError):
            with md_utils.errors_as_warnings(mask_errors=False):
                raise RuntimeError('test error')
        
        with pytest.warns(UserWarning):
            assert md_utils._shell_cmd('abracadabra') == ''
        assert md_utils._shell_cmd('echo hi') == 'hi'

        assert md_utils._mask_secret_info('s3://AWSKEYID:AWSKEYSECRET@bucketname') == 's3://bucketname'

        assert md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'a')
        assert md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'a b')
        assert md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'a b c')
        assert not md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'd')
        assert not md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'd e')
        assert not md_utils.dict_has_keys({}, 'd e')

        assert md_utils.byteify('') == ''
        assert md_utils.byteify(u'') == ''
        assert md_utils.byteify('ABC') == 'ABC'
        assert md_utils.byteify(u'ABC') == 'ABC'
        assert md_utils.byteify([u'A', u'B']) == ['A', 'B']
        assert md_utils.byteify((u'A', u'B')) == ('A', 'B')
        assert md_utils.byteify([u'ABC']) == ['ABC']
        assert md_utils.byteify({u'A': u'ABC'}) == {'A': 'ABC'}
        assert md_utils.byteify({u'A': [u'B', (u'ABC', u'D', {u'E': u'F'})]}) == {'A': ['B', ('ABC', 'D', {'E': 'F'})]}


@contextlib.contextmanager
def step_id_saver():
    """Save and yield the step ID of the command run within the `with` suite"""
    with util.file.tempfname(suffix='.stepid') as step_id_fname, util.misc.tmp_set_env('VIRAL_NGS_METADATA_SAVE_STEP_ID_TO', step_id_fname):
        yield step_id_fname

@pytest.mark.usefixtures('no_detailed_env', 'warnings_as_errors')
class TestMetadataRecording(TestCaseWithTmp):

    def key_matches(self, k, patterns):
        """Test whether a nested-dict key `k` (tuple of keys) matches one of the `patterns`.
        Each pattern is a tuple of strings giving possible values for each position."""
        return any(all(str(k_elt) in p_elt.split() for k_elt, p_elt in zip(k, p) if p_elt) for p in patterns if len(k)>=len(p))

    def test_key_matches(self):
        """Test self.key_matches"""
        assert self.key_matches(('a', 'b', 'c', 'd'), [('a', 'b f g')])
        assert self.key_matches(('a', 'b', 'c', 'd'), [('a', '', 'c')])
        assert not self.key_matches(('a', 'b', 'c', 'd'), [('a', 'a', 'c')])

    def canonicalize_step(self, step_record):
        """From a step record, either remove keys whose values change between execution, or canonicalize the values by replacing
        them with just the value type."""

        return {k: type(v) if self.key_matches(k, (('step', 'run_env run_info run_id step_id version_info'),
                                                   ('step', 'args', '', 'val'),
                                                   ('step', 'args', '', '0 1 2 3', 'val'),
                                                   ('step', 'args', '', 'files', '',
                                                   'abspath ctime device fname inode mtime owner realpath'),
                                                   ('step', 'args', '', '0 1 2 3', 'files', '', 
                                                    'abspath ctime device fname inode mtime owner realpath'),
                                                   ('step', 'metadata_from_cmd_return', 'runtime'))) \
                else v for k, v in util.misc.flatten_dict(step_record, as_dict=(tuple,list)).items() \
                if k[:3] != ('step', 'run_info', 'argv')}  # the command-line here is the py.test invocation, with variable options

    def chk_step(self, step_record, expected_fname):
        """Check the step record `step_record` against expected data from file identified by `expected_fname`"""
        expected_step_record = md_utils.byteify(json.loads(util.file.slurp_file(expected_fname)))
        assert self.canonicalize_step(step_record) == self.canonicalize_step(expected_step_record)

    def test_invalid_metadata_loc(self):
        """Test that metadata-related errors, like logging to non-existent location, result only in warnings"""
        with util.file.tempfname('.txt') as size_fname:
            data1_fname = self.input('data1.txt')

            with pytest.warns(UserWarning):
                with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', '/non/existent/dir'):
                    util.cmd.run_cmd(tst_cmds, 'get_file_size', [data1_fname, size_fname])

    def test_simple_cmd(self):
        """Test basic metadata recording"""
        data1_fname = self.input('data1.txt')
        with util.file.tempfnames(suffixes=('.txt',)) as (size_fname,), step_id_saver() as step_id_fname:
            util.cmd.run_cmd(tst_cmds, 'get_file_size', [data1_fname, size_fname])
            assert util.file.slurp_file(size_fname) == str(os.path.getsize(data1_fname))
            
            step_record = metadata_db.load_step_record(util.file.slurp_file(step_id_fname))
            expected_step1 = self.input('expected.get_file_size.data1.step.json.gz')
            self.chk_step(step_record, expected_step1)
            for arg, fname in ('in_fname', data1_fname), ('size_fname', size_fname):
                assert step_record['step']['args'][arg]['files'][0]['realpath'] == os.path.realpath(fname)

    def test_complex_cmd(self):
        """Test a complex command"""
        with util.file.tempfnames(suffixes=('.info.txt', '.cpy', '.empty')) as \
             (info_fname, cpy_fname, empty_fname), \
             util.file.fifo() as fifo, step_id_saver() as step_id_fname, \
             util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_project', 'testing1'):
            cat = subprocess.Popen(['cat', fifo])
            data1_fname, data1a_fname, ex_pfx, dir_pfx = self.inputs('data1.txt', 'data1a.txt', 'data2', 'data_dir/')

            util.cmd.run_cmd(tst_cmds, 'get_file_info', [data1_fname, os.path.relpath(data1a_fname), info_fname, '--in-fnames-pfx', ex_pfx,
                                                         '--in-fnames-dir', dir_pfx, '--factor', 3, '--make-empty', empty_fname, 
                                                         '--make-empty', fifo,
                                                         '--copy-info-to', cpy_fname, '--metadata', 'purpose', 'testing2'])
            cat.wait()
            assert cat.returncode == 0
            step_record = metadata_db.load_step_record(util.file.slurp_file(step_id_fname))
            expected_step2 = self.input('expected.get_file_info.data1.step.json.gz')
            #util.file.dump_file(self.input('expected.get_file_info.data1.step.json'), json.dumps(step_record, sort_keys=True, indent=4))
            self.chk_step(step_record, expected_step2)

    def test_failing_cmd(self):
        """Test a command that fails"""
        with util.file.tempfname(suffix='.info.txt') as info_fname:
            data1_fname = self.input('data1.txt')
            with step_id_saver() as step_id_fname:
                with pytest.raises(RuntimeError):
                    util.cmd.run_cmd(tst_cmds, 'get_file_info', [data1_fname, info_fname, '--fail'])
                step_record = metadata_db.load_step_record(util.file.slurp_file(step_id_fname))
                expected_step_fail = self.input('expected.get_file_info.fail.step.json.gz')
                self.chk_step(step_record, expected_step_fail)
                #util.file.dump_file(self.input('expected.get_file_info.fail.step.json'), json.dumps(step_record, sort_keys=True, indent=4))

    def test_infile_oneof(self):
        """Test InFile_OneOf"""
        with util.file.tempfnames(suffixes=('.out1', '.out2')) as (out1_fname, out2_fname):
            data1_fname, data2_pfx, data_dir = self.inputs('data1.txt', 'data2', 'data_dir')
            util.cmd.run_cmd(tst_cmds, 'concat_input_files', [data1_fname, out1_fname])
            self.assertEqualContents(data1_fname, out1_fname)

            util.cmd.run_cmd(tst_cmds, 'concat_input_files', [data2_pfx, out1_fname])
            util.file.concat((data2_pfx+'.ex1', data2_pfx+'.ex2'), out2_fname)
            self.assertEqualContents(out1_fname, out2_fname)

            util.cmd.run_cmd(tst_cmds, 'concat_input_files', [data_dir, out1_fname])
            util.file.concat(util.file.files_or_pipes_in_dir(data_dir), out2_fname)
            self.assertEqualContents(out1_fname, out2_fname)
                
# end: class TestMetadataRecording(TestCaseWithTmp)
