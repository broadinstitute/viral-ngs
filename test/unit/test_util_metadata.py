"""Unit tests for util.metadata.

Integration tests are implemented as a pytest plugin in test/conftest.py , which causes the entire test suite to be run
under metadata recording.
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
import concurrent.futures

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

        assert md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'a')
        assert md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'a b')
        assert md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'a b c')
        assert not md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'd')
        assert not md_utils.dict_has_keys(dict(a=1, b=2, c=3), 'd e')
        assert not md_utils.dict_has_keys({}, 'd e')

        assert md_utils.tuple_key_matches(('a', 'b', 'c', 'd'), [('a', 'b f g')])
        assert md_utils.tuple_key_matches(('a', 'b', 'c', 'd'), [('a', '', 'c')])
        assert not md_utils.tuple_key_matches(('a', 'b', 'c', 'd'), [('a', 'a', 'c')])

@contextlib.contextmanager
def step_id_saver():
    """Save and yield the step ID of the command run within the `with` suite"""
    with util.file.tempfname(suffix='.stepid') as step_id_fname, util.misc.tmp_set_env('VIRAL_NGS_METADATA_SAVE_STEP_ID_TO', step_id_fname):
        yield step_id_fname

class TestMetadataRecording(TestCaseWithTmp):

    # def chk_step(self, step_record, expected_fname):
    #     """Check the step record `step_record` against expected data from file identified by `expected_fname`"""
    #     expected_step_record = md_utils.byteify(json.loads(util.file.slurp_file(expected_fname)))
    #     assert self.canonicalize_step(step_record) == self.canonicalize_step(expected_step_record)

    def test_metadata_db(self):
        """Test the basic operations of the metadata database"""

        with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', 's3://AWSKEYID:AWSKEYSECRET@bucketname'):
            assert 'AWSKEYSECRET' in str(metadata_db.metadata_dirs())
            assert 'AWSKEYSECRET' not in str(metadata_db.metadata_dirs_sanitized())

        with util.file.tempfname('.paths.txt') as paths_fn:
            util.file.dump_file(paths_fn, '/tmp/prov')
            with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', '@' + paths_fn):
                assert metadata_db.metadata_dirs() == ['/tmp/prov']

        with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', util.file.get_test_input_path(self)):
            recs = metadata_db.load_all_records()
            assert len(recs) == 3
            with util.file.tmp_dir() as tmp_mdb, util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', tmp_mdb):
                assert len(metadata_db.load_all_records()) == 0
                for r in recs:
                    assert metadata_db.is_valid_step_record(r)
                    metadata_db.store_step_record(r)

                loaded_recs = metadata_db.load_all_records()
                def get_step_id(r): return r['step']['step_id']
                assert sorted(loaded_recs, key=get_step_id) == sorted(recs, key=get_step_id)

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
            
            # step_record = metadata_db.load_step_record(util.file.slurp_file(step_id_fname))
            # expected_step1 = self.input('expected.get_file_size.data1.step.json.gz')
            # self.chk_step(step_record, expected_step1)
            # for arg, fname in ('in_fname', data1_fname), ('size_fname', size_fname):
            #     assert step_record['step']['args'][arg]['files'][0]['realpath'] == os.path.realpath(fname)

    def test_complex_cmd(self):
        """Test a complex command"""
        with util.file.tempfnames(suffixes=('.info.txt', '.cpy', '.empty', '.str1')) as \
             (info_fname, cpy_fname, empty_fname, str1_fname), \
             util.file.fifo(num_pipes=2) as (fifo, fifo2), step_id_saver() as step_id_fname, \
             util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_project', 'testing1'), \
             concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:

            cat = executor.submit(util.file.slurp_file, fifo)
            cat2 = executor.submit(util.file.slurp_file, fifo2)
            data1_fname, data1a_fname, ex_pfx, dir_pfx = self.inputs('data1.txt', 'data1a.txt', 'data2', 'data_dir/')

            res = executor.submit(util.cmd.run_cmd, 'test.unit.tst_cmds', 'get_file_info', 
                                  [data1_fname, os.path.relpath(data1a_fname), info_fname, '--in-fnames-pfx', ex_pfx,
                                   '--in-fnames-dir', dir_pfx, '--factor', 3, '--make-empty', empty_fname, 
                                   '--make-empty', fifo, '--write-data', 'my_data', '--write-dest', fifo2,
                                   '--copy-info-to', cpy_fname, '--metadata', 'purpose', 'testing2'])
            res_result = res.result(timeout=10)
            assert cat.result(timeout=10) == ''
            assert cat2.result(timeout=10) == 'my_data'
#            step_record = metadata_db.load_step_record(util.file.slurp_file(step_id_fname))
#            expected_step2 = self.input('expected.get_file_info.data1.step.json.gz')
            #util.file.dump_file(self.input('expected.get_file_info.data1.step.json'), json.dumps(step_record, sort_keys=True, indent=4))
#            self.chk_step(step_record, expected_step2)

    def test_failing_cmd(self):
        """Test a command that fails"""
        with util.file.tempfname(suffix='.info.txt') as info_fname:
            data1_fname = self.input('data1.txt')
            with step_id_saver() as step_id_fname:
                with pytest.raises(RuntimeError):
                    util.cmd.run_cmd(tst_cmds, 'get_file_info', [data1_fname, info_fname, '--fail'])
#                step_record = metadata_db.load_step_record(util.file.slurp_file(step_id_fname))
#                expected_step_fail = self.input('expected.get_file_info.fail.step.json.gz')
#                self.chk_step(step_record, expected_step_fail)
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
