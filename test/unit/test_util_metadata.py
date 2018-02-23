"""Unit tests for util.metadata"""

__author__ = "ilya@broadinstitute.org"
__commands__ = []

import os
import os.path
import collections
import argparse
import json
import shutil

import util.cmd
import util.file
import util.misc
from util.metadata import InFile, OutFile
from util._metadata import metadata_db, file_arg, md_utils
from test import tst_inp, TestCaseWithTmp

from test.unit import tst_cmds

import pytest

@pytest.fixture
def tmp_metadata_db(tmpdir):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir.mkdir('metadata_db')
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        yield metadata_db_path

@pytest.fixture
def no_detailed_env():
    """Disable time-consuming gathering of detailed env"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_DETAILED_ENV', None):
        yield

@pytest.mark.usefixtures('tmp_metadata_db', 'no_detailed_env')
class TestMetadataRecording(TestCaseWithTmp):

    def canonicalize_step(self, step_record):
        """From a step record, either remove keys whose values change between execution, or canonicalize the values by replacing
        them with just the value type."""
        return {k: type(v) if len(k)>1 and k[1] in ('run_env', 'run_info', 'run_id', 'step_id', 'version_info') or \
                len(k)>3 and k[1]=='args' and (k[3]=='val' or k[3]=='files' and k[5] not in ('hash', 'size')) \
                else v for k, v in util.misc.flatten_dict(step_record, as_dict=(tuple,list)).items() \
                if k[:3] != ('step', 'run_info', 'argv')}

    def chk_step(self, step_record, expected_fname):
        """Check the step record `step_record` against expected data from file identified by `expected_fname`"""
        expected_step_record = md_utils.byteify(json.loads(util.file.slurp_file(expected_fname)))
        assert self.canonicalize_step(step_record) == self.canonicalize_step(expected_step_record)

    def test_simple_cmd(self):
        """Test basic metadata recording"""
        with util.file.tempfname('.txt') as size_fname:
            data1_fname = self.input('data1.txt')
            util.cmd.run_cmd(tst_cmds, 'get_file_size', [data1_fname, size_fname])
            assert util.file.slurp_file(size_fname) == str(os.path.getsize(data1_fname))
            records = metadata_db.load_all_records()
            assert len(records)==1
            step_record = records[0]

            expected_step1 = self.input('expected.get_file_size.data1.step.json.gz')

            self.chk_step(step_record, expected_step1)
            for arg, fname in ('in_fname', data1_fname), ('size_fname', size_fname):
                assert step_record['step']['args'][arg]['files'][0]['realpath'] == os.path.realpath(fname)


        with util.file.tempfnames(suffixes=('.info.txt', '.cpy', '.empty')) as (info_fname, cpy_fname, empty_fname):
            data1_fname, data1a_fname, ex_pfx, dir_pfx = self.inputs('data1.txt', 'data1a.txt', 'data2', 'data_dir/')
            util.cmd.run_cmd(tst_cmds, 'get_file_info', [data1_fname, data1a_fname, info_fname, '--in-fnames-pfx', ex_pfx,
                                                         '--in-fnames-dir', dir_pfx, '--factor', '3', '--make-empty', empty_fname,
                                                         '--copy-info-to', cpy_fname])
            records = metadata_db.load_all_records()
            assert len(records)==2
            step_record = [r for r in records if r['step']['cmd_name']=='get_file_info'][0]

# end: class TestMetadataRecording(TestCaseWithTmp)
