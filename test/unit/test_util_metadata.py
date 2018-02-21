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
from util._metadata import metadata_db, file_arg
from test import tst_inp

from test.unit import tst_cmds

import pytest

@pytest.fixture
def tmp_metadata_db(tmpdir):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir.mkdir('metadata_db')
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        yield metadata_db_path

def test_metadata_recording(tmp_metadata_db):
    """Test basic metadata recording"""
    with util.file.tempfname('.txt') as size_fname:
        data1_fname = tst_inp('TestUtilCmd/data1.txt')
        util.cmd.run_cmd(tst_cmds, 'get_file_size', [data1_fname, size_fname])
        assert util.file.slurp_file(size_fname) == str(os.path.getsize(data1_fname))
        records = metadata_db.load_all_records()
        assert len(records)==1
        step_record = records[0]

        def canonicalize_step(s):
            return {k: type(v) if len(k)>1 and k[1] in ('run_env', 'run_info', 'run_id', 'step_id', 'version_info') or \
                    len(k)>3 and k[1]=='args' and (k[3]=='val' or k[3]=='files' and k[5] not in ('hash', 'size')) \
                    else v for k, v in util.misc.flatten_dict(s, as_dict=(tuple,list)).items() \
                    if k[:3] != ('step', 'run_info', 'argv')}

        def chk_step(s, fn):
            assert canonicalize_step(s) == \
                canonicalize_step(json.loads(util.file.slurp_file(tst_inp('TestUtilCmd/expected.{}.step.json.gz'.format(fn)))))

        chk_step(step_record, 'get_file_size.data1')


def test_metadata_recording2(tmp_metadata_db):
    """Test basic metadata recording"""
    with util.file.tempfname('.tsv') as metricsFile:
        import reports
        util.cmd.run_cmd(reports, 'final_assembly_metrics', [tst_inp('TestAssembleSpades/trinity_contigs.fasta'), metricsFile])
        records = metadata_db.load_all_records()
        assert len(records)==1
        step_record = records[0]
        assert metadata_db.is_valid_step_record(step_record)
        assert set(step_record['step']['args'].keys()) == set(('assembly_fname', 'metrics_fname'))
        assert step_record['step']['cmd_module'] == 'reports'
        assert step_record['step']['cmd_name'] == 'final_assembly_metrics'
