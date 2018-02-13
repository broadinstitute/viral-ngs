"""Unit tests for util.metadata"""

__author__ = "ilya@broadinstitute.org"
__commands__ = []

import os
import os.path
import argparse
import json
import shutil

import util.cmd
import util.file
import util.misc
from util.metadata import InFile, OutFile
from util._metadata.metadata_db import is_valid_step_record
from test import tst_inp

import reports

import pytest

@pytest.fixture
def metadata_db(tmpdir):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir.mkdir('metadata_db')
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        yield metadata_db_path

def test_metadata_recording(metadata_db):
    """Test basic metadata recording"""
    with util.file.tempfname('.tsv') as metricsFile:
        util.cmd.run_cmd(reports, 'final_assembly_metrics', [tst_inp('TestAssembleSpades/trinity_contigs.fasta'), metricsFile])
        records = metadata_db.listdir()
        assert len(records)==1
        json_str = records[0].read_text(encoding='utf-8')
        step_record = json.loads(json_str)
        assert is_valid_step_record(step_record)
        assert set(step_record['step']['args'].keys()) == set(('assembly_fname', 'metrics_fname'))
        assert step_record['step']['cmd_name'] == 'final_assembly_metrics'
