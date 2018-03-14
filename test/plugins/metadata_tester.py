"""pytest plugin for testing metadata recording"""

import os
import os.path
import functools
import operator
import uuid
import collections

import util.misc
import util.file
import util._metadata.metadata_db as metadata_db
import util._metadata.md_utils as md_utils

import pytest

def canonicalize_step_record(step_record, canonicalize_keys=()):
    """Return a canonicalized flat dict of key-value pairs representing this record, for regression testing purposes.
    Canonicalization uniformizes or drops fields that may change between runs.
    """
    def canonicalize_value(v):
        if v is None: return v
        return type(v)()
    pfx = (step_record['step']['cmd_name'],)
    return {pfx+k: canonicalize_value(v) \
            if md_utils.tuple_key_matches(k, canonicalize_keys) \
            else v for k, v in util.misc.flatten_dict(step_record, as_dict=(tuple,list)).items() \
            if k[:3] != ('step', 'run_info', 'argv')}

def canonicalize_step_records(step_records, canonicalize_keys=()):
    """Canonicalize a group of step records"""
    return sorted(map(str, functools.reduce(operator.concat, 
                                            [list(canonicalize_step_record(r, canonicalize_keys).items())
                                             for r in step_records], [])))

@pytest.fixture(scope='session', autouse='true')
def tmp_metadata_db(tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = os.environ.get('VIRAL_NGS_TEST_METADATA_PATH', tmpdir_factory.mktemp('metadata_db'))
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        yield metadata_db_path        

#@pytest.fixture(autouse='true')
def per_test_metadata_db(request, tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir_factory.mktemp('metadata_db')

    # The following tests, for varying reasons, result in some captured metadata varying between test runs.
    # For now, just don't attempt to check that the metadata does not change; eventually some/most of these exceptions
    # can be removed.
    nondet_tests = (
        'test/unit/test_assembly.py::TestDeambigAndTrimFasta::test_deambig_fasta',
    )

    canonicalize_keys = (
        ('step', 'run_env run_info run_id step_id version_info'),
        ('step', 'args', '', 'val'),
        ('step', 'args', '', ' '.join(map(str, range(10))), 'val'),
        ('step', 'args', '', 'files', '',
         'abspath ctime device fname inode mtime owner realpath'),
        ('step', 'args', '', ' '.join(map(str, range(10))), 'files', '',
         'abspath ctime device fname inode mtime owner realpath'),
        ('step', 'metadata_from_cmd_return', 'runtime'),
        ('step', 'enclosing_steps'),
        ('step', 'args', 'tmp_dir'),
        ('step', 'args', 'tmp_dirKeep'),
        ('step', 'args', 'novo_params'),
        ('step', 'args', 'refDbs'),
        ('step', 'args', 'outputDirectory', 'files', '2'),
        ('step', 'args', 'inVcf', 'files', '0'),
        ('step', 'args', 'blastDbs', '0', 'files', '2'),
    )

    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        print('metadata_db_path:', metadata_db_path)
        yield metadata_db_path
        recs_canon = canonicalize_step_records(metadata_db.load_all_records(), canonicalize_keys)

        cmd_rec_fname = os.path.join(util.file.get_test_input_path(), 'cmd', util.file.string_to_file_name(request.node.nodeid))
        if os.path.isfile(cmd_rec_fname):
            expected_lines = util.file.slurp_file(cmd_rec_fname).strip().split('\n')
            if recs_canon != expected_lines:
                pytest.fail('lines do not match: \n{}'.format('\n'.join(sorted(set(recs_canon) ^ set(expected_lines)))))
        elif 'VIRAL_NGS_TEST_GATHER_CMDS' in os.environ and recs_canon and \
             request.node.nodeid not in nondet_tests:
            util.file.dump_file(cmd_rec_fname, '\n'.join(recs_canon))

@pytest.fixture(scope='session', autouse='true')
def no_detailed_env():
    """Disable time-consuming gathering of detailed env"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_DETAILED_ENV', None):
        yield

@pytest.fixture(autouse='true')
def set_nodeid_metadata(request):
    """Add pytest nodeid and a unique ID of this test invocation to the metadata"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_nodeid', request.node.nodeid), \
         util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_testid', uuid.uuid4()):
        yield

def gather_metadata_regtests():
    """Gather data for metadata regtesting."""
    
    recs = metadata_db.load_all_records()
    nodeid2testid2recs = collections.defaultdict(functools.partial(collections.defaultdict, list))
    for rec in recs:
        pfx = (rec['step']['cmd_name'],)
        rec_canon = {pfx+k: v for k, v in util.misc.flatten_dict(rec, as_dict=(tuple, list)).items()}
        mdata = rec['step']['metadata_from_cmd_line']
        nodeid2testid2recs[mdata['pytest_nodeid']][mdata['pytest_testid']] = rec_canon

    for nodeid, testid2recs in nodeid2testid2recs.items():
        common_keys = set.intersection(*[set(k2v.keys()) for k2v in testid2recs.values()])
        regtest_k2v = {}
        for k in common_keys:
            vals = { k2v[k]for k2v in testid2recs.values() }
            if len(vals) == 1:
                regtest_k2v[k] = list(vals)[0]
            else:
                def canonicalize_value(v):
                    if v is None: return v
                    return type(v)()

                val_types = set(map(canonicalize_value, vals))
                assert len(val_types) == 1
                regtest_k2v[k] 
        
