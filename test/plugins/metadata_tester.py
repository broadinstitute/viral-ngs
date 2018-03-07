"""pytest plugin for testing metadata recording"""

import os
import os.path

import util.misc
import util.file
import util._metadata.metadata_db as metadata_db

import pytest

#@pytest.fixture(scope='session', autouse='true')
def tmp_metadata_db(tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = os.environ.get('VIRAL_NGS_TEST_METADATA_PATH', tmpdir_factory.mktemp('metadata_db'))
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        print('metadata_db_path:', metadata_db_path)
        yield metadata_db_path

@pytest.fixture(autouse='true')
def per_test_metadata_db(request, tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir_factory.mktemp('metadata_db')
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        print('metadata_db_path:', metadata_db_path)
        yield metadata_db_path
        recs_canon = metadata_db.canonicalize_step_records(metadata_db.load_all_records())

        cmd_rec_fname = os.path.join(util.file.get_test_input_path(), 'cmd', util.file.string_to_file_name(request.node.nodeid))
        if os.path.isfile(cmd_rec_fname):
            expected_lines = util.file.slurp_file(cmd_rec_fname).strip().split('\n')
            if recs_canon != expected_lines:
                pytest.fail('lines do not match: \n{}'.format('\n'.join(sorted(set(recs_canon) ^ set(expected_lines)))))
        elif 'VIRAL_NGS_TEST_GATHER_CMDS' in os.environ and recs_canon:
            util.file.dump_file(cmd_rec_fname, '\n'.join(recs_canon))

@pytest.fixture(scope='session', autouse='true')
def no_detailed_env():
    """Disable time-consuming gathering of detailed env"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_DETAILED_ENV', None):
        yield

@pytest.fixture(autouse='true')
def set_nodeid_metadata(request):
    """Add pytest nodeid to the metadata"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_VALUE_pytest_nodeid', request.node.nodeid):
        yield
