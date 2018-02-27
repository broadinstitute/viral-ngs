import os
import warnings

import util.misc

import pytest

@pytest.fixture(scope='session', autouse='true')
def tmp_metadata_db(tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = os.environ.get('VIRAL_NGS_TEST_METADATA_PATH', tmpdir_factory.mktemp('metadata_db'))
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        print('metadata_db_path:', metadata_db_path)
        yield metadata_db_path

@pytest.fixture(scope='session', autouse='true')
def no_detailed_env():
    """Disable time-consuming gathering of detailed env"""
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_DETAILED_ENV', None):
        yield

@pytest.fixture(scope='session', autouse='true')
def warnings_as_errors():
    """Turn warnings into errors, so we see them during testing"""
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)
        yield
