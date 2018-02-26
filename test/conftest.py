
import shutil

import util.misc
import util._metadata.metadata_db

import pytest

@pytest.fixture(scope='session')
def tmp_metadata_db(tmpdir_factory):
    """Sets up the metadata database in a temp dir"""
    metadata_db_path = tmpdir_factory.mktemp('metadata_db')
    with util.misc.tmp_set_env('VIRAL_NGS_METADATA_PATH', metadata_db_path):
        yield metadata_db_path
