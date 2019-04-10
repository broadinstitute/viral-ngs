import gzip
import os
from os.path import join
import shutil
import tools
import tools.kaiju
import tools.picard
import tools.krona
import util.file

# @pytest.fixture(scope='module')
def krona():
    krona = tools.krona.Krona()
    krona.install()
    return krona


# @pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


# @pytest.fixture(scope='module')
def krona_db(request, tmpdir_module, krona, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = join(tmpdir_module, 'krona_db_{}'.format(db_type))
    os.mkdir(db)
    for d in ('names.dmp', 'nodes.dmp'):
        src = join(db_dir, 'taxonomy', d)
        dest = join(db, d)
        if os.path.isfile(src):
            os.symlink(src, dest)
    krona.create_db(db)
    return db


def taxonomy_db(request, tmpdir_module, db_type):
    taxonomy = join(tmpdir_module, db_type, 'taxonomy')
    shutil.copytree(join(util.file.get_test_input_path(), db_type, 'db', 'taxonomy'),
                    taxonomy)
    prot = join(taxonomy, 'accession2taxid', 'prot.accession2taxid')
    prot_gz = prot + '.gz'

    with open(prot, 'rb') as f_in:
        with gzip.open(prot_gz, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return taxonomy
