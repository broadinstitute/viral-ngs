import os.path
from os.path import join

import util.file
import metagenomics


def test_taxonomy_subset_zaire(request, tmpdir_factory):
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
    db_dir = join(data_dir, 'db', 'taxonomy')
    sub_dir = str(tmpdir_factory.mktemp('taxonomy_subset'))
    # Zaire species
    metagenomics.subset_taxonomy(db_dir, sub_dir, whitelistTaxids=[], whitelistTreeTaxids=[186538])

    tax_db = metagenomics.TaxonomyDb(sub_dir, load_nodes=True, load_names=True)
    assert 186538 in tax_db.parents  # Zaire species
    assert 186540 not in tax_db.parents  # Sudan species
    assert 2 not in tax_db.parents  # Bacteria
