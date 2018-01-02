# Integration tests for diamond

from builtins import super
import argparse
import fnmatch
import gzip
import os
from os.path import join
import sys
import shutil
import tempfile
import pytest
from Bio import SeqIO

import metagenomics
import tools
import tools.diamond
import tools.picard
import util.file

def find_files(root_dir, filt):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, filt):
            yield join(root, filename)


@pytest.fixture(scope='module')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='module')
def sam_to_fastq():
    return tools.picard.SamToFastqTool()


@pytest.fixture(scope='module')
def diamond():
    diamond = tools.diamond.Diamond()
    diamond.install()
    return diamond


@pytest.fixture(scope='module')
def krona():
    krona = tools.krona.Krona()
    krona.install()
    return krona


@pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


def input_fastq_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
    return [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]


def input_bam_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='module')
def input_bam(request, tmpdir_module, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam_name = 'zaire_ebola.bam'
        bam = os.path.join(tmpdir_module, bam_name)
        fastq_to_sam.execute(fastqs[0], fastqs[1], '', bam)
        return bam

    return input_bam_paths()


@pytest.fixture(scope='module')
def input_fastqs(request, tmpdir_module, sam_to_fastq, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        return fastqs

    bam = join(data_dir, 'test-reads.bam')
    basename = os.path.basename(bam)
    fastq1 = join(tmpdir_module, '{}.1.fastq'.format(basename))
    fastq2 = join(tmpdir_module, '{}.2.fastq'.format(basename))
    sam_to_fastq.execute(bam, fastq1, fastq2)
    return fastq1, fastq2


@pytest.fixture(scope='module')
def taxonomy_db(request, tmpdir_module, db_type):
    taxonomy = os.path.join(tmpdir_module, db_type, 'taxonomy')
    shutil.copytree(join(util.file.get_test_input_path(), db_type, 'db', 'taxonomy'),
                    taxonomy)
    prot = os.path.join(taxonomy, 'accession2taxid', 'prot.accession2taxid')
    prot_gz = prot + '.gz'

    with open(prot, 'rb') as f_in:
        with gzip.open(prot_gz, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return taxonomy


@pytest.fixture(scope='module')
def diamond_db(request, tmpdir_module, diamond, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = os.path.join(tmpdir_module, db_type + '.dmnd')
    translated = os.path.join(tmpdir_module, db_type + '.faa')

    lib_dir = join(db_dir, 'library')
    util.file.cat(translated, find_files(db_dir, '*.faa'))

    diamond.build(db, [translated])
    return os.path.splitext(db)[0]


TAXONOMY_FILES = ('gi_taxid_nucl.dmp',
                  'gi_taxid_prot.dmp',
                  'names.dmp',
                  'nodes.dmp',
                  'merged.dmp')


@pytest.fixture(scope='module')
def krona_db(request, tmpdir_module, krona, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = os.path.join(data_dir, 'db')

    db = os.path.join(tmpdir_module, 'krona_db_{}'.format(db_type))
    os.mkdir(db)
    for d in TAXONOMY_FILES:
        src = join(db_dir, 'taxonomy', d)
        dest = join(db, d)
        os.symlink(src, dest)
    krona.create_db(db)
    return db


@pytest.mark.skipif(tools.is_osx(), reason="not currently tested under OSX")
def test_diamond(diamond_db, taxonomy_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.lca.tsv')
    cmd = [input_bam, diamond_db, taxonomy_db, out_report, '--outReads', out_reads]
    parser = metagenomics.parser_diamond(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(out_reads) > 0
