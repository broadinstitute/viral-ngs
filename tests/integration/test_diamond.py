# Integration tests for diamond

from builtins import super
import argparse
import fnmatch
import os
from os.path import join
import sys
import tempfile
import pytest
from Bio import SeqIO
import metagenomics
from test.integration import snake
import tools
import tools.diamond
import tools.picard
import util.file


def find_files(root_dir, filt):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, filt):
            yield join(root, filename)


@pytest.fixture(scope='session', autouse=True)
def set_tempdir(request):
    util.file.set_tmp_dir(None)
    request.addfinalizer(util.file.destroy_tmp_dir)


@pytest.fixture(scope='session')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='session')
def sam_to_fastq():
    return tools.picard.SamToFastqTool()


@pytest.fixture(scope='session')
def diamond():
    diamond = tools.diamond.Diamond()
    diamond.install()
    return diamond


@pytest.fixture(scope='session')
def krona():
    krona = tools.krona.Krona()
    krona.install()
    return krona


@pytest.fixture(scope='session', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


def input_fastq_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
    return [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]


def input_bam_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='session')
def input_bam(request, tmpdir_factory, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam_name = 'zaire_ebola.bam'
        bam = str(tmpdir_factory.getbasetemp().join(bam_name))
        fastq_to_sam.execute(fastqs[0], fastqs[1], '', bam)
        return bam

    return input_bam_paths()


@pytest.fixture(scope='session')
def input_fastqs(request, tmpdir_factory, sam_to_fastq, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        return fastqs

    bam = join(data_dir, 'test-reads.bam')
    basename = os.path.basename(bam)
    fastq1 = join(str(tmpdir_factory.getbasetemp()), '{}.1.fastq'.format(basename))
    fastq2 = join(str(tmpdir_factory.getbasetemp()), '{}.2.fastq'.format(basename))
    sam_to_fastq.execute(bam, fastq1, fastq2)
    return fastq1, fastq2


@pytest.fixture(scope='session')
def taxonomy_db(request, tmpdir_factory, db_type):
    return join(util.file.get_test_input_path(), db_type, 'db', 'taxonomy')


@pytest.fixture(scope='session')
def diamond_db(request, tmpdir_factory, diamond, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = str(tmpdir_factory.getbasetemp().join(db_type + '.dmnd'))
    translated = str(tmpdir_factory.getbasetemp().join(db_type + '.fa'))

    lib_dir = join(db_dir, 'library')
    with open(translated, "w") as f_out:
        for fname in find_files(db_dir, '*.ffn'):
            with open(fname) as f:
                for seq_record in SeqIO.parse(f, 'fasta'):
                    seq_record.seq = seq_record.seq.translate()
                    SeqIO.write(seq_record, f_out, 'fasta')
    diamond.build(db, [translated])
    return db


TAXONOMY_FILES = ('gi_taxid_nucl.dmp',
                  'gi_taxid_prot.dmp',
                  'names.dmp',
                  'nodes.dmp',)


@pytest.fixture(scope='session')
def krona_db(request, tmpdir_factory, krona, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = os.path.join(data_dir, 'db')

    db = str(tmpdir_factory.mktemp('krona_db_{}'.format(db_type)))
    for d in TAXONOMY_FILES:
        src = join(db_dir, 'taxonomy', d)
        dest = join(db, d)
        os.symlink(src, dest)
    krona.create_db(db)
    return db


@pytest.mark.skipif(sys.version_info < (3, 2), reason="Python version is too old for snakemake.")
def test_pipes(tmpdir, diamond_db, taxonomy_db, krona_db, input_bam):
    runner = snake.SnakemakeRunner(workdir=str(tmpdir))
    override_config = {
        'diamond_db': diamond_db,
        'taxonomy_db': taxonomy_db,
        'krona_db': krona_db,
    }
    runner.set_override_config(override_config)
    runner.setup()
    runner.link_samples([input_bam], destination='per_sample', link_transform=snake.rename_raw_bam)
    runner.create_sample_files(sample_files=['samples_metagenomics'])

    krona_out = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                         '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'diamond.krona.html']))

    diamond_out = join(
        runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'diamond.report'])
    )
    runner.run([krona_out])
    assert os.path.getsize(diamond_out) > 0
    assert os.path.getsize(krona_out) > 0


def test_diamond(diamond_db, taxonomy_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_lca = util.file.mkstempfname('.lca.tsv')
    cmd = [input_bam, diamond_db, taxonomy_db, out_report, '--outLca', out_lca]
    parser = metagenomics.parser_diamond(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(out_lca) > 0
