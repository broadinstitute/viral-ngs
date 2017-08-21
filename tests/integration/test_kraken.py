# Integration tests for kraken

from builtins import super
import argparse
import os.path
from os.path import join
import sys
import tempfile
import pytest
import metagenomics
import util.file
import tools
import tools.kraken
import tools.krona
import tools.picard
from test import TestCaseWithTmp
from test.integration import snake


@pytest.fixture(scope='module')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='module')
def input_bam(request, tmpdir_module, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam = os.path.join(tmpdir_module, 'zaire_ebola.bam')
        fastq_to_sam.execute(fastqs[0], fastqs[1], '', bam)
        return bam

    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='module')
def kraken():
    kraken = tools.kraken.Kraken()
    kraken.install()
    return kraken


@pytest.fixture(scope='module')
def krona():
    krona = tools.krona.Krona()
    krona.install()
    return krona


@pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


@pytest.fixture(scope='module')
def kraken_db(request, tmpdir_module, kraken, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = os.path.join(tmpdir_module, 'kraken_db_{}'.format(db_type))
    os.mkdir(db)
    for d in ['library', 'taxonomy']:
        realpath = join(db_dir, d)
        name = join(db, d)
        os.symlink(realpath, name)

    # Minimizer len corresponds to memory/disk usage of index.
    kraken.build(
        db, options={
            '--minimizer-len': 10,
            '--build': None,
        }
    )
    return db


def test_taxonomy_subset(request, tmpdir_function):
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
    db_dir = os.path.join(data_dir, 'db', 'taxonomy')
    sub_dir = tempfile.mktemp('taxonomy_subset')
    metagenomics.subset_taxonomy(db_dir, sub_dir, whitelistTaxids=[], whitelistTreeTaxids=[186536])

    tax_db = metagenomics.TaxonomyDb(sub_dir, load_nodes=True, load_names=True)
    assert 186536 in tax_db.parents


TAXONOMY_FILES = ('gi_taxid_nucl.dmp',
                  'gi_taxid_prot.dmp',
                  'names.dmp',
                  'nodes.dmp',)


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


def test_kraken_tool(kraken, kraken_db, input_bam):
    outdir = tempfile.mkdtemp('-kraken')
    out = join(outdir, 'zaire_ebola.kraken')
    out_filtered = join(outdir, 'zaire_ebola.filtered-kraken')
    out_report = join(outdir, 'zaire_ebola.kraken-report')
    kraken.classify(input_bam, kraken_db, out)
    kraken.filter(out, kraken_db, out_filtered, 0.05)
    kraken.report(out_filtered, kraken_db, out_report)
    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(out_filtered) > 0


def test_kraken(kraken_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')
    cmd = [input_bam, kraken_db, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(out_reads) > 0


@pytest.mark.skipif(sys.version_info < (3, 2), reason="Python version is too old for snakemake.")
def test_pipes(tmpdir_function, kraken_db, krona_db, input_bam):
    runner = snake.SnakemakeRunner(workdir=tmpdir_function)
    override_config = {
        'kraken_db': kraken_db,
        'krona_db': krona_db,
    }
    runner.set_override_config(override_config)
    runner.setup()
    runner.link_samples([input_bam], destination='per_sample', link_transform=snake.rename_raw_bam)
    runner.create_sample_files(sample_files=['samples_metagenomics'])

    kraken_out = join(
        runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'kraken.report'])
    )

    krona_out = join(
        runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'kraken.krona.html'])
    )

    # runner.run(['all_metagenomics'])
    runner.run([kraken_out, krona_out])
    assert os.path.getsize(kraken_out) > 0
    assert os.path.getsize(krona_out) > 0


def test_kraken_krona(kraken_db, krona_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')

    cmd = [input_bam, kraken_db, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    out_html = util.file.mkstempfname('.krona.html')
    parser = metagenomics.parser_krona(argparse.ArgumentParser())
    args = parser.parse_args([out_reads, krona_db, out_html])
    args.func_main(args)
