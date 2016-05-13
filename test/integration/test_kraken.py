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
import tools.kraken
import tools.krona
import tools.picard
from test import TestCaseWithTmp
from test.integration.snake import SnakemakeRunner


@pytest.fixture(autouse=True)
def set_tempdir(request):
    util.file.set_tmp_dir(None)
    request.addfinalizer(util.file.destroy_tmp_dir)


@pytest.fixture(scope='session')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='session')
def sam_to_fastq():
    return tools.picard.SamToFastqTool()


def input_fastq_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
    return [os.path.join(data_dir, f)
            for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]


def input_bam_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='session')
def input_bam(request, tmpdir_factory, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [os.path.join(data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam_name = 'zaire_ebola.bam'
        bam = str(tmpdir_factory.getbasetemp().join(bam_name))
        fastq_to_sam.execute(fastqs[0], fastqs[1], '', bam)
        return bam

    return input_bam_paths()


@pytest.fixture(scope='session')
def input_fastqs(request, tmpdir_factory, sam_to_fastq, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [join(data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        return fastqs

    bam = join(data_dir, 'test-reads.bam')
    basename = os.path.basename(bam)
    fastq1 = join(str(tmpdir_factory.getbasetemp()),
                  '{}.1.fastq'.format(basename))
    fastq2 = join(str(tmpdir_factory.getbasetemp()),
                  '{}.2.fastq'.format(basename))
    sam_to_fastq.execute(bam, fastq1, fastq2)
    return fastq1, fastq2


@pytest.fixture(scope='session')
def kraken():
    kraken = tools.kraken.Kraken()
    kraken.install()
    return kraken


@pytest.fixture(scope='session')
def krona():
    krona = tools.krona.Krona()
    krona.install()
    return krona


@pytest.fixture(scope='session',
                params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


@pytest.fixture(scope='session')
def kraken_db(request, tmpdir_factory, kraken, db_type):

    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = str(tmpdir_factory.mktemp('kraken_db_{}'.format(db_type)))
    for d in ['library', 'taxonomy']:
        realpath = join(db_dir, d)
        name = join(db, d)
        os.symlink(realpath, name)

    # Minimizer len corresponds to memory/disk usage of index.
    assert kraken.build(db, options={
        '--minimizer-len': 10,
        '--build': None,
    }).returncode == 0
    return db


TAXONOMY_FILES = (
    'gi_taxid_nucl.dmp',
    'gi_taxid_prot.dmp',
    'names.dmp',
    'nodes.dmp',
)


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


def test_kraken_tool(tmpdir, kraken, kraken_db, input_fastqs):
    outdir = tempfile.mkdtemp('-kraken')
    out = join(outdir, 'zaire_ebola.kraken')
    out_filtered = join(outdir, 'zaire_ebola.filtered-kraken')
    out_report = join(outdir, 'zaire_ebola.kraken-report')
    assert kraken.classify(kraken_db, input_fastqs, out).returncode == 0
    result = kraken.execute(
        'kraken-filter', kraken_db, out_filtered, [out],
        options={'--threshold': 0.05})
    assert result.returncode == 0
    result = kraken.execute(
        'kraken-report', kraken_db, out_report, [out_filtered])
    assert result.returncode == 0

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


@pytest.mark.skipif(sys.version_info < (3,2),
                    reason="Python version is too old for snakemake.")
def test_pipes(tmpdir, kraken_db, krona_db, input_bam):
    runner = SnakemakeRunner(workdir=str(tmpdir))
    override_config = {
        'kraken_db': kraken_db,
        'krona_db': krona_db,
    }
    runner.set_override_config(override_config)
    runner.setup()
    runner.link_samples([input_bam], destination='source')
    runner.create_sample_files(sample_files=['samples_metagenomics'])

    kraken_out = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                         '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'kraken.report']))

    krona_out = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                         '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'kraken.krona.html']))

    # runner.run(['all_metagenomics'])
    runner.run([kraken_out])
    assert os.path.getsize(kraken_out) > 0
    # assert os.path.getsize(krona_out) > 0


@pytest.mark.skipif(True, reason="krona create db takes too much disk io")
def test_kraken_krona(tmpdir, kraken_db, krona_db, input_bam):
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
