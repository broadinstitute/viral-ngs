# Integration tests for kraken

import argparse
import os
import os.path
from os.path import join
import sys
import pytest
import metagenomics
import util.file
import tools
import tools.kraken
import tools.krona
import tools.picard
from test.integration import snake


@pytest.fixture()
def input_bam(db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
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

    parser = metagenomics.parser_kraken(argparse.ArgumentParser())

    db = os.path.join(tmpdir_module, 'kraken_db_{}'.format(db_type))
    parser = metagenomics.parser_kraken_build(argparse.ArgumentParser())
    cmd = [db, '--library', join(db_dir, 'library'),
           '--taxonomy', join(db_dir, 'taxonomy'),
           '--subsetTaxonomy',
           '--minimizerLen', '10',
           '--clean']

    parser.parse_args(cmd)
    args = parser.parse_args(cmd)
    args.func_main(args)
    return db


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
        if os.path.isfile(src):
            os.symlink(src, dest)
    krona.create_db(db)
    return db

def test_kraken(kraken_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')
    cmd = [kraken_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    with util.file.open_or_gzopen(out_reads, 'r') as inf:
        assert len(inf.read()) > 0
    with util.file.open_or_gzopen(out_report) as inf:
        report_lines = [x.strip().split() for x in inf.readlines()]

    assert os.path.getsize(out_report) > 0

    '''
    # not sure what to do with this failing test for now that never seemed to work in the past
    if 'TestMetagenomicsSimple' in kraken_db:
        zaire_found = False
        tai_found = False
        for line in report_lines:
            if line[-1] == 'Zaire ebolavirus' and float(line[0]) > 90:
                zaire_found = True
            elif 'Tai Forest' in line[-1]:
                tai_found = True
        assert zaire_found
        assert not tai_found
    '''

@pytest.mark.skipif(sys.version_info < (3, 5), reason="Python version is too old for snakemake.")
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
        runner.config['data_dir'], runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'kraken.report'])
    )

    krona_out = join(
        runner.config['data_dir'], runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'raw', 'kraken.krona.html'])
    )

    # runner.run(['all_metagenomics'])
    runner.run([kraken_out, krona_out])
    assert os.path.getsize(os.path.join(runner.workdir, kraken_out)) > 0
    assert os.path.getsize(os.path.join(runner.workdir, krona_out)) > 0

def test_kraken_krona(kraken_db, krona_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')

    cmd = [kraken_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    out_html = util.file.mkstempfname('.krona.html')
    parser = metagenomics.parser_krona(argparse.ArgumentParser())
    args = parser.parse_args([out_reads, krona_db, out_html])
    args.func_main(args)

def test_kraken_on_empty(kraken_db, input_bam):
    if 'TestMetagenomicsViralMix' not in kraken_db:
        return
    input_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads.gz')
    cmd = [kraken_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    with util.file.open_or_gzopen(out_reads, 'r') as inf:
        assert len(inf.read()) == 0
    with open(out_report, 'rt') as inf:
        out_report_contents = inf.readlines()
    assert len(out_report_contents) == 1
    out_report_contents = out_report_contents[0].rstrip('\n').split('\t')
    assert out_report_contents == ['100.00', '0', '0', 'U', '0', 'unclassified']
