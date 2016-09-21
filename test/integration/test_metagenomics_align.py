# Integration tests for metagenomics direct alignment

from builtins import super
import argparse
import fnmatch
import os.path
from os.path import join
import sys
import tempfile

import pytest
from Bio import SeqIO

import metagenomics
import util.file
import tools
import tools.bwa
import tools.krona
import tools.picard
from test.integration import snake


def find_files(root_dir, filt):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, filt):
            yield join(root, filename)


@pytest.fixture(autouse=True, scope='session')
def set_tempdir(request):
    util.file.set_tmp_dir(None)
    request.addfinalizer(util.file.destroy_tmp_dir)


@pytest.fixture(scope='session')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='session')
def taxonomy_db(request, tmpdir_factory, db_type):
    return join(util.file.get_test_input_path(), db_type, 'db', 'taxonomy')


@pytest.fixture(scope='session')
def input_bam(request, tmpdir_factory, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam_name = 'zaire_ebola.bam'
        bam = str(tmpdir_factory.getbasetemp().join(bam_name))
        fastq_to_sam.execute(fastqs[0], fastqs[1], '', bam)
        return bam

    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='session')
def bwa():
    bwa = tools.bwa.Bwa()
    bwa.install()
    return bwa


# @pytest.fixture(scope='session', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
@pytest.fixture(scope='session', params=['TestMetagenomicsSimple'])
def db_type(request):
    return request.param


FNA_TAXIDS = {
    'NC_014373.fna': '565995',    # Bundibugyo_ebolavirus_uid51245
    'NC_004161.fna': '186539',    # Reston_ebolavirus_uid15006
    'NC_006432.fna': '186540',    # Sudan_ebolavirus_uid15012
    'NC_014372.fna': '186541',    # Tai_Forest_ebolavirus_uid51257
    'NC_002549.fna': '186538',    # Zaire_ebolavirus_uid14703
}


@pytest.fixture(scope='session')
def bwa_db(request, tmpdir_factory, bwa, db_type):

    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    index_fa = str(tmpdir_factory.getbasetemp().join(db_type + '.bwa_index.fa'))
    db = str(tmpdir_factory.getbasetemp().join(db_type + '.bwa'))

    with open(index_fa, "w") as f_out:
        for fname in find_files(join(db_dir, 'library'), '*.fna'):
            with open(fname) as f:
                for seq_record in SeqIO.parse(f, 'fasta'):
                    seq_id = seq_record.id
                    try:
                        tax_id = FNA_TAXIDS[os.path.basename(fname)]
                    except KeyError:
                        continue
                    seq_record.id = 'taxid|{}|{}'.format(tax_id, seq_id)
                    SeqIO.write(seq_record, f_out, 'fasta')

    bwa.index(index_fa, output=db)
    return db


def test_meta_bwa(bwa_db, taxonomy_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_bam = util.file.mkstempfname('.output.bam')
    cmd = [input_bam, bwa_db, taxonomy_db, out_report, '--outBam', out_bam]
    parser = metagenomics.parser_align_rna_metagenomics(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(out_bam) > 0


@pytest.mark.skipif(sys.version_info < (3, 2), reason="Python version is too old for snakemake.")
def test_pipes(tmpdir, bwa_db, taxonomy_db, input_bam):
    runner = snake.SnakemakeRunner(workdir=str(tmpdir))
    override_config = {
        'align_rna_db': bwa_db,
        'taxonomy_db': taxonomy_db,
    }
    runner.set_override_config(override_config)
    runner.setup()
    runner.link_samples([input_bam], destination='per_sample', link_transform=snake.rename_raw_bam)
    runner.create_sample_files(sample_files=['samples_metagenomics'])

    report_out = join(
        runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'rna_bwa.report'])
    )

    bam_out = join(
        runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
        '.'.join([os.path.splitext(os.path.basename(input_bam))[0], 'rna_bwa.bam'])
    )

    runner.run([report_out])
    assert os.path.getsize(report_out) > 0
    assert os.path.getsize(bam_out) > 0
