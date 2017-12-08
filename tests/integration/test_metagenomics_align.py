# Integration tests for metagenomics direct alignment

from builtins import super
import argparse
import fnmatch
from os import listdir
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

def find_files(root_dir, filt):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, filt):
            yield join(root, filename)


@pytest.fixture(scope='module')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='module')
def taxonomy_db(request, tmpdir_factory, db_type):
    return join(util.file.get_test_input_path(), db_type, 'db', 'taxonomy')


@pytest.fixture(scope='module')
def input_bam(request, tmpdir_module, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [os.path.join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam_name = 'zaire_ebola.bam'
        bam = os.path.join(tmpdir_module, bam_name)
        fastq_to_sam.execute(fastqs[0], fastqs[1], '', bam)
        return bam

    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='module')
def bwa():
    bwa = tools.bwa.Bwa()
    bwa.install()
    return bwa


# @pytest.fixture(scope='session', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
@pytest.fixture(scope='module', params=['TestMetagenomicsSimple'])
def db_type(request):
    return request.param


FNA_TAXIDS = {
    'GCF_000889155.1_ViralProj51245_genomic.fna': '565995',    # Bundibugyo_ebolavirus
    'GCF_000854085.1_ViralProj15006_genomic.fna': '186539',    # Reston_ebolavirus
    'GCF_000855585.1_ViralProj15012_genomic.fna': '186540',    # Sudan_ebolavirus
    'GCF_000888475.1_ViralProj51257_genomic.fna': '186541',    # Tai_Forest_ebolavirus
    'GCF_000848505.1_ViralProj14703_genomic.fna': '186538',    # Zaire_ebolavirus
}


@pytest.fixture(scope='module')
def bwa_db(request, tmpdir_module, bwa, db_type):

    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    index_fa = os.path.join(tmpdir_module, db_type + '.bwa_index.fa')
    db = os.path.join(tmpdir_module, db_type + '')

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
    dupe_report = util.file.mkstempfname('.dupes.report')
    out_bam = util.file.mkstempfname('.output.bam')
    cmd = [input_bam, bwa_db, taxonomy_db, out_report, '--dupeReport', dupe_report, '--outBam', out_bam]
    parser = metagenomics.parser_align_rna_metagenomics(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(dupe_report) > 0
    assert os.path.getsize(out_bam) > 0
