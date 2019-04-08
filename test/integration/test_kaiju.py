# Integration tests for kaiju
from builtins import super
import argparse
import binascii
import fnmatch
import gzip
import os
from os.path import join
import sys
import shutil
import tempfile

import lxml.html.clean
import pytest
from Bio import SeqIO

import metagenomics
import tools
import tools.kaiju
import tools.picard
import util.file
import test.integration.fixtures


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def find_files(root_dir, filt):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, filt):
            yield join(root, filename)


krona = pytest.fixture(scope='module')(test.integration.fixtures.krona)
krona_db = pytest.fixture(scope='module')(test.integration.fixtures.krona_db)
taxonomy_db = pytest.fixture(scope='module')(test.integration.fixtures.taxonomy_db)


@pytest.fixture(scope='module')
def fastq_to_sam():
    return tools.picard.FastqToSamTool()


@pytest.fixture(scope='module')
def sam_to_fastq():
    return tools.picard.SamToFastqTool()


@pytest.fixture(scope='module')
def kaiju():
    kaiju = tools.kaiju.Kaiju()
    kaiju.install()
    return kaiju


@pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


def input_fastq_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
    return [join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]


def input_bam_paths():
    data_dir = join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
    return join(data_dir, 'test-reads.bam')


@pytest.fixture(scope='module')
def input_bam(request, tmpdir_module, fastq_to_sam, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    if db_type == 'TestMetagenomicsSimple':
        fastqs = [join(data_dir, f) for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        bam_name = 'zaire_ebola.bam'
        bam = join(tmpdir_module, bam_name)
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
def kaiju_db(request, tmpdir_module, kaiju, taxonomy_db, db_type):
    data_dir = join(util.file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db_prefix = join(tmpdir_module, db_type)
    translated = join(tmpdir_module, db_type + '.faa')

    lib_dir = join(db_dir, 'library')
    util.file.cat(translated, find_files(db_dir, '*.faa'))

    kaiju.build(db_prefix, [translated], tax_db=taxonomy_db, translate_accessions=True)
    return db_prefix + '.fmi'


def test_kaiju(kaiju_db, krona_db, taxonomy_db, input_bam):
    out_report = util.file.mkstempfname('.report')
    out_reads = util.file.mkstempfname('.reads')
    cmd = [input_bam, kaiju_db, taxonomy_db, out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kaiju(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    assert os.path.getsize(out_report) > 0
    assert os.path.getsize(out_reads) > 0

    with util.file.open_or_gzopen(out_report) as inf:
        report_lines = [x.strip().split('\t') for x in inf.readlines()]
        report_lines = [x for x in report_lines if x]

    assert is_gz_file(out_reads)
    assert os.path.getsize(out_report) > 0

    if 'TestMetagenomicsSimple' in kaiju_db:
        zaire_found = False
        tai_found = False
        for line in report_lines:
            if len(line) < 2:
                continue
            if 'Zaire ebolavirus' in line[-2] and float(line[0]) > 40:
                zaire_found = True
            elif 'Tai Forest' in line[-2]:
                tai_found = True
        assert zaire_found
        assert not tai_found

    out_html = util.file.mkstempfname('.krona.html')
    parser = metagenomics.parser_krona(argparse.ArgumentParser())
    args = parser.parse_args(['--inputType', 'kaiju', out_report, krona_db, out_html])
    args.func_main(args)

    if 'TestMetagenomicsSimple' in kaiju_db:
      ebola_found = False
      cleaner = lxml.html.clean.Cleaner(remove_unknown_tags=False, page_structure=False)
      tree = cleaner.clean_html(lxml.html.parse(out_html))
      root_node = tree.xpath('//krona/node')[0]
      for n in root_node.iterdescendants():
          if n.get('name') == 'Zaire ebolavirus':
              if int(n.xpath('magnitude/val')[0].text) > 0:
                  ebola_found = True
      assert ebola_found
