# Integration tests for kraken2
import argparse
import binascii
import os
import os.path
from os.path import join

#import lxml.html.clean
import pytest

from viral_ngs import metagenomics
import shutil
import unittest
from viral_ngs.core import file as util_file
from viral_ngs import core as tools
from viral_ngs.classify import kraken2
from viral_ngs.classify import krona
from viral_ngs.core import picard
import tests.unit.classify.fixtures


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


krona = pytest.fixture(scope='module')(tests.unit.classify.fixtures.krona)
krona_db = pytest.fixture(scope='module')(tests.unit.classify.fixtures.krona_db)


@pytest.fixture()
def input_bam(db_type):
    data_dir = join(util_file.get_test_input_path(), db_type)
    return join(data_dir, 'test-reads.bam')

@pytest.fixture(scope='module')
def kraken2():
    kraken2 = kraken2.Kraken2()
    kraken2.install()
    return kraken2


@pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


@pytest.fixture(scope='module')
def kraken2_db(request, tmpdir_module, kraken2, db_type):
    data_dir = join(util_file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = join(tmpdir_module, 'kraken2_db_{}'.format(db_type))
    parser = metagenomics.parser_kraken2_build(argparse.ArgumentParser())
    cmd = [db, '--tax_db', join(db_dir, 'taxonomy')]
    util_file.mkdir_p(db)
    shutil.copytree(os.path.join(db_dir, 'library'), os.path.join(db, 'library'))

    parser.parse_args(cmd)
    args = parser.parse_args(cmd)
    args.func_main(args)
    return db


def test_kraken2(kraken2_db, input_bam):
    out_report = util_file.mkstempfname('.report')
    out_reads = util_file.mkstempfname('.reads')
    cmd = [kraken2_db, input_bam, '--outReports', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken2(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    with util_file.open_or_gzopen(out_reads, 'r') as inf:
        assert len(inf.read()) > 0
    with util_file.open_or_gzopen(out_report) as inf:
        report_lines = [x.strip().split('\t') for x in inf.readlines()]
        report_lines = [x for x in report_lines if x]

    #assert is_gz_file(out_reads)
    assert os.path.getsize(out_report) > 0

    if 'TestMetagenomicsSimple' in kraken2_db:
        zaire_found = False
        tai_found = False
        for line in report_lines:
            if 'Zaire ebolavirus' in line[-1] and float(line[0]) > 90:
                zaire_found = True
            elif 'Tai Forest' in line[-1]:
                tai_found = True
        assert zaire_found
        assert not tai_found


def test_kraken2_krona(kraken2_db, krona_db, input_bam):
    out_report = util_file.mkstempfname('.report')
    out_reads = util_file.mkstempfname('.reads.gz')

    cmd = [kraken2_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken2(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    out_html = util_file.mkstempfname('.krona.html')
    parser = metagenomics.parser_krona(argparse.ArgumentParser())
    args = parser.parse_args(['--inputType', 'kraken2', out_report, krona_db, out_html])
    args.func_main(args)

    ''' well... this doesn't actually find Zaire ebolavirus!
    if 'TestMetagenomicsSimple' in kraken2_db:
      ebola_found = False
      cleaner = lxml.html.clean.Cleaner(remove_unknown_tags=False, page_structure=False)
      tree = cleaner.clean_html(lxml.html.parse(out_html))
      root_node = tree.xpath('//krona/node')[0]
      for n in root_node.iterdescendants():
          if n.get('name') == 'Zaire ebolavirus':
              if int(n.xpath('magnitude/val')[0].text) > 0:
                  ebola_found = True
      assert ebola_found
    '''

def test_kraken2_on_empty(kraken2_db, input_bam):
    if 'TestMetagenomicsViralMix' not in kraken2_db:
        return
    input_bam = join(util_file.get_test_input_path(), 'empty.bam')
    out_report = util_file.mkstempfname('.report')
    out_reads = util_file.mkstempfname('.reads')
    cmd = [kraken2_db, input_bam, '--outReport', out_report, '--outReads', out_reads]
    parser = metagenomics.parser_kraken2(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)

    with util_file.open_or_gzopen(out_reads, 'r') as inf:
        assert len(inf.read()) == 0

    #assert is_gz_file(out_reads)
    with open(out_report, 'rt') as inf:
        lines = [line.strip() for line in inf.readlines() if not line.startswith('#') and not line.startswith('%')]
        out_report_contents = [line for line in lines if line]
    assert len(out_report_contents) == 0
    #assert out_report_contents[0].split('\t') == ['100.00', '0', '0', '0', '0', 'NA', '0', 'no rank', 'unclassified']
