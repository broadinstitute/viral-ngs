# Integration tests for Centrifuger
import argparse
import glob
import os
import shutil

import pytest

from viral_ngs import metagenomics
from viral_ngs.classify import centrifuger
from viral_ngs.core import file as util_file


pytestmark = pytest.mark.skipif(
    not all(shutil.which(cmd) for cmd in ('centrifuger', 'centrifuger-build', 'centrifuger-quant')),
    reason='Centrifuger executables are not installed',
)


@pytest.fixture(scope='module')
def centrifuger_tool():
    tool = centrifuger.Centrifuger()
    tool.install()
    return tool


@pytest.fixture(scope='module')
def metagenomics_simple():
    base = os.path.join(util_file.get_test_input_path(), 'TestMetagenomicsSimple')
    db_dir = os.path.join(base, 'db')
    return {
        'bam': os.path.join(base, 'test-reads.bam'),
        'nodes': os.path.join(db_dir, 'taxonomy', 'nodes.dmp'),
        'names': os.path.join(db_dir, 'taxonomy', 'names.dmp'),
        'conversion': os.path.join(db_dir, 'library', 'krakenuniq.map'),
        'refs': sorted(glob.glob(os.path.join(
            db_dir,
            'library',
            'Viruses',
            '*',
            '*_genomic.fna',
        ))),
    }


@pytest.fixture(scope='module')
def centrifuger_db(tmpdir_module, centrifuger_tool, metagenomics_simple):
    db_prefix = os.path.join(tmpdir_module, 'centrifuger_ebola')
    cmd = [
        db_prefix,
        metagenomics_simple['nodes'],
        metagenomics_simple['names'],
        '--ref_fastas', *metagenomics_simple['refs'],
        '--conversion_table', metagenomics_simple['conversion'],
        '--build_mem', '256M',
        '--threads', '2',
    ]
    parser = metagenomics.parser_centrifuger_build(argparse.ArgumentParser())
    args = parser.parse_args(cmd)
    args.func_main(args)
    assert os.path.exists(db_prefix + '.1.cfr')
    return db_prefix


@pytest.fixture(scope='module')
def centrifuger_classification(tmpdir_module, centrifuger_db, metagenomics_simple):
    out_classification = os.path.join(tmpdir_module, 'centrifuger.tsv')
    classify_cmd = [
        centrifuger_db,
        metagenomics_simple['bam'],
        out_classification,
        '--k', '1',
        '--threads', '2',
    ]
    parser = metagenomics.parser_centrifuger(argparse.ArgumentParser())
    args = parser.parse_args(classify_cmd)
    args.func_main(args)
    return out_classification


def test_centrifuger_classify(centrifuger_classification):
    with open(centrifuger_classification, 'rt') as inf:
        rows = [line.rstrip('\n').split('\t') for line in inf if line.strip()]

    assert rows[0] == [
        'readID',
        'seqID',
        'taxID',
        'score',
        '2ndBestScore',
        'hitLength',
        'queryLength',
        'numMatches',
    ]
    assert any(row[2] == '186538' for row in rows[1:])


def test_centrifuger_quant(centrifuger_db, centrifuger_classification):
    out_quant = util_file.mkstempfname('.centrifuger-quant.tsv')
    quant_cmd = [centrifuger_db, centrifuger_classification, out_quant]
    parser = metagenomics.parser_centrifuger_quant(argparse.ArgumentParser())
    args = parser.parse_args(quant_cmd)
    args.func_main(args)

    with open(out_quant, 'rt') as inf:
        report_rows = [line.rstrip('\n').split('\t') for line in inf if line.strip()]

    assert report_rows[0] == [
        'name',
        'taxID',
        'taxRank',
        'genomeSize',
        'numReads',
        'numUniqueReads',
        'abundance',
    ]
    assert any(row[1] == '186538' for row in report_rows[1:])


def test_centrifuger_kreport(centrifuger_db, centrifuger_classification):
    out_kreport = util_file.mkstempfname('.centrifuger-kreport.tsv')
    kreport_cmd = [centrifuger_db, centrifuger_classification, out_kreport]
    parser = metagenomics.parser_centrifuger_kreport(argparse.ArgumentParser())
    args = parser.parse_args(kreport_cmd)
    args.func_main(args)

    with open(out_kreport, 'rt') as inf:
        rows = [line.rstrip('\n').split('\t') for line in inf if line.strip()]

    # Kraken-style report has 6 columns:
    # percent, clade_reads, taxon_reads, rank_code, taxID, name
    assert all(len(row) == 6 for row in rows), \
        'kreport rows should have 6 tab-delimited columns'

    # Zaire ebolavirus taxID should appear in the report
    assert any(row[4] == '186538' for row in rows), \
        'Expected Zaire ebolavirus (taxID 186538) in kreport output'

    # An unclassified row and a root row should always be present
    rank_codes = [row[3] for row in rows]
    assert 'U' in rank_codes, 'kreport should include an unclassified (U) row'
