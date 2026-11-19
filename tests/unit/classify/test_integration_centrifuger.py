# Integration tests for Centrifuger
import argparse
import csv
import glob
import os
import re
import shutil

import pysam
import pytest

from viral_ngs import metagenomics
from viral_ngs.classify import centrifuger
from viral_ngs.core import file as util_file


pytestmark = pytest.mark.skipif(
    not all(shutil.which(cmd) for cmd in (
        'centrifuger', 'centrifuger-build', 'centrifuger-quant', 'centrifuger-kreport', 'bbduk.sh',
    )),
    reason='Centrifuger executables are not installed or bbduk.sh is missing',
)


ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
CUSTOM_ADAPTER = "TTGACCTAGCAATCGTACCGATTCTAAGCATCGAT"


def _derive_bam(src, dst, transform):
    """Write the first 20 read names and their mates with transformed read data."""
    names = set()
    with pysam.AlignmentFile(src, 'rb', check_sq=False) as source, \
         pysam.AlignmentFile(dst, 'wb', template=source) as output:
        for rec in source:
            if len(names) < 20:
                names.add(rec.query_name)
            if rec.query_name not in names:
                continue
            rec.set_tags([])
            seq, quals = transform(rec)
            rec.query_sequence = seq
            rec.query_qualities = quals
            output.write(rec)


def _classify_bam(db, in_bam, out_classification, *options):
    parser = metagenomics.parser_centrifuger(argparse.ArgumentParser())
    args = parser.parse_args([db, str(in_bam), str(out_classification), '--threads', '2', *options])
    args.func_main(args)


def _classifications_by_read(path):
    with open(path, 'rt') as inf:
        return {
            row['readID']: (row['taxID'], int(row['queryLength']))
            for row in csv.DictReader(inf, delimiter='\t')
        }


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


def test_centrifuger_classify_then_kreport_on_empty_bam(centrifuger_db):
    '''End-to-end empty-BAM short-circuit verification.

    Drives the same classify -> kreport chain the viral-pipelines WDL uses,
    starting from the empty.bam fixture. Both short-circuits must fire
    end-to-end via the metagenomics CLI:

      - Centrifuger.classify() detects the empty BAM via samtools.isEmpty,
        skips invoking centrifuger, and writes a header-only TSV (single
        line of column headers, no data rows) so downstream parsers see a
        well-formed but empty file.
      - Centrifuger.kreport() then detects the header-only classification
        via _classification_is_empty(has_header=True), skips invoking
        centrifuger-kreport, and writes a truly empty (0-byte) kreport.
    '''
    empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
    out_classification = util_file.mkstempfname('.empty.centrifuger.tsv')
    out_kreport = util_file.mkstempfname('.empty.centrifuger.kreport')

    classify_cmd = [centrifuger_db, empty_bam, out_classification, '--threads', '2']
    parser = metagenomics.parser_centrifuger(argparse.ArgumentParser())
    args = parser.parse_args(classify_cmd)
    args.func_main(args)

    assert os.path.exists(out_classification)
    with open(out_classification, 'rt') as inf:
        lines = [line.rstrip('\n') for line in inf if line.strip()]
    assert lines == [
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches'
    ], 'classify on empty BAM should produce a header-only TSV with no data rows'

    kreport_cmd = [centrifuger_db, out_classification, out_kreport]
    parser = metagenomics.parser_centrifuger_kreport(argparse.ArgumentParser())
    args = parser.parse_args(kreport_cmd)
    args.func_main(args)

    assert os.path.exists(out_kreport)
    assert os.path.getsize(out_kreport) == 0, \
        'kreport on header-only classification should produce an empty report'


def test_centrifuger_default_trimming_matches_clean_reads(tmp_path, centrifuger_db, metagenomics_simple):
    src = metagenomics_simple['bam']
    clean_bam = tmp_path / 'clean.bam'
    raw_bam = tmp_path / 'raw.bam'
    name_indices = {}

    def clean_payload(rec):
        index = name_indices.setdefault(rec.query_name, len(name_indices))
        seq = rec.query_sequence
        # Do not let an added poly-G tail absorb Gs already at the end of the payload.
        if index % 2 == 0:
            seq = seq.rstrip('G')
        return seq, rec.query_qualities[:len(seq)]

    def add_contaminant(rec):
        seq, quals = clean_payload(rec)
        suffix = ADAPTER if name_indices[rec.query_name] % 2 else 'G' * 15
        return seq + suffix, list(quals) + [30] * len(suffix)

    _derive_bam(src, clean_bam, clean_payload)
    _derive_bam(src, raw_bam, add_contaminant)
    clean_tsv = tmp_path / 'clean.tsv'
    trimmed_tsv = tmp_path / 'trimmed.tsv'
    untrimmed_tsv = tmp_path / 'untrimmed.tsv'
    _classify_bam(centrifuger_db, clean_bam, clean_tsv)
    _classify_bam(centrifuger_db, raw_bam, trimmed_tsv)
    _classify_bam(centrifuger_db, raw_bam, untrimmed_tsv, '--no_trim')

    trimmed = _classifications_by_read(trimmed_tsv)
    untrimmed = _classifications_by_read(untrimmed_tsv)
    assert trimmed == _classifications_by_read(clean_tsv)
    assert any(taxid == '186538' for taxid, _ in trimmed.values())
    assert untrimmed.keys() == trimmed.keys()
    assert all(untrimmed[read_id][1] > length for read_id, (_, length) in trimmed.items())


def test_centrifuger_default_trimming_preserves_xt_clipping(tmp_path, centrifuger_db, metagenomics_simple):
    src = metagenomics_simple['bam']
    clean_bam = tmp_path / 'clean.bam'
    xt_bam = tmp_path / 'xt.bam'

    def add_xt_clipped_adapter(rec):
        seq = rec.query_sequence
        rec.set_tag('XT', len(seq) + 1, value_type='i')
        return seq + ADAPTER, list(rec.query_qualities) + [30] * len(ADAPTER)

    _derive_bam(src, clean_bam, lambda rec: (rec.query_sequence, rec.query_qualities))
    _derive_bam(src, xt_bam, add_xt_clipped_adapter)
    clean_tsv = tmp_path / 'clean.tsv'
    xt_tsv = tmp_path / 'xt.tsv'
    _classify_bam(centrifuger_db, clean_bam, clean_tsv)
    _classify_bam(centrifuger_db, xt_bam, xt_tsv)

    assert _classifications_by_read(xt_tsv) == _classifications_by_read(clean_tsv)


def test_centrifuger_trim_clip_db_replaces_bundled_adapters(tmp_path, centrifuger_db, metagenomics_simple):
    custom_fasta = tmp_path / 'custom.fasta'
    custom_fasta.write_text('>custom\n' + CUSTOM_ADAPTER + '\n')
    clean_bam = tmp_path / 'clean.bam'
    tailed_bam = tmp_path / 'tailed.bam'
    name_indices = {}

    def add_adapter(rec):
        index = name_indices.setdefault(rec.query_name, len(name_indices))
        tail = CUSTOM_ADAPTER if index % 2 else ADAPTER
        return rec.query_sequence + tail, list(rec.query_qualities) + [30] * len(tail)

    src = metagenomics_simple['bam']
    _derive_bam(src, clean_bam, lambda rec: (rec.query_sequence, rec.query_qualities))
    _derive_bam(src, tailed_bam, add_adapter)
    clean_tsv = tmp_path / 'clean.tsv'
    tailed_tsv = tmp_path / 'tailed.tsv'
    _classify_bam(centrifuger_db, clean_bam, clean_tsv, '--trim_clip_db', str(custom_fasta))
    _classify_bam(centrifuger_db, tailed_bam, tailed_tsv, '--trim_clip_db', str(custom_fasta))

    clean = _classifications_by_read(clean_tsv)
    tailed = _classifications_by_read(tailed_tsv)
    assert clean
    assert tailed.keys() == clean.keys()
    for read_id in clean:
        if name_indices[re.sub(r'/[12]$', '', read_id)] % 2:
            assert tailed[read_id] == clean[read_id]
        else:
            assert tailed[read_id][1] > clean[read_id][1]


def test_centrifuger_default_trimming_to_zero_reads_writes_header_only(
        tmp_path, centrifuger_db, metagenomics_simple):
    raw_bam = tmp_path / 'all_contaminant.bam'
    name_indices = {}

    def replace_with_contaminant(rec):
        index = name_indices.setdefault(rec.query_name, len(name_indices))
        seq = ADAPTER if index % 2 else 'G' * 40
        return seq, [30] * len(seq)

    _derive_bam(metagenomics_simple['bam'], raw_bam, replace_with_contaminant)
    out_classification = tmp_path / 'classification.tsv'
    out_kreport = tmp_path / 'kreport.tsv'
    _classify_bam(centrifuger_db, raw_bam, out_classification)
    parser = metagenomics.parser_centrifuger_kreport(argparse.ArgumentParser())
    args = parser.parse_args([centrifuger_db, str(out_classification), str(out_kreport)])
    args.func_main(args)

    assert out_classification.read_text().splitlines() == [
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches'
    ]
    assert out_kreport.stat().st_size == 0
