# Integration tests for kallisto
import argparse
import csv
import os
import platform
from pathlib import Path

import anndata
import pytest

from viral_ngs import metagenomics
from viral_ngs.core import file as util_file
from viral_ngs.core import picard

# Skip all tests on ARM platforms - kallisto is x86-only in this build configuration
IS_ARM = platform.machine() in ('arm64', 'aarch64')
pytestmark = pytest.mark.skipif(IS_ARM, reason="kallisto requires x86-only bioconda packages (not available on ARM)")


@pytest.fixture(scope='module')
def kallisto_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestKallisto')
    fastq = os.path.join(base, 'SRR12340077.2.sample.fastq.gz')
    index = os.path.join(base, 'palmdb.corona.idx')
    t2g = os.path.join(base, 'palmdb_clustered_t2g.txt')
    ref_fasta = os.path.join(base, 'palmdb_rdrp_seqs.corona.fa')
    for path in (fastq, index, t2g, ref_fasta):
        if not os.path.exists(path):
            pytest.skip(f"Required kallisto test input missing: {path}")
    return {'fastq': fastq, 'index': index, 't2g': t2g, 'ref_fasta': ref_fasta}


@pytest.fixture(scope='module')
def kallisto_bam(tmp_path_factory, kallisto_inputs):
    workdir = tmp_path_factory.mktemp('kallisto_input_bam')
    bam_path = Path(workdir) / 'input.bam'
    try:
        picard_tool = picard.FastqToSamTool()
        picard_tool.execute(kallisto_inputs['fastq'], None, 'kallisto_sample', str(bam_path))
    except Exception as exc:
        pytest.skip(f"Unable to create BAM for kallisto tests: {exc}")
    return str(bam_path)


@pytest.fixture(scope='module')
def kallisto_count_result(tmp_path_factory, kallisto_inputs, kallisto_bam):
    workdir = tmp_path_factory.mktemp('kallisto_count')
    out_dir = Path(workdir) / 'counts'
    out_dir.mkdir()
    argv = [
        '--index', kallisto_inputs['index'],
        '--t2g', kallisto_inputs['t2g'],
        '--kmer_len', '31',
        '--technology', 'bulk',
        '--out_dir', str(out_dir),
        '--threads', '4',
        kallisto_bam,
    ]
    _run_metagenomics(metagenomics.parser_kallisto, argv, cwd=str(workdir))
    h5ad_file = out_dir / 'adata.h5ad'
    counts_tsv = out_dir / 'counts.tsv'
    assert h5ad_file.exists(), f"No top-level h5ad produced at {h5ad_file}"
    assert counts_tsv.exists(), f"No counts TSV produced at {counts_tsv}"
    return {
        'workdir': Path(workdir),
        'out_dir': out_dir,
        'h5ad': h5ad_file,
        'counts_tsv': counts_tsv,
    }


@pytest.fixture(scope='module')
def kallisto_extract_result(tmp_path_factory, kallisto_inputs, kallisto_bam):
    workdir = tmp_path_factory.mktemp('kallisto_extract')
    out_dir = Path(workdir) / 'extract'
    out_dir.mkdir()
    argv = [
        '--index', kallisto_inputs['index'],
        '--t2g', kallisto_inputs['t2g'],
        '--out_dir', str(out_dir),
        '--targets', 'u100031',
        '--protein',
        '--threads', '4',
        kallisto_bam
    ]
    _run_metagenomics(metagenomics.parser_kallisto_extract, argv, cwd=str(workdir))
    extracted = out_dir / 'u100031' / '1.fastq.gz'
    read_hits = out_dir / 'read_hits.tsv'
    assert extracted.exists(), f"Expected extracted reads at {extracted}"
    assert read_hits.exists(), f"Expected read hit manifest at {read_hits}"
    count = util_file.count_fastq_reads(str(extracted))
    return {'workdir': Path(workdir), 'reads': extracted, 'read_hits': read_hits, 'count': count}


@pytest.fixture(scope='module')
def kallisto_ref_result(tmp_path_factory, kallisto_inputs):
    workdir = tmp_path_factory.mktemp('kallisto_ref')
    out_dir = Path(workdir) / 'ref'
    out_dir.mkdir()
    idx_path = out_dir / 'index.idx'
    fasta_path = Path(kallisto_inputs['ref_fasta'])
    argv = [
        '--index', str(idx_path),
        '--workflow', 'custom',
        '--kmer_len', '31',
        '--protein',
        '--threads', '2',
        str(fasta_path),
    ]
    _run_metagenomics(metagenomics.parser_kallisto_build, argv, cwd=str(workdir))
    assert idx_path.exists(), f"Index not created at {idx_path}"
    return {'workdir': Path(workdir), 'index': idx_path}


def _run_metagenomics(parser_fn, argv, cwd=None):
    parser = parser_fn(argparse.ArgumentParser())
    args = parser.parse_args(argv)
    prev_cwd = os.getcwd()
    try:
        if cwd:
            os.chdir(cwd)
        args.func_main(args)
    finally:
        if cwd:
            os.chdir(prev_cwd)


def test_kallisto_help_reports_usage():
    parser = metagenomics.parser_kallisto(argparse.ArgumentParser())
    help_str = parser.format_help()
    assert '--index' in help_str
    assert '--t2g' in help_str
    assert 'kallisto_out' in help_str


def test_kallisto_count_produces_h5ad_and_counts_tsv(kallisto_count_result):
    adata = anndata.read_h5ad(str(kallisto_count_result['h5ad']))
    assert adata.shape[0] > 0
    assert float(adata.X.sum()) >= 0
    with open(kallisto_count_result['counts_tsv']) as inf:
        rows = list(csv.DictReader(inf, delimiter='\t'))
    assert rows
    assert set(rows[0]) == {'sample_id', 'db_hit_id', 'count'}
    assert {row['sample_id'] for row in rows} == {'input'}
    assert sum(int(row['count']) for row in rows) == int(adata.X.sum())


def test_kallisto_extract_yields_expected_reads(kallisto_extract_result):
    assert kallisto_extract_result['count'] == 3
    with open(kallisto_extract_result['read_hits']) as inf:
        rows = list(csv.DictReader(inf, delimiter='\t'))
    assert len(rows) == 3
    assert {row['db_hit_id'] for row in rows} == {'u100031'}
    assert all(row['read_id'] for row in rows)


def test_kallisto_ref_builds_index(kallisto_ref_result):
    assert kallisto_ref_result['index'].stat().st_size > 0
