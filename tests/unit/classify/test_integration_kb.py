# Integration tests for kb-python
import argparse
import os
import platform
from pathlib import Path

import anndata
import pytest

from viral_ngs import metagenomics
from viral_ngs.core import file as util_file
from viral_ngs.core import picard

# Skip all tests on ARM platforms - kallisto/kb-python are x86-only (no ARM64 builds)
IS_ARM = platform.machine() in ('arm64', 'aarch64')
pytestmark = pytest.mark.skipif(IS_ARM, reason="kallisto/kb-python require x86-only bioconda packages (not available on ARM)")


@pytest.fixture(scope='module')
def kb_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestKbPython')
    fastq = os.path.join(base, 'SRR12340077.2.sample.fastq.gz')
    index = os.path.join(base, 'palmdb.corona.idx')
    t2g = os.path.join(base, 'palmdb_clustered_t2g.txt')
    ref_fasta = os.path.join(base, 'palmdb_rdrp_seqs.corona.fa')
    for path in (fastq, index, t2g, ref_fasta):
        if not os.path.exists(path):
            pytest.skip(f"Required kb test input missing: {path}")
    return {'fastq': fastq, 'index': index, 't2g': t2g, 'ref_fasta': ref_fasta}


@pytest.fixture(scope='module')
def kb_bam(tmp_path_factory, kb_inputs):
    workdir = tmp_path_factory.mktemp('kb_input_bam')
    bam_path = Path(workdir) / 'input.bam'
    try:
        picard = picard.FastqToSamTool()
        picard.execute(kb_inputs['fastq'], None, 'kb_sample', str(bam_path))
    except Exception as exc:
        pytest.skip(f"Unable to create BAM for kb tests: {exc}")
    return str(bam_path)


@pytest.fixture(scope='module')
def kb_count_result(tmp_path_factory, kb_inputs, kb_bam):
    workdir = tmp_path_factory.mktemp('kb_count')
    out_dir = Path(workdir) / 'counts'
    out_dir.mkdir()
    argv = [
        '--index', kb_inputs['index'],
        '--t2g', kb_inputs['t2g'],
        '--kmer_len', '31',
        '--technology', 'bulk',
        '--h5ad',
        '--out_dir', str(out_dir),
        '--threads', '4',
        kb_bam,
    ]
    _run_metagenomics(metagenomics.parser_kb, argv, cwd=str(workdir))
    h5ad_dir = out_dir / 'counts_unfiltered'
    h5ad_files = list(h5ad_dir.glob('*.h5ad'))
    assert h5ad_files, f"No h5ad produced under {h5ad_dir}"
    return {
        'workdir': Path(workdir),
        'out_dir': out_dir,
        'h5ad': h5ad_files[0],
    }


@pytest.fixture(scope='module')
def kb_extract_result(tmp_path_factory, kb_inputs, kb_bam):
    workdir = tmp_path_factory.mktemp('kb_extract')
    out_dir = Path(workdir) / 'extract'
    out_dir.mkdir()
    argv = [
        '--index', kb_inputs['index'],
        '--t2g', kb_inputs['t2g'],
        '--out_dir', str(out_dir),
        '--targets', 'u100031',
        '--protein',
        '--threads', '4',
        kb_bam
    ]
    _run_metagenomics(metagenomics.parser_kb_extract, argv, cwd=str(workdir))
    extracted = out_dir / 'u100031' / '1.fastq.gz'
    assert extracted.exists(), f"Expected extracted reads at {extracted}"
    count = util_file.count_fastq_reads(str(extracted))
    return {'workdir': Path(workdir), 'reads': extracted, 'count': count}


@pytest.fixture(scope='module')
def kb_ref_result(tmp_path_factory, kb_inputs):
    workdir = tmp_path_factory.mktemp('kb_ref')
    out_dir = Path(workdir) / 'ref'
    out_dir.mkdir()
    idx_path = out_dir / 'index.idx'
    fasta_path = Path(kb_inputs['ref_fasta'])
    argv = [
        '--index', str(idx_path),
        '--workflow', 'custom',
        '--kmer_len', '31',
        '--protein',
        '--threads', '2',
        str(fasta_path),
    ]
    _run_metagenomics(metagenomics.parser_kb_build, argv, cwd=str(workdir))
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


def test_kb_help_reports_usage():
    parser = metagenomics.parser_kb(argparse.ArgumentParser())
    help_str = parser.format_help()
    assert '--index' in help_str
    assert '--t2g' in help_str
    assert 'kb_out' in help_str


def test_kb_count_produces_h5ad(kb_count_result):
    adata = anndata.read_h5ad(str(kb_count_result['h5ad']))
    assert adata.shape[0] > 0
    assert float(adata.X.sum()) >= 0


def test_kb_extract_yields_expected_reads(kb_extract_result):
    assert kb_extract_result['count'] == 3


def test_kb_ref_builds_index(kb_ref_result):
    assert kb_ref_result['index'].stat().st_size > 0
