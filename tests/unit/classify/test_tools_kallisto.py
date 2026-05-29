# Unit tests for kallisto
import csv
import gzip
import os
import subprocess
from unittest.mock import patch

import anndata
import numpy as np
import pytest
import scipy.sparse

from viral_ngs.classify import kallisto
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc


@pytest.fixture
def kallisto_tool():
    with patch('viral_ngs.classify.kallisto.shutil.which', return_value='/usr/bin/kb'):
        yield kallisto.Kallisto()


@pytest.fixture
def kallisto_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestKallisto')
    paths = {
        'fastq': os.path.join(base, 'SRR12340077.2.sample.fastq.gz'),
        'bam': os.path.join(util_file.get_test_input_path(), 'G5012.3.testreads.bam'),
        'index': os.path.join(base, 'palmdb.corona.idx'),
        't2g': os.path.join(base, 'palmdb_clustered_t2g.txt'),
    }
    for p in paths.values():
        assert os.path.exists(p), f'Missing expected kallisto test input: {p}'
    return paths


def test_build_invokes_ref_with_expected_arguments(kallisto_tool, kallisto_inputs):
    with patch('viral_ngs.classify.kallisto.subprocess.Popen', autospec=True) as mock_popen:
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kallisto_tool.build(kallisto_inputs['fastq'], kallisto_inputs['index'], workflow='custom', kmer_len=27, protein=True, num_threads=9)

        args = mock_popen.call_args[0][0]
        assert ['kb', 'ref'] == args[:2]
        assert util_misc.list_contains(['-i', kallisto_inputs['index']], args)
        assert util_misc.list_contains(['-k', '27'], args)
        assert util_misc.list_contains(['--workflow', 'custom'], args)
        assert '--aa' in args
        expected_threads = str(util_misc.sanitize_thread_count(9))
        assert util_misc.list_contains(['-t', expected_threads], args)
        assert args[-1] == kallisto_inputs['fastq']


def test_version_parses_kb_python_version_from_micromamba_json(kallisto_tool):
    package_json = (
        '[\n'
        '  {\n'
        '    "name": "kb-python",\n'
        '    "version": "0.30.2",\n'
        '    "dist_name": "kb-python-0.30.2-pyh106432d_0"\n'
        '  }\n'
        ']\n'
    )

    with patch(
        'viral_ngs.classify.kallisto.subprocess.run',
        return_value=subprocess.CompletedProcess(
            ['micromamba', 'list', 'kb-python', '--json'],
            0,
            stdout=package_json,
            stderr='',
        ),
    ) as mock_run:
        assert kallisto_tool.version() == '0.30.2'

    mock_run.assert_called_once_with(
        ['micromamba', 'list', 'kb-python', '--json'],
        capture_output=True,
        text=True,
        check=True,
    )


def test_version_returns_unknown_when_micromamba_fails(kallisto_tool):
    with patch(
        'viral_ngs.classify.kallisto.subprocess.run',
        side_effect=subprocess.CalledProcessError(1, ['micromamba', 'list', 'kb-python', '--json']),
    ):
        assert kallisto_tool.version() == 'unknown'


def test_execute_raises_when_kb_writes_traceback_with_zero_return_code(kallisto_tool):
    with patch('viral_ngs.classify.kallisto.subprocess.Popen', autospec=True) as mock_popen:
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', 'Traceback (most recent call last)\nboom\n')
        mock_process.returncode = 0

        with pytest.raises(subprocess.CalledProcessError) as exc_info:
            kallisto_tool.execute('kb count', None)

    assert exc_info.value.returncode == 1


def test_classify_runs_count_single_end_from_bam(kallisto_tool, kallisto_inputs):
    """Test classify with BAM input - should convert to FASTQ via Picard"""
    mkstemp_vals = ['single.1.fastq', 'single.2.fastq', 'single.s.fastq']
    size_map = {
        'single.2.fastq': 1,
        'single.s.fastq': 10,
    }

    with patch('viral_ngs.classify.kallisto.file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kallisto.os.unlink'), \
         patch('viral_ngs.classify.kallisto.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.kallisto.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.kallisto.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kallisto.subprocess.Popen', autospec=True) as mock_popen, \
         patch('viral_ngs.classify.kallisto.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000)), \
         patch.object(kallisto_tool, '_finalize_count_outputs') as finalize_outputs:

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kallisto_tool.classify(kallisto_inputs['bam'], kallisto_inputs['index'], 'out_dir', kallisto_inputs['t2g'], num_threads=3, loom=True)

        args = mock_popen.call_args[0][0]
        assert ['kb', 'count'] == args[:2]
        assert args[-1] == 'single.s.fastq'
        assert util_misc.list_contains(['--parity', 'single'], args)
        assert '--loom' in args
        assert '--h5ad' in args
        expected_threads = str(util_misc.sanitize_thread_count(3))
        assert util_misc.list_contains(['-t', expected_threads], args)
        finalize_outputs.assert_called_once_with('out_dir', 'G5012.3.testreads')
        picard.execute.assert_called_once_with(
            kallisto_inputs['bam'],
            'single.1.fastq',
            'single.2.fastq',
            outFastq0='single.s.fastq',
            picardOptions='clip-opts',
            JVMmemory='4G',
        )


def test_classify_runs_count_with_fastq_input(kallisto_tool, kallisto_inputs):
    """Test classify with FASTQ input - should skip Picard and use file directly"""
    with patch('viral_ngs.classify.kallisto.os.path.exists', return_value=True), \
         patch('viral_ngs.classify.kallisto.subprocess.Popen', autospec=True) as mock_popen, \
         patch('viral_ngs.classify.kallisto.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch.object(kallisto_tool, '_finalize_count_outputs') as finalize_outputs:

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kallisto_tool.classify(kallisto_inputs['fastq'], kallisto_inputs['index'], 'out_dir', kallisto_inputs['t2g'], num_threads=3)

        args = mock_popen.call_args[0][0]
        assert ['kb', 'count'] == args[:2]
        assert args[-1] == kallisto_inputs['fastq']
        assert util_misc.list_contains(['--parity', 'single'], args)
        assert '--h5ad' in args
        expected_threads = str(util_misc.sanitize_thread_count(3))
        assert util_misc.list_contains(['-t', expected_threads], args)
        finalize_outputs.assert_called_once_with('out_dir', 'SRR12340077.2.sample')


@pytest.mark.parametrize('missing_name', ['missing.bam', 'missing.fastq'])
def test_classify_raises_when_input_is_missing(kallisto_tool, kallisto_inputs, tmp_path, missing_name):
    missing_input = tmp_path / missing_name

    with pytest.raises(FileNotFoundError, match=missing_name):
        kallisto_tool.classify(str(missing_input), kallisto_inputs['index'], str(tmp_path / 'out'), kallisto_inputs['t2g'])


def test_classify_returns_early_when_bam_is_empty(kallisto_tool, kallisto_inputs, tmp_path):
    mkstemp_vals = ['ignored.1.fastq', 'ignored.2.fastq', 'ignored.s.fastq']
    out_dir = str(tmp_path / 'empty-kallisto-out')

    with patch('viral_ngs.classify.kallisto.file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kallisto.os.unlink'), \
         patch('viral_ngs.classify.kallisto.samtools.SamtoolsTool', autospec=True) as samtools_cls:

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = True

        kallisto_tool.classify(kallisto_inputs['bam'], kallisto_inputs['index'], out_dir, kallisto_inputs['t2g'])

    h5ad_path = os.path.join(out_dir, 'adata.h5ad')
    counts_tsv = os.path.join(out_dir, 'counts.tsv')
    assert os.path.exists(h5ad_path)
    assert os.path.exists(counts_tsv)

    adata = anndata.read_h5ad(h5ad_path)
    assert adata.shape == (1, 0)
    assert adata.obs_names.tolist() == ['G5012.3.testreads']
    assert adata.obs['sample'].tolist() == ['G5012.3.testreads']
    assert adata.obs['batch_name'].tolist() == ['G5012.3.testreads']

    with open(counts_tsv) as inf:
        assert inf.read() == 'sample_id\tdb_hit_id\tcount\n'


@pytest.mark.parametrize('missing_name', ['missing.bam', 'missing.fastq'])
def test_extract_raises_when_input_is_missing(kallisto_tool, kallisto_inputs, tmp_path, missing_name):
    missing_input = tmp_path / missing_name

    with pytest.raises(FileNotFoundError, match=missing_name):
        kallisto_tool.extract(
            str(missing_input),
            kallisto_inputs['index'],
            ['hit1'],
            str(tmp_path / 'out'),
            kallisto_inputs['t2g'],
        )


def test_extract_writes_empty_tsvs_when_bam_is_empty(kallisto_tool, kallisto_inputs, tmp_path):
    out_dir = str(tmp_path / 'empty-extract-out')

    with patch('viral_ngs.classify.kallisto.samtools.SamtoolsTool', autospec=True) as samtools_cls:
        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = True

        kallisto_tool.extract(kallisto_inputs['bam'], kallisto_inputs['index'], ['hit1'], out_dir, kallisto_inputs['t2g'])

    with open(os.path.join(out_dir, 'read_hits.tsv')) as inf:
        assert inf.read() == 'read_id\tdb_hit_id\n'
    with open(os.path.join(out_dir, 'summary.tsv')) as inf:
        assert inf.read() == 'SAMPLE_ID\tREAD_ID\tDB_ID\tTAXONOMY_LINEAGE\tTAXONOMY_NAME\tSEQUENCE_LENGTH\n'


def test_write_empty_count_outputs_creates_h5ad_and_header_only_counts(kallisto_tool, tmp_path):
    kallisto_tool._write_empty_count_outputs(str(tmp_path), 'empty-sample')

    h5ad_path = tmp_path / 'adata.h5ad'
    counts_tsv = tmp_path / 'counts.tsv'

    assert h5ad_path.exists()
    assert counts_tsv.exists()

    adata = anndata.read_h5ad(h5ad_path)
    assert adata.shape == (1, 0)
    assert adata.obs_names.tolist() == ['empty-sample']
    assert adata.obs['sample'].tolist() == ['empty-sample']
    assert adata.obs['batch_name'].tolist() == ['empty-sample']

    with open(counts_tsv) as inf:
        assert inf.read() == 'sample_id\tdb_hit_id\tcount\n'


def test_extract_with_no_targets_writes_empty_tsvs(kallisto_tool, kallisto_inputs, tmp_path):
    kallisto_tool.extract(kallisto_inputs['bam'], kallisto_inputs['index'], [], str(tmp_path), kallisto_inputs['t2g'])

    with open(tmp_path / 'read_hits.tsv') as inf:
        assert inf.read() == 'read_id\tdb_hit_id\n'
    with open(tmp_path / 'summary.tsv') as inf:
        assert inf.read() == 'SAMPLE_ID\tREAD_ID\tDB_ID\tTAXONOMY_LINEAGE\tTAXONOMY_NAME\tSEQUENCE_LENGTH\n'


def read_tsv(path):
    with open(path) as inf:
        return list(csv.DictReader(inf, delimiter='\t'))


def test_sample_name_from_input_strips_known_suffixes():
    assert kallisto.Kallisto._sample_name_from_input('/path/sample.fastq.gz') == 'sample'
    assert kallisto.Kallisto._sample_name_from_input('/path/sample.fq.gz') == 'sample'
    assert kallisto.Kallisto._sample_name_from_input('/path/sample.fastq') == 'sample'
    assert kallisto.Kallisto._sample_name_from_input('/path/sample.fq') == 'sample'
    assert kallisto.Kallisto._sample_name_from_input('/path/sample.bam') == 'sample'
    assert kallisto.Kallisto._sample_name_from_input('/path/sample.txt') == 'sample'


@pytest.mark.parametrize('matrix', [
    np.array([[1, 0, 2], [0, 3, 0]]),
    scipy.sparse.csr_matrix([[1, 0, 2], [0, 3, 0]]),
])
def test_write_counts_tsv_from_h5ad_writes_long_form_counts(kallisto_tool, tmp_path, matrix):
    h5ad_path = tmp_path / 'adata.h5ad'
    out_tsv = tmp_path / 'counts.tsv'
    adata = anndata.AnnData(matrix)
    adata.obs_names = ['barcode1', 'barcode2']
    adata.obs['sample'] = ['sample1', 'sample2']
    adata.var_names = ['hit1', 'hit2', 'hit3']
    adata.write_h5ad(h5ad_path)

    kallisto_tool.write_counts_tsv_from_h5ad(str(h5ad_path), str(out_tsv))

    assert read_tsv(out_tsv) == [
        {'sample_id': 'sample1', 'db_hit_id': 'hit1', 'count': '1'},
        {'sample_id': 'sample1', 'db_hit_id': 'hit3', 'count': '2'},
        {'sample_id': 'sample2', 'db_hit_id': 'hit2', 'count': '3'},
    ]


def test_write_counts_tsv_from_empty_h5ad_writes_header_only(kallisto_tool, tmp_path):
    h5ad_path = tmp_path / 'empty.h5ad'
    out_tsv = tmp_path / 'counts.tsv'
    adata = anndata.AnnData(np.zeros((0, 2)))
    adata.var_names = ['hit1', 'hit2']
    adata.write_h5ad(h5ad_path)

    kallisto_tool.write_counts_tsv_from_h5ad(str(h5ad_path), str(out_tsv))

    with open(out_tsv) as inf:
        assert inf.read() == 'sample_id\tdb_hit_id\tcount\n'


def test_extract_hit_ids_from_h5ad_rejects_multi_sample_h5ad(kallisto_tool, tmp_path):
    h5ad_path = tmp_path / 'multi_sample.h5ad'
    adata = anndata.AnnData(np.array([[1, 0], [0, 2]]))
    adata.obs_names = ['sample1', 'sample2']
    adata.var_names = ['hit1', 'hit2']
    adata.write_h5ad(h5ad_path)

    with pytest.raises(ValueError, match='Expected single-sample h5ad'):
        kallisto_tool.extract_hit_ids_from_h5ad(str(h5ad_path))


def test_write_top_taxa_report_from_counts_tsv_filters_and_ranks_focal_taxon(kallisto_tool, tmp_path):
    counts_tsv = tmp_path / 'counts.tsv'
    counts_tsv.write_text(
        'sample_id\tdb_hit_id\tcount\n'
        'sampleA\thit1\t3\n'
        'sampleA\thit2\t2\n'
        'sampleA\thit3\t5\n'
        'sampleB\thit1\t1\n'
        'sampleB\thit4\t0\n'
    )
    taxonomy_map = tmp_path / 'taxonomy.csv'
    taxonomy_map.write_text(
        'palmDB_ID,palmDB_ID,tax_level_1,tax_level_2,tax_level_3,strand\n'
        'hit1,hit1,Viruses,.,Coronaviridae,+\n'
        'hit2,hit2,Viruses,Filoviridae,.,+\n'
        'hit3,hit3,Bacteria,Firmicutes,.,+\n'
    )
    out_report = tmp_path / 'top_taxa.tsv'

    kallisto_tool.write_top_taxa_report_from_counts_tsv(
        str(counts_tsv),
        str(out_report),
        id_to_tax_map=str(taxonomy_map),
        target_taxon='Viruses',
    )

    assert read_tsv(out_report) == [
        {
            'focal_taxon_name': 'Viruses',
            'focal_taxon_count': '6',
            'palmdb_id': 'hit1',
            'hit_id': 'Coronaviridae',
            'hit_lowest_taxa_name': 'Coronaviridae',
            'hit_reads': '4',
            'pct_of_focal': '66.66666666666667',
        },
        {
            'focal_taxon_name': 'Viruses',
            'focal_taxon_count': '6',
            'palmdb_id': 'hit2',
            'hit_id': 'Filoviridae',
            'hit_lowest_taxa_name': 'Filoviridae',
            'hit_reads': '2',
            'pct_of_focal': '33.333333333333336',
        },
    ]


def test_write_top_taxa_report_from_headerless_taxonomy_map_ignores_trailing_strand(kallisto_tool, tmp_path):
    counts_tsv = tmp_path / 'counts.tsv'
    counts_tsv.write_text('sample_id\tdb_hit_id\tcount\nsampleA\tu104347\t429\n')
    taxonomy_map = tmp_path / 'taxonomy.csv'
    taxonomy_map.write_text(
        'u104347,Viruses,u1,Pisuviricota,Pisoniviricetes,Nidovirales,'
        'Coronaviridae,Betacoronavirus,'
        'Severe acute respiratory syndrome-related coronavirus,-ssRNA\n'
    )
    out_report = tmp_path / 'top_taxa.tsv'

    kallisto_tool.write_top_taxa_report_from_counts_tsv(
        str(counts_tsv),
        str(out_report),
        id_to_tax_map=str(taxonomy_map),
        target_taxon='Viruses',
    )

    assert read_tsv(out_report) == [
        {
            'focal_taxon_name': 'Viruses',
            'focal_taxon_count': '429',
            'palmdb_id': 'u104347',
            'hit_id': 'Severe acute respiratory syndrome-related coronavirus',
            'hit_lowest_taxa_name': 'Severe acute respiratory syndrome-related coronavirus',
            'hit_reads': '429',
            'pct_of_focal': '100.0',
        },
    ]


def test_load_id_to_tax_map_includes_first_headerless_row(kallisto_tool, tmp_path):
    taxonomy_map = tmp_path / 'taxonomy.csv'
    taxonomy_map.write_text(
        'u1,Viruses,u1,Pisuviricota,Pisoniviricetes,Nidovirales,'
        'Coronaviridae,Betacoronavirus,'
        'Severe acute respiratory syndrome-related coronavirus,+ssRNA\n'
    )

    assert kallisto_tool._load_id_to_tax_map(str(taxonomy_map), taxonomy_level='deepest') == {
        'u1': (
            'Viruses;u1;Pisuviricota;Pisoniviricetes;Nidovirales;Coronaviridae;'
            'Betacoronavirus;Severe acute respiratory syndrome-related coronavirus',
            'Severe acute respiratory syndrome-related coronavirus',
        )
    }


def test_write_top_taxa_report_from_counts_tsv_without_taxonomy_reports_all_hits(kallisto_tool, tmp_path):
    counts_tsv = tmp_path / 'counts.tsv'
    counts_tsv.write_text(
        'sample_id\tdb_hit_id\tcount\n'
        'sampleA\thit2\t2\n'
        'sampleA\thit1\t3\n'
    )
    out_report = tmp_path / 'top_taxa.tsv'

    kallisto_tool.write_top_taxa_report_from_counts_tsv(str(counts_tsv), str(out_report))

    assert read_tsv(out_report) == [
        {
            'focal_taxon_name': 'Viruses',
            'focal_taxon_count': '5',
            'palmdb_id': 'hit1',
            'hit_id': 'hit1',
            'hit_lowest_taxa_name': 'hit1',
            'hit_reads': '3',
            'pct_of_focal': '60.0',
        },
        {
            'focal_taxon_name': 'Viruses',
            'focal_taxon_count': '5',
            'palmdb_id': 'hit2',
            'hit_id': 'hit2',
            'hit_lowest_taxa_name': 'hit2',
            'hit_reads': '2',
            'pct_of_focal': '40.0',
        },
    ]


def test_write_top_taxa_report_from_counts_tsv_writes_zero_row_for_no_focal_hits(kallisto_tool, tmp_path):
    counts_tsv = tmp_path / 'counts.tsv'
    counts_tsv.write_text('sample_id\tdb_hit_id\tcount\nsampleA\thit1\t3\n')
    taxonomy_map = tmp_path / 'taxonomy.tsv'
    taxonomy_map.write_text('id\ttax_level_1\ttax_level_2\nhit1\tBacteria\tFirmicutes\n')
    out_report = tmp_path / 'top_taxa.tsv'

    kallisto_tool.write_top_taxa_report_from_counts_tsv(
        str(counts_tsv),
        str(out_report),
        id_to_tax_map=str(taxonomy_map),
        target_taxon='Viruses',
    )

    assert read_tsv(out_report) == [
        {
            'focal_taxon_name': 'Viruses',
            'focal_taxon_count': '0',
            'palmdb_id': '',
            'hit_id': '',
            'hit_lowest_taxa_name': '',
            'hit_reads': '0',
            'pct_of_focal': '0.0',
        },
    ]


def test_write_top_taxa_report_from_counts_tsv_validates_counts_schema(kallisto_tool, tmp_path):
    counts_tsv = tmp_path / 'bad_counts.tsv'
    counts_tsv.write_text('sample_id\tdb_hit_id\nsampleA\thit1\n')

    with pytest.raises(ValueError, match='missing required columns'):
        kallisto_tool.write_top_taxa_report_from_counts_tsv(str(counts_tsv), str(tmp_path / 'top_taxa.tsv'))


def test_write_top_taxa_report_from_counts_tsv_rejects_invalid_count(kallisto_tool, tmp_path):
    counts_tsv = tmp_path / 'bad_counts.tsv'
    counts_tsv.write_text('sample_id\tdb_hit_id\tcount\nsampleA\thit1\t1.5\n')

    with pytest.raises(ValueError, match='Non-integer count'):
        kallisto_tool.write_top_taxa_report_from_counts_tsv(str(counts_tsv), str(tmp_path / 'top_taxa.tsv'))


def test_write_extract_tsvs_writes_targeted_gz_fastq_hits(kallisto_tool, tmp_path):
    target_dir = tmp_path / 'hit1'
    target_dir.mkdir()
    with gzip.open(target_dir / '1.fastq.gz', 'wt') as outf:
        outf.write('@read one\nACGT\n+\n!!!!\n@read2\nTGCA\n+\n!!!!\n')
    unrelated_dir = tmp_path / 'metadata'
    unrelated_dir.mkdir()
    with open(unrelated_dir / 'ignored.fastq', 'wt') as outf:
        outf.write('@ignored\nACGT\n+\n!!!!\n')

    kallisto_tool._write_extract_tsvs(str(tmp_path), ['hit1'], sample_name='')

    assert read_tsv(tmp_path / 'read_hits.tsv') == [
        {'read_id': 'read', 'db_hit_id': 'hit1'},
        {'read_id': 'read2', 'db_hit_id': 'hit1'},
    ]


def test_write_extract_tsvs_writes_summary_with_taxonomy(kallisto_tool, tmp_path):
    target_dir = tmp_path / 'hit1'
    target_dir.mkdir()
    with gzip.open(target_dir / '1.fastq.gz', 'wt') as outf:
        outf.write('@read1/1\nACGT\n+\n!!!!\n@read2/2\nTGCAT\n+\n!!!!!\n')
    taxonomy_map = tmp_path / 'taxonomy.csv'
    taxonomy_map.write_text(
        'palmDB_ID,palmDB_ID,tax_level_1,tax_level_2,strand\n'
        'hit1,hit1,Viruses,Coronaviridae,+\n'
    )

    kallisto_tool._write_extract_tsvs(
        str(tmp_path),
        ['hit1'],
        sample_name='sampleA',
        id_to_tax_map=str(taxonomy_map),
        taxonomy_level='deepest',
    )

    assert read_tsv(tmp_path / 'summary.tsv') == [
        {
            'SAMPLE_ID': 'sampleA',
            'READ_ID': 'read1/1',
            'DB_ID': 'hit1',
            'TAXONOMY_LINEAGE': 'Viruses;Coronaviridae',
            'TAXONOMY_NAME': 'Coronaviridae',
            'SEQUENCE_LENGTH': '4',
        },
        {
            'SAMPLE_ID': 'sampleA',
            'READ_ID': 'read2/2',
            'DB_ID': 'hit1',
            'TAXONOMY_LINEAGE': 'Viruses;Coronaviridae',
            'TAXONOMY_NAME': 'Coronaviridae',
            'SEQUENCE_LENGTH': '5',
        },
    ]


def test_write_extract_tsvs_requires_targets(kallisto_tool, tmp_path):
    with pytest.raises(ValueError, match='target_ids must be provided'):
        kallisto_tool._write_extract_tsvs(str(tmp_path), [], sample_name='')


def test_iter_fastq_records_rejects_malformed_fastq(tmp_path):
    malformed = tmp_path / 'bad.fastq'
    malformed.write_text('not-a-header\nACGT\n+\n!!!!\n')

    with pytest.raises(ValueError, match='Expected FASTQ header'):
        list(kallisto.Kallisto._iter_fastq_records(str(malformed)))

    truncated = tmp_path / 'truncated.fastq'
    truncated.write_text('@read\nACGT\n+\n')

    with pytest.raises(ValueError, match='Unexpected end of FASTQ record'):
        list(kallisto.Kallisto._iter_fastq_records(str(truncated)))
