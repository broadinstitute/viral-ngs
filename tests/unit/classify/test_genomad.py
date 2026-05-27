# Unit tests for Genomad
import argparse
import os
from unittest.mock import patch

import pytest

from viral_ngs.classify import genomad
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc


@pytest.fixture
def genomad_tool():
    with patch('viral_ngs.classify.genomad.shutil.which', return_value='/usr/bin/genomad'):
        yield genomad.Genomad()


@pytest.fixture
def genomad_inputs(tmp_path):
    base = os.path.join(util_file.get_test_input_path(), 'TestGenomad')
    paths = {
        'fasta': os.path.join(base, 'small.fasta'),
        'db_path': base,  # Use the TestGenomad dir itself as a stand-in for a db directory
        'out_dir': str(tmp_path / 'genomad_output'),
    }
    return paths


def test_end_to_end_invokes_genomad_with_correct_arguments(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True), \
         patch('viral_ngs.classify.genomad.file.mkdir_p'):

        genomad_tool.end_to_end(genomad_inputs['fasta'], genomad_inputs['db_path'], genomad_inputs['out_dir'], num_threads=8)

        mock_check_call.assert_called_once()
        args = mock_check_call.call_args[0][0]

        assert args[0] == 'genomad'
        assert args[1] == 'end-to-end'
        assert genomad_inputs['fasta'] in args
        assert genomad_inputs['db_path'] in args
        assert genomad_inputs['out_dir'] in args

        # Verify correct argument order: fasta before out_dir, out_dir before db_path
        assert args.index(genomad_inputs['fasta']) < args.index(genomad_inputs['out_dir']) < args.index(genomad_inputs['db_path'])

        # Verify threads
        assert util_misc.list_contains(['--threads', str(util_misc.sanitize_thread_count(8))], args)


def test_end_to_end_invokes_genomad_with_workflow_options(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True), \
         patch('viral_ngs.classify.genomad.file.mkdir_p'):

        genomad_tool.end_to_end(
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            cleanup=True,
            restart=True,
            enable_score_calibration=True,
            composition='virome',
            min_score=0.8,
            max_fdr=0.05,
            min_number_genes=2,
            max_uscg=3,
            splits=4,
        )

        args = mock_check_call.call_args[0][0]
        assert '--cleanup' in args
        assert '--restart' in args
        assert '--enable-score-calibration' in args
        assert util_misc.list_contains(['--composition', 'virome'], args)
        assert util_misc.list_contains(['--min-score', '0.8'], args)
        assert util_misc.list_contains(['--max-fdr', '0.05'], args)
        assert util_misc.list_contains(['--min-number-genes', '2'], args)
        assert util_misc.list_contains(['--max-uscg', '3'], args)
        assert util_misc.list_contains(['--splits', '4'], args)


def test_end_to_end_invokes_genomad_with_filter_preset(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True), \
         patch('viral_ngs.classify.genomad.file.mkdir_p'):

        genomad_tool.end_to_end(
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            filter_preset='conservative',
        )

        args = mock_check_call.call_args[0][0]
        assert '--conservative' in args
        assert '--relaxed' not in args


def test_end_to_end_rejects_filter_preset_with_explicit_filters(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True):
        with pytest.raises(ValueError, match="filter_preset cannot be combined"):
            genomad_tool.end_to_end(
                genomad_inputs['fasta'],
                genomad_inputs['db_path'],
                genomad_inputs['out_dir'],
                filter_preset='relaxed',
                min_score=0.5,
            )


def test_end_to_end_creates_output_directory(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.subprocess.check_call', autospec=True), \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True), \
         patch('viral_ngs.classify.genomad.file.mkdir_p') as mock_mkdir:

        genomad_tool.end_to_end(genomad_inputs['fasta'], genomad_inputs['db_path'], genomad_inputs['out_dir'])

        mock_mkdir.assert_called_with(genomad_inputs['out_dir'])


def test_end_to_end_raises_on_invalid_database_path(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.os.path.isdir', return_value=False):
        with pytest.raises(ValueError, match="Database path"):
            genomad_tool.end_to_end(genomad_inputs['fasta'], '/nonexistent/db', genomad_inputs['out_dir'])


def test_end_to_end_raises_on_missing_fasta(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.os.path.isfile', return_value=False), \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True):
        with pytest.raises(FileNotFoundError):
            genomad_tool.end_to_end('/nonexistent/input.fasta', genomad_inputs['db_path'], genomad_inputs['out_dir'])


def test_end_to_end_skips_subprocess_on_empty_fasta(genomad_tool, genomad_inputs):
    with patch('viral_ngs.classify.genomad.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True), \
         patch('viral_ngs.classify.genomad.file.mkdir_p'), \
         patch('viral_ngs.classify.genomad.os.path.isfile', return_value=True), \
         patch('viral_ngs.classify.genomad.os.path.getsize', return_value=0):

        genomad_tool.end_to_end(genomad_inputs['fasta'], genomad_inputs['db_path'], genomad_inputs['out_dir'])

        mock_check_call.assert_not_called()


def test_end_to_end_skips_subprocess_on_whitespace_only_fasta(
        genomad_tool, genomad_inputs, tmp_path):
    """Exercise the second branch of _is_fasta_empty.

    The size-0 path is covered above. This covers the file-with-content-but-
    no-actual-data path: file size > 0 but <= 1024 bytes and content strips
    to empty (e.g., a FASTA containing only blank lines and whitespace).
    Uses a real tmp file rather than mocks so the open_or_gzopen + strip()
    branch is actually exercised end-to-end.
    """
    whitespace_fasta = tmp_path / 'whitespace_only.fasta'
    whitespace_fasta.write_text('\n   \n\t\n   \n')

    with patch('viral_ngs.classify.genomad.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.genomad.os.path.isdir', return_value=True), \
         patch('viral_ngs.classify.genomad.file.mkdir_p'):

        genomad_tool.end_to_end(
            str(whitespace_fasta),
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
        )

        mock_check_call.assert_not_called()


def test_version_parses_genomad_output(genomad_tool):
    with patch('viral_ngs.classify.genomad.subprocess.run', autospec=True) as mock_run:
        mock_run.return_value.stdout = 'geNomad, version 1.11.2\n'
        mock_run.return_value.stderr = ''

        version = genomad_tool.version()

        assert version == '1.11.2'


def test_main_genomad_single_file(genomad_inputs):
    """Verify main_genomad delegates one FASTA to the geNomad wrapper."""
    from viral_ngs import metagenomics

    with patch('viral_ngs.metagenomics.genomad.Genomad') as mock_genomad:
        metagenomics.main_genomad(
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            cleanup=True,
            restart=True,
            filter_preset='conservative',
            enable_score_calibration=True,
            composition='metagenome',
            splits=4,
            threads=4,
        )

        mock_genomad.assert_called_once_with()
        mock_genomad.return_value.end_to_end.assert_called_once_with(
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            num_threads=4,
            cleanup=True,
            restart=True,
            filter_preset='conservative',
            enable_score_calibration=True,
            composition='metagenome',
            min_score=None,
            max_fdr=None,
            min_number_genes=None,
            max_uscg=None,
            splits=4,
        )


def test_genomad_parser_invokes_tool(genomad_inputs):
    """Drive the metagenomics CLI parser end-to-end and verify dispatch.

    Unlike test_main_genomad_single_file (which calls main_genomad
    directly), this exercises the parser -> args.func_main -> wrapper
    path, so it locks in CLI argument names, order, and types. Catches
    regressions where parser_genomad and main_genomad drift apart.
    """
    from viral_ngs import metagenomics

    with patch('viral_ngs.metagenomics.genomad.Genomad') as mock_genomad:
        argv = [
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            '--threads', '4',
            '--cleanup',
            '--restart',
            '--enableScoreCalibration',
            '--composition', 'virome',
            '--minScore', '0.8',
            '--maxFdr', '0.05',
            '--minNumberGenes', '2',
            '--maxUscg', '3',
            '--splits', '4',
        ]
        parser = metagenomics.parser_genomad(argparse.ArgumentParser())
        args = parser.parse_args(argv)
        args.func_main(args)

        mock_genomad.assert_called_once_with()
        mock_genomad.return_value.end_to_end.assert_called_once_with(
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            num_threads=4,
            cleanup=True,
            restart=True,
            filter_preset=None,
            enable_score_calibration=True,
            composition='virome',
            min_score=0.8,
            max_fdr=0.05,
            min_number_genes=2,
            max_uscg=3,
            splits=4,
        )
