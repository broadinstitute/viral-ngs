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
def genomad_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestGenomad')
    paths = {
        'fasta': os.path.join(base, 'small.fasta'),
        'db_path': base,  # Use the TestGenomad dir itself as a stand-in for a db directory
        'out_dir': '/tmp/genomad_test_output',
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
            threads=4,
        )

        mock_genomad.assert_called_once_with()
        mock_genomad.return_value.end_to_end.assert_called_once_with(
            genomad_inputs['fasta'],
            genomad_inputs['db_path'],
            genomad_inputs['out_dir'],
            num_threads=4,
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
        )
