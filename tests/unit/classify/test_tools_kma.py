# Unit tests for KMA
import os
from unittest.mock import patch, mock_open

import pytest

from viral_ngs.classify import kma
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc


@pytest.fixture
def kma_tool():
    with patch('viral_ngs.classify.kma.shutil.which', return_value='/usr/bin/kma'):
        yield kma.KMA()


@pytest.fixture
def kma_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestKMA')
    paths = {
        'bam': os.path.join(util_file.get_test_input_path(), 'G5012.3.testreads.bam'),
        'db_prefix': os.path.join(base, 'test_db'),
        'ref_fasta': os.path.join(base, 'ref.fasta'),
    }
    return paths


def test_build_invokes_kma_index_with_expected_arguments(kma_tool, kma_inputs):
    with patch('viral_ngs.classify.kma.subprocess.check_call', autospec=True) as mock_check_call:
        kma_tool.build(kma_inputs['ref_fasta'], kma_inputs['db_prefix'], num_threads=9)

        args = mock_check_call.call_args[0][0]
        assert 'kma_index' == args[0]
        assert util_misc.list_contains(['-i', kma_inputs['ref_fasta']], args)
        assert util_misc.list_contains(['-o', kma_inputs['db_prefix']], args)
        expected_threads = str(util_misc.sanitize_thread_count(9))
        assert util_misc.list_contains(['-t', expected_threads], args)


def test_classify_single_end_from_bam(kma_tool, kma_inputs):
    """Test classify with BAM input in single-end mode"""
    mkstemp_vals = ['single.1.fastq', 'single.2.fastq', 'single.s.fastq']
    size_map = {
        'single.2.fastq': 1,
        'single.s.fastq': 10,
    }

    with patch('viral_ngs.classify.kma.file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kma.os.unlink'), \
         patch('viral_ngs.classify.kma.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.kma.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.kma.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kma.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.kma.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000)):

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        kma_tool.classify(kma_inputs['bam'], kma_inputs['db_prefix'], 'out_prefix', num_threads=3)

        args = mock_check_call.call_args[0][0]
        assert 'kma' == args[0]
        assert util_misc.list_contains(['-t_db', kma_inputs['db_prefix']], args)
        assert util_misc.list_contains(['-o', 'out_prefix'], args)
        assert util_misc.list_contains(['-i', 'single.s.fastq'], args)
        expected_threads = str(util_misc.sanitize_thread_count(3))
        assert util_misc.list_contains(['-t', expected_threads], args)
        picard.execute.assert_called_once_with(
            kma_inputs['bam'],
            'single.1.fastq',
            'single.2.fastq',
            outFastq0='single.s.fastq',
            picardOptions='clip-opts',
            JVMmemory='4G',
        )


def test_classify_paired_end_from_bam(kma_tool, kma_inputs):
    """Test classify with BAM input in paired-end mode"""
    mkstemp_vals = ['paired.1.fastq', 'paired.2.fastq', 'paired.s.fastq']
    size_map = {
        'paired.2.fastq': 10,
        'paired.s.fastq': 1,
    }

    with patch('viral_ngs.classify.kma.file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kma.os.unlink'), \
         patch('viral_ngs.classify.kma.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.kma.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.kma.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kma.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.kma.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000)):

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        kma_tool.classify(kma_inputs['bam'], kma_inputs['db_prefix'], 'out_prefix', num_threads=4)

        args = mock_check_call.call_args[0][0]
        assert 'kma' == args[0]
        assert util_misc.list_contains(['-t_db', kma_inputs['db_prefix']], args)
        assert util_misc.list_contains(['-o', 'out_prefix'], args)
        assert util_misc.list_contains(['-ipe', 'paired.1.fastq paired.2.fastq'], args)
        expected_threads = str(util_misc.sanitize_thread_count(4))
        assert util_misc.list_contains(['-t', expected_threads], args)
        picard.execute.assert_called_once_with(
            kma_inputs['bam'],
            'paired.1.fastq',
            'paired.2.fastq',
            outFastq0='paired.s.fastq',
            picardOptions='clip-opts',
            JVMmemory='4G',
        )


def test_classify_returns_early_when_bam_is_empty(kma_tool, kma_inputs):
    mkstemp_vals = ['ignored.1.fastq', 'ignored.2.fastq', 'ignored.s.fastq']

    with patch('viral_ngs.classify.kma.file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kma.os.unlink'), \
         patch('viral_ngs.classify.kma.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.kma.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.kma.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kma.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.kma.open', mock_open()) as mocked_open:

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = True

        kma_tool.classify(kma_inputs['bam'], kma_inputs['db_prefix'], 'out_prefix')

        mock_check_call.assert_not_called()
        picard.execute.assert_not_called()
        mocked_open.assert_called_once_with('out_prefix.res', 'wt')
