# Unit tests for KMA
import os

import pytest

import classify.kma
import util.file
import util.misc


@pytest.fixture
def kma_tool(mocker):
    mocker.patch('classify.kma.shutil.which', return_value='/usr/bin/kma')
    return classify.kma.KMA()


@pytest.fixture
def kma_inputs():
    base = os.path.join(util.file.get_test_input_path(), 'TestKMA')
    paths = {
        'bam': os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam'),
        'db_prefix': os.path.join(base, 'test_db'),
        'ref_fasta': os.path.join(base, 'ref.fasta'),
    }
    return paths


def _mock_common_io(mocker, mkstemp_paths, *, is_empty=False):
    mocker.patch('classify.kma.util.file.mkstempfname', side_effect=mkstemp_paths)
    mocker.patch('classify.kma.os.unlink')

    picard_cls = mocker.patch('classify.kma.tools.picard.SamToFastqTool', autospec=True)
    picard_cls.illumina_clipping_attribute = 'XT'
    picard = picard_cls.return_value
    picard.jvmMemDefault = '4G'
    picard.execute.return_value = None
    mocker.patch('classify.kma.tools.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts')

    samtools_cls = mocker.patch('classify.kma.tools.samtools.SamtoolsTool', autospec=True)
    samtools = samtools_cls.return_value
    samtools.isEmpty.return_value = is_empty

    mock_check_call = mocker.patch('classify.kma.subprocess.check_call', autospec=True)

    return {
        'check_call': mock_check_call,
        'picard_execute': picard.execute,
        'samtools': samtools,
    }


def test_build_invokes_kma_index_with_expected_arguments(mocker, kma_tool, kma_inputs):
    mock_check_call = mocker.patch('classify.kma.subprocess.check_call', autospec=True)

    kma_tool.build(kma_inputs['ref_fasta'], kma_inputs['db_prefix'], num_threads=9)

    args = mock_check_call.call_args[0][0]
    assert 'kma_index' == args[0]
    assert util.misc.list_contains(['-i', kma_inputs['ref_fasta']], args)
    assert util.misc.list_contains(['-o', kma_inputs['db_prefix']], args)
    expected_threads = str(util.misc.sanitize_thread_count(9))
    assert util.misc.list_contains(['-t', expected_threads], args)


def test_classify_single_end_from_bam(mocker, kma_tool, kma_inputs):
    """Test classify with BAM input in single-end mode"""
    mkstemp_vals = ['single.1.fastq', 'single.2.fastq', 'single.s.fastq']
    mocks = _mock_common_io(mocker, mkstemp_vals)

    size_map = {
        'single.2.fastq': 1,
        'single.s.fastq': 10,
    }
    mocker.patch('classify.kma.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000))

    kma_tool.classify(kma_inputs['bam'], kma_inputs['db_prefix'], 'out_prefix', num_threads=3)

    args = mocks['check_call'].call_args[0][0]
    assert 'kma' == args[0]
    assert util.misc.list_contains(['-t_db', kma_inputs['db_prefix']], args)
    assert util.misc.list_contains(['-o', 'out_prefix'], args)
    assert util.misc.list_contains(['-i', 'single.s.fastq'], args)
    expected_threads = str(util.misc.sanitize_thread_count(3))
    assert util.misc.list_contains(['-t', expected_threads], args)
    mocks['picard_execute'].assert_called_once_with(
        kma_inputs['bam'],
        'single.1.fastq',
        'single.2.fastq',
        outFastq0='single.s.fastq',
        picardOptions='clip-opts',
        JVMmemory='4G',
    )


def test_classify_paired_end_from_bam(mocker, kma_tool, kma_inputs):
    """Test classify with BAM input in paired-end mode"""
    mkstemp_vals = ['paired.1.fastq', 'paired.2.fastq', 'paired.s.fastq']
    mocks = _mock_common_io(mocker, mkstemp_vals)

    size_map = {
        'paired.2.fastq': 10,
        'paired.s.fastq': 1,
    }
    mocker.patch('classify.kma.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000))

    kma_tool.classify(kma_inputs['bam'], kma_inputs['db_prefix'], 'out_prefix', num_threads=4)

    args = mocks['check_call'].call_args[0][0]
    assert 'kma' == args[0]
    assert util.misc.list_contains(['-t_db', kma_inputs['db_prefix']], args)
    assert util.misc.list_contains(['-o', 'out_prefix'], args)
    assert util.misc.list_contains(['-ipe', 'paired.1.fastq paired.2.fastq'], args)
    expected_threads = str(util.misc.sanitize_thread_count(4))
    assert util.misc.list_contains(['-t', expected_threads], args)
    mocks['picard_execute'].assert_called_once_with(
        kma_inputs['bam'],
        'paired.1.fastq',
        'paired.2.fastq',
        outFastq0='paired.s.fastq',
        picardOptions='clip-opts',
        JVMmemory='4G',
    )


def test_classify_returns_early_when_bam_is_empty(mocker, kma_tool, kma_inputs):
    mkstemp_vals = ['ignored.1.fastq', 'ignored.2.fastq', 'ignored.s.fastq']
    mocks = _mock_common_io(mocker, mkstemp_vals, is_empty=True)

    mock_open = mocker.patch('classify.kma.open', mocker.mock_open())

    kma_tool.classify(kma_inputs['bam'], kma_inputs['db_prefix'], 'out_prefix')

    mocks['check_call'].assert_not_called()
    mocks['picard_execute'].assert_not_called()
    mock_open.assert_called_once_with('out_prefix.res', 'wt')
