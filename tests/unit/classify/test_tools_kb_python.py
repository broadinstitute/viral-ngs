# Unit tests for kb_python
import os

import pytest

import classify.kb
import util.file
import util.misc


@pytest.fixture
def kb_tool(mocker):
    mocker.patch('classify.kb.shutil.which', return_value='/usr/bin/kb')
    return classify.kb.kb()


@pytest.fixture
def kb_inputs():
    base = os.path.join(util.file.get_test_input_path(), 'TestKbPython')
    paths = {
        'fastq': os.path.join(base, 'SRR12340077.2.sample.fastq.gz'),
        'index': os.path.join(base, 'palmdb.corona.idx'),
        't2g': os.path.join(base, 'palmdb_clustered_t2g.txt'),
    }
    for p in paths.values():
        assert os.path.exists(p), f'Missing expected kb test input: {p}'
    return paths


def _mock_common_io(mocker, mkstemp_paths, *, is_empty=False):
    mocker.patch('classify.kb.util.file.mkstempfname', side_effect=mkstemp_paths)
    mocker.patch('classify.kb.os.unlink')

    picard_cls = mocker.patch('classify.kb.tools.picard.SamToFastqTool', autospec=True)
    picard_cls.illumina_clipping_attribute = 'XT'
    picard = picard_cls.return_value
    picard.jvmMemDefault = '4G'
    picard.execute.return_value = None
    mocker.patch('classify.kb.tools.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts')

    samtools_cls = mocker.patch('classify.kb.tools.samtools.SamtoolsTool', autospec=True)
    samtools = samtools_cls.return_value
    samtools.isEmpty.return_value = is_empty

    mock_popen = mocker.patch('classify.kb.subprocess.Popen', autospec=True)
    mock_process = mock_popen.return_value
    mock_process.communicate.return_value = ('', '')
    mock_process.returncode = 0
    return {
        'popen': mock_popen,
        'picard_execute': picard.execute,
        'samtools': samtools,
    }


def test_build_invokes_kb_ref_with_expected_arguments(mocker, kb_tool, kb_inputs):
    mock_popen = mocker.patch('classify.kb.subprocess.Popen', autospec=True)
    mock_process = mock_popen.return_value
    mock_process.communicate.return_value = ('', '')
    mock_process.returncode = 0

    kb_tool.build(kb_inputs['fastq'], kb_inputs['index'], workflow='custom', kmer_len=27, protein=True, num_threads=9)

    args = mock_popen.call_args[0][0]
    assert ['kb', 'ref'] == args[:2]
    assert util.misc.list_contains(['-i', kb_inputs['index']], args)
    assert util.misc.list_contains(['-k', '27'], args)
    assert util.misc.list_contains(['--workflow', 'custom'], args)
    assert '--aa' in args
    expected_threads = str(util.misc.sanitize_thread_count(9))
    assert util.misc.list_contains(['-t', expected_threads], args)
    assert args[-1] == kb_inputs['fastq']


def test_classify_runs_kb_count_single_end(mocker, kb_tool, kb_inputs):
    mkstemp_vals = ['single.1.fastq', 'single.2.fastq', 'single.s.fastq']
    mocks = _mock_common_io(mocker, mkstemp_vals)

    size_map = {
        'single.2.fastq': 1,
        'single.s.fastq': 10,
    }
    mocker.patch('classify.kb.os.path.getsize', side_effect=lambda path: size_map[path])

    kb_tool.classify(kb_inputs['fastq'], kb_inputs['index'], 'out_dir', kb_inputs['t2g'], num_threads=3, loom=True)

    args = mocks['popen'].call_args[0][0]
    assert ['kb', 'count'] == args[:2]
    assert args[-1] == 'single.s.fastq'
    assert util.misc.list_contains(['--parity', 'single'], args)
    assert '--loom' in args
    assert '--h5ad' not in args
    expected_threads = str(util.misc.sanitize_thread_count(3))
    assert util.misc.list_contains(['--threads', expected_threads], args)
    mocks['picard_execute'].assert_called_once_with(
        kb_inputs['fastq'],
        'single.1.fastq',
        'single.2.fastq',
        outFastq0='single.s.fastq',
        picardOptions='clip-opts',
        JVMmemory='4G',
    )


def test_classify_returns_early_when_bam_is_empty(mocker, kb_tool, kb_inputs):
    mkstemp_vals = ['ignored.1.fastq', 'ignored.2.fastq', 'ignored.s.fastq']
    mocks = _mock_common_io(mocker, mkstemp_vals, is_empty=True)

    kb_tool.classify(kb_inputs['fastq'], kb_inputs['index'], 'out_dir', kb_inputs['t2g'])

    mocks['popen'].assert_not_called()
    mocks['picard_execute'].assert_not_called()
