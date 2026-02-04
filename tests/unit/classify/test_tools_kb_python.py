# Unit tests for kb_python
import os
from unittest.mock import patch, MagicMock

import pytest

from viral_ngs.classify import kb
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc


@pytest.fixture
def kb_tool():
    with patch('viral_ngs.classify.kb.shutil.which', return_value='/usr/bin/kb'):
        yield kb.kb()


@pytest.fixture
def kb_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestKbPython')
    paths = {
        'fastq': os.path.join(base, 'SRR12340077.2.sample.fastq.gz'),
        'bam': os.path.join(util_file.get_test_input_path(), 'G5012.3.testreads.bam'),
        'index': os.path.join(base, 'palmdb.corona.idx'),
        't2g': os.path.join(base, 'palmdb_clustered_t2g.txt'),
    }
    for p in paths.values():
        assert os.path.exists(p), f'Missing expected kb test input: {p}'
    return paths


def test_build_invokes_kb_ref_with_expected_arguments(kb_tool, kb_inputs):
    with patch('viral_ngs.classify.kb.subprocess.Popen', autospec=True) as mock_popen:
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kb_tool.build(kb_inputs['fastq'], kb_inputs['index'], workflow='custom', kmer_len=27, protein=True, num_threads=9)

        args = mock_popen.call_args[0][0]
        assert ['kb', 'ref'] == args[:2]
        assert util_misc.list_contains(['-i', kb_inputs['index']], args)
        assert util_misc.list_contains(['-k', '27'], args)
        assert util_misc.list_contains(['--workflow', 'custom'], args)
        assert '--aa' in args
        expected_threads = str(util_misc.sanitize_thread_count(9))
        assert util_misc.list_contains(['-t', expected_threads], args)
        assert args[-1] == kb_inputs['fastq']


def test_classify_runs_kb_count_single_end_from_bam(kb_tool, kb_inputs):
    """Test classify with BAM input - should convert to FASTQ via Picard"""
    mkstemp_vals = ['single.1.fastq', 'single.2.fastq', 'single.s.fastq']
    size_map = {
        'single.2.fastq': 1,
        'single.s.fastq': 10,
    }

    with patch('viral_ngs.classify.kb.util_file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kb.os.unlink'), \
         patch('viral_ngs.classify.kb.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.kb.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.kb.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kb.subprocess.Popen', autospec=True) as mock_popen, \
         patch('viral_ngs.classify.kb.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000)), \
         patch('viral_ngs.classify.kb.glob.glob', return_value=[]):

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kb_tool.classify(kb_inputs['bam'], kb_inputs['index'], 'out_dir', kb_inputs['t2g'], num_threads=3, loom=True)

        args = mock_popen.call_args[0][0]
        assert ['kb', 'count'] == args[:2]
        assert args[-1] == 'single.s.fastq'
        assert util_misc.list_contains(['--parity', 'single'], args)
        assert '--loom' in args
        assert '--h5ad' not in args
        expected_threads = str(util_misc.sanitize_thread_count(3))
        assert util_misc.list_contains(['-t', expected_threads], args)
        picard.execute.assert_called_once_with(
            kb_inputs['bam'],
            'single.1.fastq',
            'single.2.fastq',
            outFastq0='single.s.fastq',
            picardOptions='clip-opts',
            JVMmemory='4G',
        )


def test_classify_runs_kb_count_with_fastq_input(kb_tool, kb_inputs):
    """Test classify with FASTQ input - should skip Picard and use file directly"""
    with patch('viral_ngs.classify.kb.os.path.exists', return_value=True), \
         patch('viral_ngs.classify.kb.subprocess.Popen', autospec=True) as mock_popen, \
         patch('viral_ngs.classify.kb.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kb.glob.glob', return_value=[]):

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kb_tool.classify(kb_inputs['fastq'], kb_inputs['index'], 'out_dir', kb_inputs['t2g'], num_threads=3, h5ad=True)

        args = mock_popen.call_args[0][0]
        assert ['kb', 'count'] == args[:2]
        assert args[-1] == kb_inputs['fastq']
        assert util_misc.list_contains(['--parity', 'single'], args)
        assert '--h5ad' in args
        expected_threads = str(util_misc.sanitize_thread_count(3))
        assert util_misc.list_contains(['-t', expected_threads], args)


def test_classify_returns_early_when_bam_is_empty(kb_tool, kb_inputs):
    mkstemp_vals = ['ignored.1.fastq', 'ignored.2.fastq', 'ignored.s.fastq']

    with patch('viral_ngs.classify.kb.util_file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.kb.os.unlink'), \
         patch('viral_ngs.classify.kb.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.kb.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.kb.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.kb.subprocess.Popen', autospec=True) as mock_popen:

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = True

        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ('', '')
        mock_process.returncode = 0

        kb_tool.classify(kb_inputs['bam'], kb_inputs['index'], 'out_dir', kb_inputs['t2g'])

        mock_popen.assert_not_called()
        picard.execute.assert_not_called()
