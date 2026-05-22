# Unit tests for Centrifuger
import os
from unittest.mock import mock_open, patch

import pytest

from viral_ngs.classify import centrifuger
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc


@pytest.fixture
def centrifuger_tool():
    with patch('viral_ngs.classify.centrifuger.shutil.which', return_value='/usr/bin/centrifuger'):
        yield centrifuger.Centrifuger()


@pytest.fixture
def centrifuger_inputs():
    base = os.path.join(util_file.get_test_input_path(), 'TestMetagenomicsSimple')
    db_dir = os.path.join(base, 'db')
    return {
        'bam': os.path.join(base, 'test-reads.bam'),
        'db_prefix': 'centrifuger_db',
        'nodes': os.path.join(db_dir, 'taxonomy', 'nodes.dmp'),
        'names': os.path.join(db_dir, 'taxonomy', 'names.dmp'),
        'conversion': os.path.join(db_dir, 'library', 'krakenuniq.map'),
        'ref1': os.path.join(
            db_dir,
            'library',
            'Viruses',
            'Zaire_ebolavirus',
            'GCF_000848505.1_ViralProj14703_genomic.fna',
        ),
        'ref2': os.path.join(
            db_dir,
            'library',
            'Viruses',
            'Sudan_ebolavirus',
            'GCF_000855585.1_ViralProj15012_genomic.fna',
        ),
    }


def test_build_invokes_centrifuger_build_with_expected_arguments(
        centrifuger_tool, centrifuger_inputs):
    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call:
        centrifuger_tool.build(
            centrifuger_inputs['db_prefix'],
            centrifuger_inputs['nodes'],
            centrifuger_inputs['names'],
            ref_fastas=[centrifuger_inputs['ref1'], centrifuger_inputs['ref2']],
            conversion_table=centrifuger_inputs['conversion'],
            build_mem='256M',
            num_threads=9,
        )

        args = mock_check_call.call_args[0][0]
        assert args[0] == 'centrifuger-build'
        assert util_misc.list_contains(['-o', centrifuger_inputs['db_prefix']], args)
        assert util_misc.list_contains(['--taxonomy-tree', centrifuger_inputs['nodes']], args)
        assert util_misc.list_contains(['--name-table', centrifuger_inputs['names']], args)
        assert util_misc.list_contains(['-r', centrifuger_inputs['ref1']], args)
        assert util_misc.list_contains(['-r', centrifuger_inputs['ref2']], args)
        assert util_misc.list_contains(['--conversion-table', centrifuger_inputs['conversion']], args)
        assert util_misc.list_contains(['--build-mem', '256M'], args)
        expected_threads = str(util_misc.sanitize_thread_count(9))
        assert util_misc.list_contains(['-t', expected_threads], args)


def test_build_requires_exactly_one_reference_input(
        centrifuger_tool, centrifuger_inputs):
    with pytest.raises(ValueError, match='exactly one'):
        centrifuger_tool.build(
            centrifuger_inputs['db_prefix'],
            centrifuger_inputs['nodes'],
            centrifuger_inputs['names'],
        )

    with pytest.raises(ValueError, match='exactly one'):
        centrifuger_tool.build(
            centrifuger_inputs['db_prefix'],
            centrifuger_inputs['nodes'],
            centrifuger_inputs['names'],
            ref_fastas=[centrifuger_inputs['ref1']],
            ref_list='refs.txt',
        )


def test_classify_single_end_from_bam(centrifuger_tool, centrifuger_inputs):
    mkstemp_vals = ['single.1.fastq', 'single.2.fastq', 'single.s.fastq']
    size_map = {
        'single.2.fastq': 1,
        'single.s.fastq': 10,
    }

    with patch('viral_ngs.classify.centrifuger.util_file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.centrifuger.os.path.exists', return_value=True), \
         patch('viral_ngs.classify.centrifuger.os.unlink'), \
         patch('viral_ngs.classify.centrifuger.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.centrifuger.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.centrifuger.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.centrifuger.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000)), \
         patch('viral_ngs.classify.centrifuger.open', mock_open(), create=True) as mocked_open:

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        centrifuger_tool.classify(
            centrifuger_inputs['bam'],
            centrifuger_inputs['db_prefix'],
            'out.tsv',
            k=3,
            unclassified_prefix='unclassified',
            classified_prefix='classified',
            min_hitlen=17,
            hitk_factor=5,
            merge_readpair=True,
            num_threads=4,
        )

        args = mock_check_call.call_args[0][0]
        assert args[0] == 'centrifuger'
        assert util_misc.list_contains(['-x', centrifuger_inputs['db_prefix']], args)
        assert util_misc.list_contains(['-u', 'single.s.fastq'], args)
        assert '-1' not in args
        assert '-2' not in args
        assert util_misc.list_contains(['-k', '3'], args)
        assert util_misc.list_contains(['--un', 'unclassified'], args)
        assert util_misc.list_contains(['--cl', 'classified'], args)
        assert util_misc.list_contains(['--min-hitlen', '17'], args)
        assert util_misc.list_contains(['--hitk-factor', '5'], args)
        assert '--merge-readpair' in args
        expected_threads = str(util_misc.sanitize_thread_count(4))
        assert util_misc.list_contains(['-t', expected_threads], args)
        mocked_open.assert_called_once_with('out.tsv', 'wt')
        picard.execute.assert_called_once_with(
            centrifuger_inputs['bam'],
            'single.1.fastq',
            'single.2.fastq',
            outFastq0='single.s.fastq',
            picardOptions='clip-opts',
            JVMmemory='4G',
        )


def test_classify_paired_end_from_bam(centrifuger_tool, centrifuger_inputs):
    mkstemp_vals = ['paired.1.fastq', 'paired.2.fastq', 'paired.s.fastq']
    size_map = {
        'paired.2.fastq': 10,
        'paired.s.fastq': 1,
    }

    with patch('viral_ngs.classify.centrifuger.util_file.mkstempfname', side_effect=mkstemp_vals), \
         patch('viral_ngs.classify.centrifuger.os.path.exists', return_value=True), \
         patch('viral_ngs.classify.centrifuger.os.unlink'), \
         patch('viral_ngs.classify.centrifuger.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.centrifuger.picard.PicardTools.dict_to_picard_opts', return_value='clip-opts'), \
         patch('viral_ngs.classify.centrifuger.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.centrifuger.os.path.getsize', side_effect=lambda path: size_map.get(path, 1000)), \
         patch('viral_ngs.classify.centrifuger.open', mock_open(), create=True):

        picard_cls.illumina_clipping_attribute = 'XT'
        picard = picard_cls.return_value
        picard.jvmMemDefault = '4G'
        picard.execute.return_value = None

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = False

        centrifuger_tool.classify(
            centrifuger_inputs['bam'],
            centrifuger_inputs['db_prefix'],
            'out.tsv',
            num_threads=2,
        )

        args = mock_check_call.call_args[0][0]
        assert args[0] == 'centrifuger'
        assert util_misc.list_contains(['-x', centrifuger_inputs['db_prefix']], args)
        assert util_misc.list_contains(['-1', 'paired.1.fastq'], args)
        assert util_misc.list_contains(['-2', 'paired.2.fastq'], args)
        assert '-u' not in args
        assert util_misc.list_contains(['-k', '1'], args)
        expected_threads = str(util_misc.sanitize_thread_count(2))
        assert util_misc.list_contains(['-t', expected_threads], args)
        picard.execute.assert_called_once_with(
            centrifuger_inputs['bam'],
            'paired.1.fastq',
            'paired.2.fastq',
            outFastq0='paired.s.fastq',
            picardOptions='clip-opts',
            JVMmemory='4G',
        )


def test_classify_returns_early_when_bam_is_empty(
        centrifuger_tool, centrifuger_inputs):
    with patch('viral_ngs.classify.centrifuger.samtools.SamtoolsTool', autospec=True) as samtools_cls, \
         patch('viral_ngs.classify.centrifuger.picard.SamToFastqTool', autospec=True) as picard_cls, \
         patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.centrifuger.open', mock_open(), create=True) as mocked_open:

        samtools = samtools_cls.return_value
        samtools.isEmpty.return_value = True

        centrifuger_tool.classify(
            centrifuger_inputs['bam'],
            centrifuger_inputs['db_prefix'],
            'out.tsv',
        )

        mock_check_call.assert_not_called()
        picard_cls.assert_not_called()
        mocked_open.assert_called_once_with('out.tsv', 'wt')
        handle = mocked_open.return_value.__enter__.return_value
        handle.write.assert_called_once_with(
            'readID\tseqID\ttaxID\tscore\t2ndBestScore\t'
            'hitLength\tqueryLength\tnumMatches\n'
        )


def test_classify_missing_bam_raises_file_not_found(
        centrifuger_tool, centrifuger_inputs):
    with patch('viral_ngs.classify.centrifuger.samtools.SamtoolsTool', autospec=True) as samtools_cls:
        with pytest.raises(FileNotFoundError):
            centrifuger_tool.classify(
                'missing.bam',
                centrifuger_inputs['db_prefix'],
                'out.tsv',
            )

        samtools_cls.assert_not_called()


def test_quant_invokes_centrifuger_quant_with_expected_arguments(
        centrifuger_tool, centrifuger_inputs):
    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.centrifuger.open', mock_open(), create=True) as mocked_open:
        centrifuger_tool.quant(
            centrifuger_inputs['db_prefix'],
            'classification.tsv',
            'quant.tsv',
            min_score=10,
            min_length=50,
            output_format=1,
        )

        args = mock_check_call.call_args[0][0]
        assert args[0] == 'centrifuger-quant'
        assert util_misc.list_contains(['-x', centrifuger_inputs['db_prefix']], args)
        assert util_misc.list_contains(['-c', 'classification.tsv'], args)
        assert util_misc.list_contains(['--min-score', '10'], args)
        assert util_misc.list_contains(['--min-length', '50'], args)
        assert util_misc.list_contains(['--output-format', '1'], args)
        mocked_open.assert_called_once_with('quant.tsv', 'wt')


def test_classification_to_kraken2_maps_classified_rows(tmp_path):
    classification = tmp_path / 'classification.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        'read1\tNC_002549.1\t186538\t8692\t0\t160\t204\t1\n'
        'read2\tunclassified\t0\t0\t0\t0\t150\t1\n'
        'read3\tNC_006432.1\t186540\t1234\t0\t95\t100\t1\n'
    )
    output = tmp_path / 'kraken2.tsv'

    centrifuger.Centrifuger.classification_to_kraken2(
        str(classification),
        str(output),
    )

    assert output.read_text() == (
        'C\tread1\t186538\t204\tN/A\n'
        'C\tread3\t186540\t100\tN/A\n'
    )


def test_classification_to_kraken2_writes_gzipped_output(tmp_path):
    classification = tmp_path / 'classification.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        'read1\tNC_002549.1\t186538\t8692\t0\t160\t204\t1\n'
    )
    output = tmp_path / 'kraken2.tsv.gz'

    centrifuger.Centrifuger.classification_to_kraken2(
        str(classification),
        str(output),
    )

    with util_file.open_or_gzopen(str(output), 'rt') as inf:
        assert inf.read() == 'C\tread1\t186538\t204\tN/A\n'


def test_classification_to_kraken2_writes_zstd_output(tmp_path):
    classification = tmp_path / 'classification.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        'read1\tNC_002549.1\t186538\t8692\t0\t160\t204\t1\n'
    )
    output = tmp_path / 'kraken2.tsv.zst'

    centrifuger.Centrifuger.classification_to_kraken2(
        str(classification),
        str(output),
    )

    with util_file.open_or_gzopen(str(output), 'rt') as inf:
        assert inf.read() == 'C\tread1\t186538\t204\tN/A\n'


def test_classification_to_kraken2_header_only_writes_empty_output(tmp_path):
    classification = tmp_path / 'classification.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
    )
    output = tmp_path / 'kraken2.tsv'

    centrifuger.Centrifuger.classification_to_kraken2(
        str(classification),
        str(output),
    )

    assert output.read_text() == ''


def test_classification_to_kraken2_rejects_malformed_rows(tmp_path):
    classification = tmp_path / 'classification.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        'read1\tNC_002549.1\t186538\n'
    )
    output = tmp_path / 'kraken2.tsv'

    with pytest.raises(ValueError, match='Malformed Centrifuger classification row'):
        centrifuger.Centrifuger.classification_to_kraken2(
            str(classification),
            str(output),
        )


def test_classification_to_kraken2_rejects_duplicate_classified_read_ids(tmp_path):
    classification = tmp_path / 'classification.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        'read1\tNC_002549.1\t186538\t8692\t8692\t160\t204\t2\n'
        'read1\tNC_006432.1\t186540\t8692\t8692\t160\t204\t2\n'
    )
    output = tmp_path / 'kraken2.tsv'

    with pytest.raises(ValueError, match='Expected k=1 input'):
        centrifuger.Centrifuger.classification_to_kraken2(
            str(classification),
            str(output),
        )


def test_kreport_invokes_centrifuger_kreport_with_expected_arguments(
        centrifuger_tool, centrifuger_inputs):
    # Bypass the empty-classification short-circuit (covered by its own tests
    # below); this test exercises CLI argv construction on a non-empty input.
    with patch.object(centrifuger.Centrifuger, '_classification_is_empty', return_value=False), \
         patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.centrifuger.open', mock_open(), create=True) as mocked_open:
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            'classification.tsv',
            'kreport.tsv',
            no_lca=True,
            show_zeros=True,
            is_count_table=True,
            min_score=10,
            min_length=50,
            report_score_data=True,
        )

        args = mock_check_call.call_args[0][0]
        assert args[0] == 'centrifuger-kreport'
        assert util_misc.list_contains(['-x', centrifuger_inputs['db_prefix']], args)
        assert '--no-lca' in args
        assert '--show-zeros' in args
        assert '--is-count-table' in args
        assert '--report-score-data' in args
        assert util_misc.list_contains(['--min-score', '10'], args)
        assert util_misc.list_contains(['--min-length', '50'], args)
        # classification file is a positional arg, appended after options
        assert args[-1] == 'classification.tsv'
        mocked_open.assert_called_once_with('kreport.tsv', 'wt')


def test_kreport_omits_unset_optional_flags(centrifuger_tool, centrifuger_inputs):
    with patch.object(centrifuger.Centrifuger, '_classification_is_empty', return_value=False), \
         patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call, \
         patch('viral_ngs.classify.centrifuger.open', mock_open(), create=True):
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            'classification.tsv',
            'kreport.tsv',
        )

        args = mock_check_call.call_args[0][0]
        assert args[0] == 'centrifuger-kreport'
        assert util_misc.list_contains(['-x', centrifuger_inputs['db_prefix']], args)
        assert args[-1] == 'classification.tsv'
        for flag in (
            '--no-lca', '--show-zeros', '--is-count-table',
            '--report-score-data', '--min-score', '--min-length',
        ):
            assert flag not in args, f'{flag} unexpectedly present with default kwargs'


def test_kreport_short_circuits_on_empty_classification(
        tmp_path, centrifuger_tool, centrifuger_inputs):
    empty_classification = tmp_path / 'empty.tsv'
    empty_classification.touch()
    output = tmp_path / 'kreport.tsv'

    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call:
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            str(empty_classification),
            str(output),
        )

    mock_check_call.assert_not_called()
    assert output.exists()
    assert output.read_text() == ''


def test_kreport_short_circuits_on_header_only_classification(
        tmp_path, centrifuger_tool, centrifuger_inputs):
    header_only = tmp_path / 'header_only.tsv'
    header_only.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
    )
    output = tmp_path / 'kreport.tsv'

    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call:
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            str(header_only),
            str(output),
        )

    mock_check_call.assert_not_called()
    assert output.exists()
    assert output.read_text() == ''


def test_kreport_short_circuits_on_header_and_blank_lines(
        tmp_path, centrifuger_tool, centrifuger_inputs):
    header_only = tmp_path / 'header_only.tsv'
    header_only.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        '\n'
        '   \n'
    )
    output = tmp_path / 'kreport.tsv'

    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call:
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            str(header_only),
            str(output),
        )

    mock_check_call.assert_not_called()
    assert output.exists()
    assert output.read_text() == ''


def test_kreport_runs_when_classification_has_data(
        tmp_path, centrifuger_tool, centrifuger_inputs):
    classification = tmp_path / 'with_data.tsv'
    classification.write_text(
        'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
        'r1\tNC_002549.1\t186538\t8692\t0\t160\t204\t1\n'
    )
    output = tmp_path / 'kreport.tsv'

    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call:
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            str(classification),
            str(output),
        )

    mock_check_call.assert_called_once()
    args = mock_check_call.call_args[0][0]
    assert args[0] == 'centrifuger-kreport'
    assert args[-1] == str(classification)


def test_kreport_runs_for_single_row_count_table(
        tmp_path, centrifuger_tool, centrifuger_inputs):
    count_table = tmp_path / 'count_table.tsv'
    count_table.write_text('186538\t7\n')
    output = tmp_path / 'kreport.tsv'

    with patch('viral_ngs.classify.centrifuger.subprocess.check_call', autospec=True) as mock_check_call:
        centrifuger_tool.kreport(
            centrifuger_inputs['db_prefix'],
            str(count_table),
            str(output),
            is_count_table=True,
        )

    mock_check_call.assert_called_once()
    args = mock_check_call.call_args[0][0]
    assert args[0] == 'centrifuger-kreport'
    assert '--is-count-table' in args
    assert args[-1] == str(count_table)
