# Unit tests for Genomad
import argparse
import csv
import importlib.util
import os
from unittest.mock import patch

import pytest

from tests.unit.classify.fixtures import contig_segment, write_contig_bam
from viral_ngs.classify import genomad
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc


requires_duckdb = pytest.mark.skipif(
    importlib.util.find_spec("duckdb") is None,
    reason="DuckDB is required for geNomad read labeling tests",
)


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


def _write_virus_summary(path, calls):
    with path.open("w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t", lineterminator="\n")
        writer.writerow(genomad.VIRUS_SUMMARY_COLUMNS)
        for seq_name, coords, score in calls:
            writer.writerow([seq_name, 1000, "linear", coords, 5, 11, score, "NA", 2, 1.5, "Viruses"])


@pytest.fixture
def one_read_bam(tmp_path):
    bam = tmp_path / "one.bam"
    write_contig_bam(bam, [contig_segment("r1")])
    return str(bam)


@requires_duckdb
def test_label_reads_by_contig_whole_contig_and_provirus_overlap(tmp_path):
    bam, summary, out = (tmp_path / name for name in ("reads.bam", "virus.tsv", "labels.tsv"))
    provirus = "NODE_2|provirus_201_400"
    _write_virus_summary(summary, [("NODE_1", "NA", 0.9), (provirus, "201-400", 0.8)])
    write_contig_bam(bam, [
        contig_segment("r1", "NODE_1", 0, mapq=1),
        contig_segment("r2", "NODE_2", 250),
        contig_segment("r3", "NODE_2", 500),
        contig_segment("r4", "NODE_3"),
        contig_segment("r5", "NODE_2", 100),
        contig_segment("r6", "NODE_2", 101),
        contig_segment("r7", "NODE_1", 0, mapq=10, flag=0x41),
        contig_segment("r7", "NODE_2", 250, mapq=60, flag=0x81),
    ])

    genomad.label_reads_by_contig(str(bam), str(summary), str(out))
    with out.open() as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        assert reader.fieldnames == genomad.READ_LABEL_COLUMNS
        labels = sorted((r["read_id"], r["contig_id"], r["seq_name"], r["mapped_well"]) for r in reader)
    assert labels == [
        ("r1", "NODE_1", "NODE_1", "false"),
        ("r2", "NODE_2", provirus, "true"),
        ("r6", "NODE_2", provirus, "true"),
        ("r7", "NODE_2", provirus, "true"),
    ]


@requires_duckdb
def test_label_reads_by_contig_zero_byte_summary_writes_header(tmp_path, one_read_bam):
    summary, out = tmp_path / "empty.tsv", tmp_path / "labels.tsv"
    summary.touch()
    genomad.label_reads_by_contig(one_read_bam, str(summary), str(out))
    assert out.read_text() == "\t".join(genomad.READ_LABEL_COLUMNS) + "\n"


@requires_duckdb
@pytest.mark.parametrize("contigs", [[], ["NODE_3"]], ids=["empty-bam", "no-virus-alignment"])
def test_label_reads_by_contig_no_matching_alignment_writes_header(tmp_path, contigs):
    bam, summary, out = (tmp_path / name for name in ("reads.bam", "virus.tsv", "labels.tsv"))
    write_contig_bam(bam, [contig_segment("r1", contig) for contig in contigs])
    _write_virus_summary(summary, [("NODE_1", "NA", 0.9)])
    genomad.label_reads_by_contig(str(bam), str(summary), str(out))
    assert out.read_text() == "\t".join(genomad.READ_LABEL_COLUMNS) + "\n"


@requires_duckdb
def test_label_reads_by_contig_rejects_missing_summary_columns(tmp_path, one_read_bam):
    summary = tmp_path / "virus.tsv"
    summary.write_text("\t".join(c for c in genomad.VIRUS_SUMMARY_COLUMNS if c != "taxonomy") + "\n")
    with pytest.raises(ValueError, match="missing required columns: taxonomy"):
        genomad.label_reads_by_contig(one_read_bam, str(summary), str(tmp_path / "out.tsv"))


@requires_duckdb
def test_label_reads_by_contig_rejects_duplicate_seq_name(tmp_path, one_read_bam):
    summary = tmp_path / "virus.tsv"
    _write_virus_summary(summary, [("NODE_1", "NA", 0.9), ("NODE_1", "NA", 0.8)])
    with pytest.raises(ValueError, match="duplicate seq_name: NODE_1"):
        genomad.label_reads_by_contig(one_read_bam, str(summary), str(tmp_path / "out.tsv"))


def test_genomad_label_reads_by_contig_parser_invokes_function():
    from viral_ngs import metagenomics

    with patch.object(genomad, "label_reads_by_contig", autospec=True) as label:
        parser = metagenomics.parser_genomad_label_reads_by_contig(argparse.ArgumentParser())
        args = parser.parse_args(["in.bam", "virus.tsv", "out.tsv", "--min-mapq", "10"])
        args.func_main(args)
        label.assert_called_once_with(
            "in.bam", "virus.tsv", "out.tsv", min_mapq=10, min_identity=90.0,
            min_query_cov=80.0, duckdb_memory_limit=None, work_dir=None,
        )
