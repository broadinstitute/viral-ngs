import csv

import pytest

from viral_ngs.classify import lucavirus


def _read_csv(path):
    with open(str(path), "rt", newline="") as fh:
        return list(csv.reader(fh))


def _read_tsv(path):
    with open(str(path), "rt", newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_prepare_sequences_writes_lucavirus_csv_and_stats(tmp_path):
    input_fasta = tmp_path / "sequences.fasta"
    lucavirus_input_csv = tmp_path / "lucavirus_input.csv"
    stats_tsv = tmp_path / "lucavirus_prepare_stats.tsv"
    input_fasta.write_text(
        ">sequence_a description\n"
        "MKTIIALSYIFCLVFADYKDDDDK\n"
        ">sequence_b\n"
        "GAVLIPFYWSTCMNQDEKRH\n"
    )

    stats = lucavirus.prepare_sequences(
        str(input_fasta),
        str(lucavirus_input_csv),
        str(stats_tsv),
    )

    assert stats == {"n_sequences": 2, "has_lucavirus_input": True}
    assert _read_csv(lucavirus_input_csv) == [
        lucavirus.LUCAVIRUS_INPUT_HEADER,
        ["sequence_a", "prot", "MKTIIALSYIFCLVFADYKDDDDK"],
        ["sequence_b", "prot", "GAVLIPFYWSTCMNQDEKRH"],
    ]
    assert _read_tsv(stats_tsv) == [
        lucavirus.PREPARE_STATS_HEADER,
        ["2", "true"],
    ]


def test_prepare_sequences_filters_by_min_length(tmp_path):
    input_fasta = tmp_path / "sequences.fasta"
    lucavirus_input_csv = tmp_path / "lucavirus_input.csv"
    stats_tsv = tmp_path / "lucavirus_prepare_stats.tsv"
    input_fasta.write_text(
        ">short_sequence\nMKT\n>long_sequence\nMKTIIALSYIFCLVFADYKDDDDK\n"
    )

    stats = lucavirus.prepare_sequences(
        str(input_fasta),
        str(lucavirus_input_csv),
        str(stats_tsv),
        min_length=10,
    )

    assert stats == {"n_sequences": 1, "has_lucavirus_input": True}
    assert _read_csv(lucavirus_input_csv) == [
        lucavirus.LUCAVIRUS_INPUT_HEADER,
        ["long_sequence", "prot", "MKTIIALSYIFCLVFADYKDDDDK"],
    ]
    assert _read_tsv(stats_tsv) == [
        lucavirus.PREPARE_STATS_HEADER,
        ["1", "true"],
    ]


def test_prepare_sequences_min_length_can_filter_all_input(tmp_path):
    input_fasta = tmp_path / "sequences.fasta"
    lucavirus_input_csv = tmp_path / "lucavirus_input.csv"
    stats_tsv = tmp_path / "lucavirus_prepare_stats.tsv"
    input_fasta.write_text(">short_sequence\nMKT\n>another_short_sequence\nGAVL\n")

    stats = lucavirus.prepare_sequences(
        str(input_fasta),
        str(lucavirus_input_csv),
        str(stats_tsv),
        min_length=10,
    )

    assert stats == {"n_sequences": 0, "has_lucavirus_input": False}
    assert _read_csv(lucavirus_input_csv) == [lucavirus.LUCAVIRUS_INPUT_HEADER]
    assert _read_tsv(stats_tsv) == [
        lucavirus.PREPARE_STATS_HEADER,
        ["0", "false"],
    ]


def test_prepare_sequences_empty_fasta_writes_header_only_input(tmp_path):
    input_fasta = tmp_path / "empty.fasta"
    lucavirus_input_csv = tmp_path / "lucavirus_input.csv"
    stats_tsv = tmp_path / "lucavirus_prepare_stats.tsv"
    input_fasta.write_text("")

    stats = lucavirus.prepare_sequences(
        str(input_fasta),
        str(lucavirus_input_csv),
        str(stats_tsv),
    )

    assert stats == {"n_sequences": 0, "has_lucavirus_input": False}
    assert _read_csv(lucavirus_input_csv) == [lucavirus.LUCAVIRUS_INPUT_HEADER]
    assert _read_tsv(stats_tsv) == [
        lucavirus.PREPARE_STATS_HEADER,
        ["0", "false"],
    ]


def test_prepare_sequences_whitespace_only_fasta_writes_header_only_input(tmp_path):
    input_fasta = tmp_path / "empty.fasta"
    lucavirus_input_csv = tmp_path / "lucavirus_input.csv"
    stats_tsv = tmp_path / "lucavirus_prepare_stats.tsv"
    input_fasta.write_text("\n \t\n")

    stats = lucavirus.prepare_sequences(
        str(input_fasta),
        str(lucavirus_input_csv),
        str(stats_tsv),
    )

    assert stats == {"n_sequences": 0, "has_lucavirus_input": False}
    assert _read_csv(lucavirus_input_csv) == [lucavirus.LUCAVIRUS_INPUT_HEADER]
    assert _read_tsv(stats_tsv) == [
        lucavirus.PREPARE_STATS_HEADER,
        ["0", "false"],
    ]


def test_prepare_sequences_rejects_unsupported_seq_type(tmp_path):
    input_fasta = tmp_path / "sequences.fasta"
    input_fasta.write_text(">sequence_a\nACGT\n")

    with pytest.raises(ValueError, match="Unsupported LucaVirus seq_type"):
        lucavirus.prepare_sequences(
            str(input_fasta),
            str(tmp_path / "lucavirus_input.csv"),
            str(tmp_path / "lucavirus_prepare_stats.tsv"),
            seq_type="gene",
        )


def test_prepare_sequences_rejects_negative_min_length(tmp_path):
    input_fasta = tmp_path / "sequences.fasta"
    input_fasta.write_text(">sequence_a\nACGT\n")

    with pytest.raises(ValueError, match="min_length must be non-negative"):
        lucavirus.prepare_sequences(
            str(input_fasta),
            str(tmp_path / "lucavirus_input.csv"),
            str(tmp_path / "lucavirus_prepare_stats.tsv"),
            min_length=-1,
        )


def test_write_empty_predictions(tmp_path):
    output_tsv = tmp_path / "empty_predictions.tsv"

    lucavirus.write_empty_predictions(str(output_tsv))

    assert _read_tsv(output_tsv) == [lucavirus.OUTPUT_HEADER]


def test_normalize_output_accepts_header_only_output(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    output_tsv = tmp_path / "normalized.tsv"
    raw_tsv.write_text("\t".join(lucavirus.OUTPUT_HEADER) + "\n")

    lucavirus.normalize_output(str(raw_tsv), str(output_tsv), task_profile="rdrp")

    assert _read_tsv(output_tsv) == [lucavirus.OUTPUT_HEADER]


def test_normalize_output_preserves_duplicate_seq_ids(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    output_tsv = tmp_path / "normalized.tsv"
    raw_tsv.write_text(
        "\t".join(lucavirus.OUTPUT_HEADER) + "\n"
        "seq1\tMKT\t0.8\t1\tviral\n"
        "seq1\tMKT\t0.7\t2\tcapsid\n"
    )

    lucavirus.normalize_output(
        str(raw_tsv),
        str(output_tsv),
        task_profile="virus_ec4",
    )

    assert _read_tsv(output_tsv) == [
        lucavirus.OUTPUT_HEADER,
        ["seq1", "MKT", "0.8", "1", "viral"],
        ["seq1", "MKT", "0.7", "2", "capsid"],
    ]


def test_normalize_output_rejects_unknown_profile(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    raw_tsv.write_text("\t".join(lucavirus.OUTPUT_HEADER) + "\n")

    with pytest.raises(ValueError, match="Unsupported LucaVirus task profile"):
        lucavirus.normalize_output(
            str(raw_tsv),
            str(tmp_path / "normalized.tsv"),
            task_profile="not_a_profile",
        )


def test_normalize_output_rejects_zero_byte_file(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    output_tsv = tmp_path / "normalized.tsv"
    raw_tsv.write_text("")

    with pytest.raises(ValueError, match="zero bytes"):
        lucavirus.normalize_output(str(raw_tsv), str(output_tsv))


def test_normalize_output_rejects_bad_header(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    output_tsv = tmp_path / "normalized.tsv"
    raw_tsv.write_text("seq_id\tseq\tprob\n")

    with pytest.raises(ValueError, match="Unexpected LucaVirus output header"):
        lucavirus.normalize_output(str(raw_tsv), str(output_tsv))


def test_normalize_output_rejects_malformed_row(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    output_tsv = tmp_path / "normalized.tsv"
    raw_tsv.write_text("\t".join(lucavirus.OUTPUT_HEADER) + "\nseq1\tMKT\t0.8\t1\n")

    with pytest.raises(ValueError, match="line 2"):
        lucavirus.normalize_output(str(raw_tsv), str(output_tsv))


def test_normalize_output_rejects_blank_data_row(tmp_path):
    raw_tsv = tmp_path / "raw.tsv"
    output_tsv = tmp_path / "normalized.tsv"
    raw_tsv.write_text("\t".join(lucavirus.OUTPUT_HEADER) + "\n\n")

    with pytest.raises(ValueError, match="Blank LucaVirus output row at line 2"):
        lucavirus.normalize_output(str(raw_tsv), str(output_tsv))
