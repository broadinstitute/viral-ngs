import math
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from viral_ngs import metagenomics
from viral_ngs.classify import lyra


SCORE_RECORDS = 1_000_000
BLOCK_SIZE = 10
BLOCKS = 100_000
SCORE_STRIDE = 7_919
SCORE_OFFSET = 17
BAM_STRIDE = 104_729
BAM_OFFSET = 31

NORMALIZED_HEADER = (
    "SAMPLE_ID",
    "READ_ID",
    "LYRA_N_SCORES",
    "LYRA_PAIRING",
    "LYRA_MIN_SCORE",
    "LYRA_MAX_SCORE",
    "LYRA_THRESHOLD",
    "LYRA_CALL",
)
SUMMARY_HEADER = (
    "SAMPLE_ID",
    "LYRA_THRESHOLD",
    "LYRA_INPUT_BAM_RECORDS",
    "LYRA_ELIGIBLE_BAM_RECORDS",
    "LYRA_SCORE_RECORDS",
    "LYRA_FRAGMENTS",
    "LYRA_SINGLE_END_FRAGMENTS",
    "LYRA_COMPLETE_PAIR_FRAGMENTS",
    "LYRA_INCOMPLETE_PAIR_FRAGMENTS",
    "LYRA_VIRAL_FRAGMENT_CALLS",
    "LYRA_NONVIRAL_FRAGMENT_CALLS",
    "LYRA_OUTPUT_BAM_RECORDS",
)

_SEQUENCE = "A" * 50
_QUALITY_STRING = "I" * 50


def _record_spec(ordinal):
    block, local = divmod(ordinal, BLOCK_SIZE)
    prefix = "block-{:06d}-".format(block)
    if local == 0:
        return prefix + "single-viral", "0.9", "0", 0x4
    if local == 1:
        return prefix + "single-nonviral", "0.1", "1", 0x4
    if local == 2:
        return prefix + "single-threshold", "0.8", "0", 0x4
    if local == 3:
        return prefix + "single-below", "0.79", "1", 0x4
    if local == 4:
        return prefix + "pair-viral", "0.9", "0", 0x4 | 0x1 | 0x40
    if local == 5:
        return prefix + "pair-viral", "0.8", "0", 0x4 | 0x1 | 0x80
    if local == 6:
        return prefix + "pair-nonviral", "0.9", "0", 0x4 | 0x1 | 0x40
    if local == 7:
        return prefix + "pair-nonviral", "0.79", "1", 0x4 | 0x1 | 0x80
    if local == 8:
        return prefix + "incomplete-r1", "1", "1", 0x4 | 0x1 | 0x40
    assert local == 9
    return prefix + "incomplete-r2", "0", "1", 0x4 | 0x1 | 0x80


def _write_scores(path):
    with path.open("w", encoding="utf-8", newline="") as output:
        output.write("read_id\tscore\tcall\n")
        for index in range(SCORE_RECORDS):
            ordinal = (index * SCORE_STRIDE + SCORE_OFFSET) % SCORE_RECORDS
            read_id, score, native_call, _ = _record_spec(ordinal)
            output.write("{}\t{}\t{}\n".format(read_id, score, native_call))


def _write_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        for index in range(SCORE_RECORDS):
            ordinal = (index * BAM_STRIDE + BAM_OFFSET) % SCORE_RECORDS
            read_id, _, _, flag = _record_spec(ordinal)
            record = pysam.AlignedSegment()
            record.query_name = read_id
            record.flag = flag
            record.query_sequence = _SEQUENCE
            record.query_qualities = pysam.qualitystring_to_array(_QUALITY_STRING)
            record.reference_id = -1
            record.reference_start = -1
            output.write(record)


def _stream_normalized(path):
    value_counts = {
        "Single-end": 0,
        "Paired-complete": 0,
        "Paired-incomplete": 0,
        "Viral": 0,
        "Non-viral": 0,
    }
    previous_read_id = None
    row_count = 0
    with path.open("rb") as input_stream:
        expected_header = ("\t".join(NORMALIZED_HEADER) + "\n").encode("utf-8")
        assert input_stream.readline() == expected_header
        for raw_line in input_stream:
            assert raw_line.endswith(b"\n")
            raw_fields = raw_line[:-1].split(b"\t")
            assert len(raw_fields) == len(NORMALIZED_HEADER)
            raw_line.decode("utf-8", errors="strict")
            assert raw_fields[0] == b"scale-sample"

            read_id = raw_fields[1]
            if previous_read_id is not None:
                assert previous_read_id <= read_id
            previous_read_id = read_id

            pairing = raw_fields[3].decode("ascii")
            call = raw_fields[7].decode("ascii")
            assert pairing in value_counts
            assert call in value_counts
            value_counts[pairing] += 1
            value_counts[call] += 1
            row_count += 1
    return row_count, value_counts


def _read_summary(path):
    with path.open("rb") as input_stream:
        raw_header = input_stream.readline()
        raw_row = input_stream.readline()
        assert input_stream.readline() == b""
    assert raw_header.endswith(b"\n")
    assert raw_row.endswith(b"\n")
    header = tuple(raw_header[:-1].decode("utf-8", errors="strict").split("\t"))
    row = tuple(raw_row[:-1].decode("utf-8", errors="strict").split("\t"))
    assert header == SUMMARY_HEADER
    assert len(row) == len(SUMMARY_HEADER)
    return dict(zip(header, row))


def _stream_viral_bam(path):
    record_count = 0
    approved_roles = ("-single-viral", "-single-threshold", "-pair-viral")
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as input_bam:
        for record in input_bam.fetch(until_eof=True):
            assert record.query_name.endswith(approved_roles)
            record_count += 1
    return record_count


@pytest.mark.slow
def test_lyra_postprocess_one_million_records_full_parser(
    tmp_path,
    monkeypatch,
    capsys,
):
    assert SCORE_RECORDS == BLOCK_SIZE * BLOCKS
    assert math.gcd(SCORE_STRIDE, SCORE_RECORDS) == 1
    assert math.gcd(BAM_STRIDE, SCORE_RECORDS) == 1

    score_path = tmp_path / "read_scores.tsv"
    bam_path = tmp_path / "source.bam"
    normalized_path = tmp_path / "normalized.tsv"
    summary_path = tmp_path / "summary.tsv"
    viral_bam_path = tmp_path / "viral.bam"
    _write_scores(score_path)
    _write_bam(bam_path)

    observations = SimpleNamespace(
        database_path=None,
        fragment_iterator_requests=0,
        viral_iterator_requests=0,
        close_calls=0,
        cleaned=False,
    )
    original_store = lyra.LyraFragmentStore

    class ObservedLyraFragmentStore(original_store):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            observations.database_path = Path(self.database_path)

        def _observe_iterator(self, iterator):
            database_path = Path(self.database_path)
            assert database_path.is_file()
            assert database_path.stat().st_size > 0
            assert not isinstance(iterator, (list, set, dict, tuple))
            return iterator

        def iter_fragments(self):
            observations.fragment_iterator_requests += 1
            return self._observe_iterator(super().iter_fragments())

        def iter_viral_read_ids(self):
            observations.viral_iterator_requests += 1
            return self._observe_iterator(super().iter_viral_read_ids())

        def close(self):
            database_path = Path(self.database_path)
            try:
                return super().close()
            finally:
                observations.close_calls += 1
                observations.cleaned = not database_path.exists()

    monkeypatch.setattr(lyra, "LyraFragmentStore", ObservedLyraFragmentStore)

    args = metagenomics.full_parser().parse_args(
        [
            "lyra_postprocess",
            str(score_path),
            str(bam_path),
            "scale-sample",
            str(normalized_path),
            str(summary_path),
            str(viral_bam_path),
            "--threshold",
            "0.8",
            "--tmp_dir",
            str(tmp_path),
        ]
    )
    args.func_main(args)

    assert capsys.readouterr().out == ""
    assert observations.database_path.name == "reconciliation.sqlite3"
    assert observations.fragment_iterator_requests >= 2
    assert observations.viral_iterator_requests >= 1
    assert observations.close_calls == 1
    assert observations.cleaned
    assert not observations.database_path.exists()

    assert lyra.NORMALIZED_HEADER == NORMALIZED_HEADER
    normalized_rows, normalized_counts = _stream_normalized(normalized_path)
    assert normalized_rows == 800_000
    assert normalized_counts == {
        "Single-end": 400_000,
        "Paired-complete": 200_000,
        "Paired-incomplete": 200_000,
        "Viral": 300_000,
        "Non-viral": 500_000,
    }

    assert lyra.SUMMARY_HEADER == SUMMARY_HEADER
    summary = _read_summary(summary_path)
    assert summary == {
        "SAMPLE_ID": "scale-sample",
        "LYRA_THRESHOLD": "0.8",
        "LYRA_INPUT_BAM_RECORDS": "1000000",
        "LYRA_ELIGIBLE_BAM_RECORDS": "1000000",
        "LYRA_SCORE_RECORDS": "1000000",
        "LYRA_FRAGMENTS": "800000",
        "LYRA_SINGLE_END_FRAGMENTS": "400000",
        "LYRA_COMPLETE_PAIR_FRAGMENTS": "200000",
        "LYRA_INCOMPLETE_PAIR_FRAGMENTS": "200000",
        "LYRA_VIRAL_FRAGMENT_CALLS": "300000",
        "LYRA_NONVIRAL_FRAGMENT_CALLS": "500000",
        "LYRA_OUTPUT_BAM_RECORDS": "400000",
    }

    input_bam_records = int(summary["LYRA_INPUT_BAM_RECORDS"])
    eligible_bam_records = int(summary["LYRA_ELIGIBLE_BAM_RECORDS"])
    score_records = int(summary["LYRA_SCORE_RECORDS"])
    fragments = int(summary["LYRA_FRAGMENTS"])
    single_end = int(summary["LYRA_SINGLE_END_FRAGMENTS"])
    complete_pair = int(summary["LYRA_COMPLETE_PAIR_FRAGMENTS"])
    incomplete_pair = int(summary["LYRA_INCOMPLETE_PAIR_FRAGMENTS"])
    viral_calls = int(summary["LYRA_VIRAL_FRAGMENT_CALLS"])
    nonviral_calls = int(summary["LYRA_NONVIRAL_FRAGMENT_CALLS"])
    output_bam_records = int(summary["LYRA_OUTPUT_BAM_RECORDS"])

    assert normalized_rows == fragments
    assert eligible_bam_records == score_records
    assert score_records == single_end + incomplete_pair + 2 * complete_pair
    assert fragments == single_end + complete_pair + incomplete_pair
    assert fragments == viral_calls + nonviral_calls
    assert 0 <= output_bam_records <= input_bam_records
    assert (output_bam_records == 0) == (viral_calls == 0)

    actual_output_bam_records = _stream_viral_bam(viral_bam_path)
    assert actual_output_bam_records == output_bam_records == 400_000
