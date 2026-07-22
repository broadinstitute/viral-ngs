from contextlib import contextmanager
from decimal import Decimal

import pytest
import pysam

import viral_ngs.core.file as util_file
from viral_ngs.classify import lyra


def _segment(
    query_name,
    *,
    flag=0x4,
    query_length=50,
    sequence_present=True,
):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = flag
    if sequence_present:
        record.query_sequence = "A" * query_length
        record.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    record.reference_id = -1
    record.reference_start = -1
    return record


def _write_bam(tmp_path, records, name="reads.bam"):
    path = tmp_path / name
    with pysam.AlignmentFile(
        str(path),
        "wb",
        header={"HD": {"VN": "1.6"}},
    ) as bam:
        for record in records:
            bam.write(record)
    return path


def _write_score_table(tmp_path, rows, name="scores.tsv"):
    path = tmp_path / name
    path.write_text(
        "read_id\tscore\tcall\n"
        + "".join("{}\t{}\t{}\n".format(*row) for row in rows),
        encoding="utf-8",
    )
    return path


@contextmanager
def _artifact_store(
    tmp_path,
    score_rows,
    bam_records,
    *,
    sample_id="sample",
    threshold="0.8",
    name="artifacts",
):
    score_path = _write_score_table(tmp_path, score_rows, name=name + ".tsv")
    bam_path = _write_bam(tmp_path, bam_records, name=name + ".bam")
    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        sample_id,
        threshold,
        work_dir=tmp_path,
    ) as store:
        yield store, bam_path


def _representative_store(tmp_path):
    return _artifact_store(
        tmp_path,
        [
            ("réad", "1.000", "0"),
            ("read-2", "0.1000", "1"),
            ("read-10", "0.9000", "0"),
            ("read-10", "0.8000", "1"),
        ],
        [
            _segment("read-2"),
            _segment("read-10", flag=0x4 | 0x1 | 0x80),
            _segment("réad", flag=0x4 | 0x1 | 0x40),
            _segment("read-10", flag=0x4 | 0x1 | 0x40),
        ],
        sample_id="échantillon exact",
        threshold=Decimal("8e-1"),
        name="representative",
    )


def test_normalized_output_has_exact_utf8_lf_schema_order_and_values(tmp_path):
    output = tmp_path / "normalized.tsv"
    with _representative_store(tmp_path) as (store, _):
        assert store.threshold == Decimal("8e-1")
        assert lyra._write_normalized(store, output) == 3

    expected = (
        "SAMPLE_ID\tREAD_ID\tLYRA_N_SCORES\tLYRA_PAIRING\t"
        "LYRA_MIN_SCORE\tLYRA_MAX_SCORE\tLYRA_THRESHOLD\tLYRA_CALL\n"
        "échantillon exact\tread-10\t2\tPaired-complete\t0.8\t0.9\t0.8\tViral\n"
        "échantillon exact\tread-2\t1\tSingle-end\t0.1\t0.1\t0.8\tNon-viral\n"
        "échantillon exact\tréad\t1\tPaired-incomplete\t1\t1\t0.8\tNon-viral\n"
    ).encode("utf-8")
    assert output.read_bytes() == expected
    assert not output.read_bytes().startswith(b"\xef\xbb\xbf")
    assert output.read_bytes().endswith(b"\n")


def test_normalized_compression_forms_decompress_to_identical_bytes(tmp_path):
    outputs = [
        tmp_path / "normalized.tsv",
        tmp_path / "normalized.tsv.gz",
        tmp_path / "normalized.tsv.zst",
    ]
    with _representative_store(tmp_path) as (store, _):
        assert [lyra._write_normalized(store, path) for path in outputs] == [3, 3, 3]

    decompressed = []
    for path in outputs:
        with util_file.open_or_gzopen(str(path), "rb") as stream:
            decompressed.append(stream.read())
    assert decompressed[0] == decompressed[1] == decompressed[2]


def test_normalized_empty_store_writes_header_and_preserves_threshold(tmp_path):
    output = tmp_path / "empty.tsv"
    with _artifact_store(
        tmp_path,
        [],
        [],
        threshold=Decimal("8e-1"),
        name="empty",
    ) as (store, _):
        assert store.threshold == Decimal("8e-1")
        assert lyra._write_normalized(store, output) == 0

    assert output.read_bytes() == ("\t".join(lyra.NORMALIZED_HEADER) + "\n").encode(
        "utf-8"
    )


def test_threshold_property_is_available_only_on_a_live_finalized_store(tmp_path):
    store = lyra.LyraFragmentStore("sample", work_dir=tmp_path)
    try:
        with pytest.raises(RuntimeError, match="not finalized"):
            _ = store.threshold
    finally:
        store.close()

    with pytest.raises(RuntimeError, match="closed"):
        _ = store.threshold


def test_normalized_no_hit_store_retains_every_nonviral_row(tmp_path):
    output = tmp_path / "no-hit.tsv"
    with _artifact_store(
        tmp_path,
        [("one", "0.1", "1"), ("two", "0.2", "1")],
        [_segment("one"), _segment("two")],
        name="no-hit",
    ) as (store, _):
        assert lyra._write_normalized(store, output) == 2

    rows = output.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3
    assert all(row.endswith("\tNon-viral") for row in rows[1:])


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (Decimal("0"), "0"),
        (Decimal("-0"), "0"),
        (Decimal("0.8000"), "0.8"),
        (Decimal("8e-1"), "0.8"),
        (Decimal("1e-7"), "0.0000001"),
        (Decimal("1.230000"), "1.23"),
    ],
)
def test_normalized_decimal_rendering_is_canonical(value, expected):
    assert lyra._canonical_output_decimal(value) == expected


@pytest.mark.parametrize(
    ("column", "invalid_value", "category"),
    [
        ("pairing", "single-end", "pairing"),
        ("call", "viral", "call"),
        ("threshold_text", "0.7", "threshold"),
    ],
)
def test_normalized_writer_rejects_unapproved_finalized_vocabulary(
    tmp_path,
    column,
    invalid_value,
    category,
):
    output = tmp_path / "invalid.tsv"
    with _artifact_store(
        tmp_path,
        [("read", "0.9", "1")],
        [_segment("read")],
        name="invalid-vocabulary-" + column,
    ) as (store, _):
        store._connection.execute(
            "UPDATE fragments SET {} = ?".format(column),
            (invalid_value,),
        )
        with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
            lyra._write_normalized(store, output)

    assert exc_info.value.category == category


@pytest.mark.parametrize(
    ("normalized", "summary", "viral_bam", "field"),
    [
        ("normalized.txt", "summary.tsv", "viral.bam", "normalized_output"),
        ("normalized.TSV", "summary.tsv", "viral.bam", "normalized_output"),
        ("normalized.tsv.GZ", "summary.tsv", "viral.bam", "normalized_output"),
        ("normalized.tsv", "summary.txt", "viral.bam", "summary_output"),
        ("normalized.tsv", "summary.tsv.gz", "viral.bam", "summary_output"),
        ("normalized.tsv", "summary.TSV", "viral.bam", "summary_output"),
        ("normalized.tsv", "summary.tsv", "viral.BAM", "viral_bam_output"),
        ("normalized.tsv", "summary.tsv", "viral.sam", "viral_bam_output"),
        ("normalized.tsv", "summary.tsv", "viral.bam.tmp", "viral_bam_output"),
    ],
)
def test_artifact_suffix_gate_is_case_sensitive_and_opens_nothing(
    tmp_path,
    normalized,
    summary,
    viral_bam,
    field,
):
    paths = [tmp_path / normalized, tmp_path / summary, tmp_path / viral_bam]
    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra._validate_artifact_output_suffixes(*paths)

    assert exc_info.value.category == "output_extension"
    assert exc_info.value.field == field
    assert not any(path.exists() for path in paths)


@pytest.mark.parametrize(
    "normalized_suffix",
    [".tsv", ".tsv.gz", ".tsv.zst"],
)
def test_artifact_suffix_gate_accepts_only_locked_destinations(
    tmp_path,
    normalized_suffix,
):
    paths = (
        tmp_path / ("normalized" + normalized_suffix),
        tmp_path / "summary.tsv",
        tmp_path / "viral.bam",
    )
    assert lyra._validate_artifact_output_suffixes(*paths) == tuple(
        str(path) for path in paths
    )
    assert not any(path.exists() for path in paths)
