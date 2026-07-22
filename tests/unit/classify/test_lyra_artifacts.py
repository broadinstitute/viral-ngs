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


def _write_bam(tmp_path, records, name="reads.bam", header=None):
    path = tmp_path / name
    with pysam.AlignmentFile(
        str(path),
        "wb",
        header=header or {"HD": {"VN": "1.6"}},
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


def _rich_segment(
    query_name,
    *,
    flag=0,
    sequence="A" * 60,
    reference_start=10,
    tag_value=1,
):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = flag
    record.query_sequence = sequence
    if sequence is not None:
        record.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    if flag & 0x4:
        record.reference_id = -1
        record.reference_start = -1
    else:
        record.reference_id = 0
        record.reference_start = reference_start
        record.mapping_quality = 37
        record.cigarstring = "{}M".format(len(sequence))
    record.next_reference_id = -1
    record.next_reference_start = -1
    record.template_length = 0
    record.set_tags(
        [("RG", "rg1"), ("NM", tag_value), ("XY", float(tag_value))]
    )
    return record


def _record_snapshot(record):
    return (
        record.query_name,
        record.flag,
        record.reference_id,
        record.reference_start,
        record.mapping_quality,
        record.cigarstring,
        record.next_reference_id,
        record.next_reference_start,
        record.template_length,
        record.query_sequence,
        None if record.query_qualities is None else tuple(record.query_qualities),
        record.get_tags(with_value_type=True),
    )


def test_viral_bam_fidelity_preserves_all_exact_qname_records(tmp_path):
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "ref", "LN": 1000}],
        "RG": [{"ID": "rg1", "SM": "sample"}],
        "PG": [{"ID": "existing", "PN": "source", "VN": "1.0"}],
    }
    records = [
        _rich_segment("viral", tag_value=1),
        _rich_segment("nonviral", reference_start=20, tag_value=2),
        _rich_segment("viral", sequence="C" * 49, reference_start=30, tag_value=3),
        _rich_segment("other", reference_start=40, tag_value=4),
        _rich_segment("viral", flag=0x4, sequence=None, tag_value=5),
        _rich_segment("viral", flag=0x100, reference_start=50, tag_value=6),
        _rich_segment("viral", flag=0x800, reference_start=60, tag_value=7),
    ]
    source_bam = _write_bam(
        tmp_path,
        records,
        name="fidelity-source.bam",
        header=header,
    )
    score_path = _write_score_table(
        tmp_path,
        [("viral", "0.9", "0"), ("nonviral", "0.1", "1"), ("other", "0.1", "1")],
        name="fidelity.tsv",
    )
    output_bam = tmp_path / "viral-output.bam"

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_bam,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        assert store.counts.score_records == 3
        emitted_count = lyra._write_viral_bam(
            store,
            source_bam,
            output_bam,
            work_dir=tmp_path,
        )

    with pysam.AlignmentFile(source_bam, "rb", check_sq=False) as source:
        source_header = source.header.to_dict()
        expected_records = [
            _record_snapshot(record)
            for record in source.fetch(until_eof=True)
            if record.query_name == "viral"
        ]
    with pysam.AlignmentFile(output_bam, "rb", check_sq=False) as output:
        output_header = output.header.to_dict()
        output_records = [
            _record_snapshot(record)
            for record in output.fetch(until_eof=True)
        ]

    assert emitted_count == len(output_records) == 5
    assert emitted_count > 3
    assert output_header == source_header
    assert output_header["PG"] == header["PG"]
    assert output_records == expected_records
    assert not (tmp_path / "viral-output.bam.bai").exists()
    assert not (tmp_path / "viral-output.bai").exists()


def test_viral_bam_empty_call_set_preserves_header_and_returns_zero(tmp_path):
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "ref", "LN": 1000}],
        "PG": [{"ID": "existing", "PN": "source"}],
    }
    source_bam = _write_bam(
        tmp_path,
        [_rich_segment("nonviral")],
        name="empty-viral-source.bam",
        header=header,
    )
    score_path = _write_score_table(
        tmp_path,
        [("nonviral", "0.1", "1")],
        name="empty-viral.tsv",
    )
    output_bam = tmp_path / "empty-viral-output.bam"

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_bam,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        assert lyra._write_viral_bam(
            store,
            source_bam,
            output_bam,
            work_dir=tmp_path,
        ) == 0

    with pysam.AlignmentFile(source_bam, "rb", check_sq=False) as source:
        source_header = source.header.to_dict()
    with pysam.AlignmentFile(output_bam, "rb", check_sq=False) as output:
        assert output.header.to_dict() == source_header
        assert list(output.fetch(until_eof=True)) == []


def test_viral_bam_rejects_a_mismatched_exact_id_cursor(tmp_path, monkeypatch):
    with _artifact_store(
        tmp_path,
        [("viral", "0.9", "1")],
        [_segment("viral")],
        name="viral-id-mismatch",
    ) as (store, source_bam):
        monkeypatch.setattr(store, "iter_viral_read_ids", lambda: iter(()))
        with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
            lyra._write_viral_bam(
                store,
                source_bam,
                tmp_path / "mismatched-ids.bam",
                work_dir=tmp_path,
            )

    assert exc_info.value.category == "viral_read_id_count"
    assert exc_info.value.expected == 1
    assert exc_info.value.actual == 0
