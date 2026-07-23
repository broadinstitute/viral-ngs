from contextlib import contextmanager
from dataclasses import replace
from decimal import Decimal
import os
from types import SimpleNamespace

import pytest
import pysam

import viral_ngs.core.file as util_file
from viral_ngs.classify import lyra


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


def _consistent_counts():
    return lyra.LyraReconciliationCounts(
        input_bam_records=8,
        eligible_bam_records=4,
        score_records=4,
        fragments=3,
        single_end_fragments=1,
        complete_pair_fragments=1,
        incomplete_pair_fragments=1,
        viral_fragment_calls=1,
        nonviral_fragment_calls=2,
    )


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


def _logical_output_bytes(path):
    with util_file.open_or_gzopen(str(path), "rb") as stream:
        return stream.read()


def _summary_row(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == list(SUMMARY_HEADER)
    assert len(lines) == 2
    values = lines[1].split("\t")
    assert len(values) == len(SUMMARY_HEADER)
    return dict(zip(SUMMARY_HEADER, values))


def _assert_summary_equations(row, normalized_rows, actual_bam_records):
    values = {
        name: int(row[name])
        for name in SUMMARY_HEADER[2:]
    }
    assert normalized_rows == values["LYRA_FRAGMENTS"]
    assert (
        values["LYRA_ELIGIBLE_BAM_RECORDS"]
        == values["LYRA_SCORE_RECORDS"]
    )
    assert values["LYRA_SCORE_RECORDS"] == (
        values["LYRA_SINGLE_END_FRAGMENTS"]
        + values["LYRA_INCOMPLETE_PAIR_FRAGMENTS"]
        + 2 * values["LYRA_COMPLETE_PAIR_FRAGMENTS"]
    )
    assert values["LYRA_FRAGMENTS"] == (
        values["LYRA_SINGLE_END_FRAGMENTS"]
        + values["LYRA_COMPLETE_PAIR_FRAGMENTS"]
        + values["LYRA_INCOMPLETE_PAIR_FRAGMENTS"]
    )
    assert values["LYRA_FRAGMENTS"] == (
        values["LYRA_VIRAL_FRAGMENT_CALLS"]
        + values["LYRA_NONVIRAL_FRAGMENT_CALLS"]
    )
    assert 0 <= values["LYRA_OUTPUT_BAM_RECORDS"] <= values[
        "LYRA_INPUT_BAM_RECORDS"
    ]
    assert (values["LYRA_OUTPUT_BAM_RECORDS"] == 0) == (
        values["LYRA_VIRAL_FRAGMENT_CALLS"] == 0
    )
    assert values["LYRA_OUTPUT_BAM_RECORDS"] == actual_bam_records


def test_summary_has_exact_schema_and_one_canonical_data_row(tmp_path):
    output = tmp_path / "summary.tsv"
    store = SimpleNamespace(
        sample_id="sample exact",
        threshold=Decimal("0.8000"),
        counts=_consistent_counts(),
    )

    lyra._write_summary(store, output, output_bam_records=5)

    assert lyra.SUMMARY_HEADER == SUMMARY_HEADER
    assert output.read_bytes() == (
        "\t".join(SUMMARY_HEADER)
        + "\n"
        + "sample exact\t0.8\t8\t4\t4\t3\t1\t1\t1\t1\t2\t5\n"
    ).encode("utf-8")


@pytest.mark.parametrize(
    ("counts", "normalized_rows", "output_bam_records", "category"),
    [
        (_consistent_counts(), 2, 5, "normalized_row_count"),
        (
            replace(_consistent_counts(), eligible_bam_records=3),
            3,
            5,
            "eligible_score_record_count",
        ),
        (
            replace(
                _consistent_counts(),
                eligible_bam_records=5,
                score_records=5,
            ),
            3,
            5,
            "score_pairing_record_count",
        ),
        (
            replace(_consistent_counts(), fragments=4),
            4,
            5,
            "fragment_pairing_count",
        ),
        (
            replace(_consistent_counts(), nonviral_fragment_calls=1),
            3,
            5,
            "fragment_call_count",
        ),
        (_consistent_counts(), 3, -1, "output_bam_record_range"),
        (_consistent_counts(), 3, 9, "output_bam_record_range"),
        (_consistent_counts(), 3, 0, "output_bam_viral_equivalence"),
        (
            replace(
                _consistent_counts(),
                viral_fragment_calls=0,
                nonviral_fragment_calls=3,
            ),
            3,
            5,
            "output_bam_viral_equivalence",
        ),
    ],
)
def test_artifact_invariant_category_identifies_each_failed_equation(
    counts,
    normalized_rows,
    output_bam_records,
    category,
):
    with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
        lyra._validate_artifact_counts(
            counts,
            normalized_rows,
            output_bam_records,
        )

    assert exc_info.value.category == category


def test_artifact_invariants_accept_actual_bam_records_above_score_records():
    lyra._validate_artifact_counts(
        _consistent_counts(),
        normalized_rows=3,
        output_bam_records=5,
    )


def test_coordinator_calls_collaborators_in_summary_last_order(
    tmp_path,
    monkeypatch,
):
    calls = []
    counts = _consistent_counts()
    store = SimpleNamespace(counts=counts)
    paths = (
        tmp_path / "normalized.tsv",
        tmp_path / "summary.tsv",
        tmp_path / "viral.bam",
    )
    validated_paths = (
        str(tmp_path / "validated-normalized.tsv"),
        str(tmp_path / "validated-summary.tsv"),
        str(tmp_path / "validated-viral.bam"),
    )

    monkeypatch.setattr(
        lyra,
        "_validate_artifact_output_suffixes",
        lambda *args: calls.append(("suffixes", args)) or validated_paths,
    )
    monkeypatch.setattr(
        lyra,
        "_write_normalized",
        lambda current_store, path: calls.append(
            ("normalized", current_store, path)
        ) or 3,
    )
    monkeypatch.setattr(
        lyra,
        "_write_viral_bam",
        lambda current_store, path, work_dir=None: calls.append(
            ("bam", current_store, path, work_dir)
        ) or 5,
    )
    monkeypatch.setattr(
        lyra,
        "_assert_source_bam_identity",
        lambda current_store, stage: calls.append(
            ("source_identity", current_store, stage)
        ),
    )
    monkeypatch.setattr(
        lyra,
        "_validate_artifact_counts",
        lambda current_counts, normalized_rows, output_bam_records: calls.append(
            ("invariants", current_counts, normalized_rows, output_bam_records)
        ),
    )
    monkeypatch.setattr(
        lyra,
        "_write_summary",
        lambda current_store, path, output_bam_records: calls.append(
            ("summary", current_store, path, output_bam_records)
        ),
    )

    lyra.write_lyra_artifacts(
        store,
        *paths,
        work_dir=tmp_path / "work",
    )

    assert [call[0] for call in calls] == [
        "suffixes",
        "source_identity",
        "normalized",
        "bam",
        "invariants",
        "summary",
    ]
    assert calls[0][1] == paths
    assert calls[1][1:] == (store, "pre_generation")
    assert calls[2][1:] == (store, validated_paths[0])
    assert calls[3][1:] == (
        store,
        validated_paths[2],
        tmp_path / "work",
    )
    assert calls[4][1:] == (counts, 3, 5)
    assert calls[5][1:] == (store, validated_paths[1], 5)


def test_coordinator_converts_each_stateful_output_pathlike_once(
    tmp_path,
    monkeypatch,
):
    class StatefulPath(os.PathLike):
        def __init__(self, validated_path, bypass_path):
            self.validated_path = str(validated_path)
            self.bypass_path = str(bypass_path)
            self.fspath_calls = 0

        def __fspath__(self):
            self.fspath_calls += 1
            if self.fspath_calls == 1:
                return self.validated_path
            return self.bypass_path

    normalized = StatefulPath(
        tmp_path / "normalized.tsv",
        tmp_path / "normalized-bypass.txt",
    )
    summary = StatefulPath(
        tmp_path / "summary.tsv",
        tmp_path / "summary-bypass.txt",
    )
    viral_bam = StatefulPath(
        tmp_path / "viral.bam",
        tmp_path / "viral-bypass.sam",
    )
    outputs = (normalized, summary, viral_bam)
    validator_calls = []
    producer_paths = {}
    original_validator = lyra._validate_artifact_output_suffixes

    def count_validator(*paths):
        validator_calls.append(paths)
        return original_validator(*paths)

    def write_normalized(current_store, path):
        producer_paths["normalized"] = (path, os.fspath(path))
        return 3

    def write_viral_bam(current_store, path, work_dir=None):
        producer_paths["viral_bam"] = (path, os.fspath(path))
        return 5

    def write_summary(current_store, path, output_bam_records):
        producer_paths["summary"] = (path, os.fspath(path))

    monkeypatch.setattr(lyra, "_validate_artifact_output_suffixes", count_validator)
    monkeypatch.setattr(lyra, "_assert_source_bam_identity", lambda *args: None)
    monkeypatch.setattr(lyra, "_write_normalized", write_normalized)
    monkeypatch.setattr(lyra, "_write_viral_bam", write_viral_bam)
    monkeypatch.setattr(lyra, "_validate_artifact_counts", lambda *args: None)
    monkeypatch.setattr(lyra, "_write_summary", write_summary)

    lyra.write_lyra_artifacts(
        SimpleNamespace(counts=_consistent_counts()),
        *outputs,
        work_dir=tmp_path,
    )

    assert validator_calls == [outputs]
    for name, output in zip(("normalized", "summary", "viral_bam"), outputs):
        received_path, opened_path = producer_paths[name]
        assert received_path == output.validated_path
        assert opened_path == output.validated_path
        assert output.fspath_calls == 1
        assert not os.path.exists(output.bypass_path)


@pytest.mark.parametrize("producer", ["normalized", "bam"])
def test_coordinator_producer_failure_leaves_summary_unopened(
    tmp_path,
    monkeypatch,
    producer,
):
    summary = tmp_path / "summary.tsv"
    store = SimpleNamespace(counts=_consistent_counts())
    monkeypatch.setattr(lyra, "_assert_source_bam_identity", lambda *args: None)

    if producer == "normalized":
        monkeypatch.setattr(
            lyra,
            "_write_normalized",
            lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("normalized")),
        )
    else:
        monkeypatch.setattr(lyra, "_write_normalized", lambda *args, **kwargs: 3)
        monkeypatch.setattr(
            lyra,
            "_write_viral_bam",
            lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("bam")),
        )

    with pytest.raises(RuntimeError, match=producer):
        lyra.write_lyra_artifacts(
            store,
            tmp_path / "normalized.tsv",
            summary,
            tmp_path / "viral.bam",
            work_dir=tmp_path,
        )

    assert not summary.exists()


@pytest.mark.parametrize(
    ("normalized_rows", "output_bam_records", "category"),
    [
        (2, 5, "normalized_row_count"),
        (3, 0, "output_bam_viral_equivalence"),
    ],
)
def test_coordinator_count_mismatch_leaves_summary_unopened(
    tmp_path,
    monkeypatch,
    normalized_rows,
    output_bam_records,
    category,
):
    summary = tmp_path / "summary.tsv"
    store = SimpleNamespace(counts=_consistent_counts())
    monkeypatch.setattr(lyra, "_assert_source_bam_identity", lambda *args: None)
    monkeypatch.setattr(
        lyra,
        "_write_normalized",
        lambda *args, **kwargs: normalized_rows,
    )
    monkeypatch.setattr(
        lyra,
        "_write_viral_bam",
        lambda *args, **kwargs: output_bam_records,
    )

    with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
        lyra.write_lyra_artifacts(
            store,
            tmp_path / "normalized.tsv",
            summary,
            tmp_path / "viral.bam",
            work_dir=tmp_path,
        )

    assert exc_info.value.category == category
    assert not summary.exists()


def test_coordinator_invalid_suffix_opens_no_output(tmp_path):
    paths = (
        tmp_path / "normalized.tsv",
        tmp_path / "summary.txt",
        tmp_path / "viral.bam",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra.write_lyra_artifacts(
            SimpleNamespace(),
            *paths,
            work_dir=tmp_path,
        )

    assert exc_info.value.category == "output_extension"
    assert not any(path.exists() for path in paths)


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
        (Decimal("8e-1"), "0.8"),
        (Decimal("1e-7"), "0.0000001"),
        (Decimal("0e-100000"), "0"),
        (Decimal("-0e+100000"), "0"),
        (Decimal("1e-158"), "0." + "0" * 157 + "1"),
    ],
)
def test_threshold_validation_accepts_bounded_canonical_output(value, expected):
    validated = lyra.validate_lyra_threshold(value)

    assert validated is value
    assert lyra._canonical_output_decimal(validated) == expected
    assert len(expected) <= 160


@pytest.mark.parametrize("value", ["1e-159", "1e-100000"])
def test_threshold_validation_rejects_excessive_canonical_output_with_bounded_context(
    value,
):
    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra.validate_lyra_threshold(value)

    error = exc_info.value
    assert error.category == "threshold"
    assert error.field == "threshold"
    assert error.reason == "canonical threshold text exceeds 160 characters"
    assert error.offending_value == repr(value)
    assert len(str(error)) <= 320


def test_normalized_writer_canonicalizes_finalized_threshold_once(
    tmp_path,
    monkeypatch,
):
    threshold = Decimal("0.8000")
    fragments = [
        lyra.LyraFragment(
            read_id="one",
            n_scores=1,
            pairing=lyra.PAIRING_SINGLE_END,
            min_score=Decimal("0.1"),
            max_score=Decimal("0.1"),
            threshold=threshold,
            call=lyra.CALL_NONVIRAL,
        ),
        lyra.LyraFragment(
            read_id="two",
            n_scores=1,
            pairing=lyra.PAIRING_SINGLE_END,
            min_score=Decimal("0.9"),
            max_score=Decimal("0.9"),
            threshold=threshold,
            call=lyra.CALL_VIRAL,
        ),
    ]
    store = SimpleNamespace(
        sample_id="sample",
        threshold=threshold,
        iter_fragments=lambda: iter(fragments),
    )
    canonicalized = []
    original_canonicalizer = lyra._canonical_output_decimal

    def record_canonicalization(value):
        canonicalized.append(value)
        return original_canonicalizer(value)

    monkeypatch.setattr(lyra, "_canonical_output_decimal", record_canonicalization)

    assert lyra._write_normalized(store, tmp_path / "normalized.tsv") == 2

    assert canonicalized.count(threshold) == 1


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
                tmp_path / "mismatched-ids.bam",
                work_dir=tmp_path,
            )

    assert exc_info.value.category == "viral_read_id_count"
    assert exc_info.value.expected == 1
    assert exc_info.value.actual == 0


@pytest.mark.parametrize(
    ("case", "input_bam_records"),
    [("empty", 0), ("all_ineligible", 4)],
)
def test_coordinated_empty_and_all_ineligible_artifacts(
    tmp_path,
    case,
    input_bam_records,
):
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [{"ID": "source", "PN": "fixture"}],
    }
    records = []
    if case == "all_ineligible":
        records = [
            _segment("short", query_length=49),
            _segment("sequence-less", sequence_present=False),
            _segment("secondary", flag=0x4 | 0x100),
            _segment("supplementary", flag=0x4 | 0x800),
        ]
    source_bam = _write_bam(
        tmp_path,
        records,
        name=case + "-source.bam",
        header=header,
    )
    score_path = _write_score_table(tmp_path, [], name=case + "-scores.tsv")
    normalized = tmp_path / (case + "-normalized.tsv")
    summary = tmp_path / (case + "-summary.tsv")
    viral_bam = tmp_path / (case + "-viral.bam")

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_bam,
        "sample empty",
        Decimal("8e-1"),
        work_dir=tmp_path,
    ) as store:
        lyra.write_lyra_artifacts(
            store,
            normalized,
            summary,
            viral_bam,
            work_dir=tmp_path,
        )

    normalized_bytes = _logical_output_bytes(normalized)
    assert normalized_bytes == (
        "\t".join(lyra.NORMALIZED_HEADER) + "\n"
    ).encode("utf-8")
    expected_counts = [input_bam_records] + [0] * 9
    assert summary.read_bytes() == (
        "\t".join(SUMMARY_HEADER)
        + "\n"
        + "sample empty\t0.8\t"
        + "\t".join(str(value) for value in expected_counts)
        + "\n"
    ).encode("utf-8")

    with pysam.AlignmentFile(source_bam, "rb", check_sq=False) as source:
        source_header = source.header.to_dict()
    with pysam.AlignmentFile(viral_bam, "rb", check_sq=False) as output:
        assert output.header.to_dict() == source_header
        output_records = list(output.fetch(until_eof=True))
    assert output_records == []
    _assert_summary_equations(
        _summary_row(summary),
        normalized_rows=0,
        actual_bam_records=0,
    )


def test_coordinated_no_hit_compressed_normalized_retains_every_fragment(tmp_path):
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [{"ID": "source", "PN": "fixture"}],
    }
    records = [
        _segment("pair", flag=0x4 | 0x1 | 0x80),
        _segment("single"),
        _segment("incomplete", flag=0x4 | 0x1 | 0x40),
        _segment("pair", flag=0x4 | 0x1 | 0x40),
    ]
    source_bam = _write_bam(
        tmp_path,
        records,
        name="no-hit-source.bam",
        header=header,
    )
    score_path = _write_score_table(
        tmp_path,
        [
            ("pair", "0.9", "1"),
            ("single", "0.1", "1"),
            ("incomplete", "1.0", "0"),
            ("pair", "0.2", "0"),
        ],
        name="no-hit-scores.tsv",
    )
    normalized = tmp_path / "no-hit-normalized.tsv.zst"
    summary = tmp_path / "no-hit-summary.tsv"
    viral_bam = tmp_path / "no-hit-viral.bam"

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_bam,
        "sample no hit",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        lyra.write_lyra_artifacts(
            store,
            normalized,
            summary,
            viral_bam,
            work_dir=tmp_path,
        )

    expected_normalized = (
        "\t".join(lyra.NORMALIZED_HEADER)
        + "\n"
        + "sample no hit\tincomplete\t1\tPaired-incomplete\t1\t1\t0.8\tNon-viral\n"
        + "sample no hit\tpair\t2\tPaired-complete\t0.2\t0.9\t0.8\tNon-viral\n"
        + "sample no hit\tsingle\t1\tSingle-end\t0.1\t0.1\t0.8\tNon-viral\n"
    ).encode("utf-8")
    assert _logical_output_bytes(normalized) == expected_normalized

    row = _summary_row(summary)
    assert row == {
        "SAMPLE_ID": "sample no hit",
        "LYRA_THRESHOLD": "0.8",
        "LYRA_INPUT_BAM_RECORDS": "4",
        "LYRA_ELIGIBLE_BAM_RECORDS": "4",
        "LYRA_SCORE_RECORDS": "4",
        "LYRA_FRAGMENTS": "3",
        "LYRA_SINGLE_END_FRAGMENTS": "1",
        "LYRA_COMPLETE_PAIR_FRAGMENTS": "1",
        "LYRA_INCOMPLETE_PAIR_FRAGMENTS": "1",
        "LYRA_VIRAL_FRAGMENT_CALLS": "0",
        "LYRA_NONVIRAL_FRAGMENT_CALLS": "3",
        "LYRA_OUTPUT_BAM_RECORDS": "0",
    }
    with pysam.AlignmentFile(source_bam, "rb", check_sq=False) as source:
        source_header = source.header.to_dict()
    with pysam.AlignmentFile(viral_bam, "rb", check_sq=False) as output:
        assert output.header.to_dict() == source_header
        output_records = list(output.fetch(until_eof=True))
    assert output_records == []
    _assert_summary_equations(row, normalized_rows=3, actual_bam_records=0)


def test_coordinated_bam_fidelity_counts_every_same_qname_companion(tmp_path):
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
        name="coordinated-fidelity-source.bam",
        header=header,
    )
    score_path = _write_score_table(
        tmp_path,
        [
            ("viral", "0.9", "0"),
            ("nonviral", "0.1", "1"),
            ("other", "0.1", "1"),
        ],
        name="coordinated-fidelity.tsv",
    )
    normalized = tmp_path / "coordinated-fidelity-normalized.tsv"
    summary = tmp_path / "coordinated-fidelity-summary.tsv"
    viral_bam = tmp_path / "coordinated-fidelity-viral.bam"

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_bam,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        lyra.write_lyra_artifacts(
            store,
            normalized,
            summary,
            viral_bam,
            work_dir=tmp_path,
        )

    normalized_rows = _logical_output_bytes(normalized).decode("utf-8").splitlines()
    assert normalized_rows[0].split("\t") == list(lyra.NORMALIZED_HEADER)
    assert [row.split("\t")[1] for row in normalized_rows[1:]] == [
        "nonviral",
        "other",
        "viral",
    ]
    assert len(normalized_rows[1:]) == 3 < len(records)

    with pysam.AlignmentFile(source_bam, "rb", check_sq=False) as source:
        source_header = source.header.to_dict()
        expected_records = [
            _record_snapshot(record)
            for record in source.fetch(until_eof=True)
            if record.query_name == "viral"
        ]
    with pysam.AlignmentFile(viral_bam, "rb", check_sq=False) as output:
        output_header = output.header.to_dict()
        output_records = [
            _record_snapshot(record)
            for record in output.fetch(until_eof=True)
        ]

    assert output_header == source_header
    assert output_header["PG"] == header["PG"]
    assert output_records == expected_records
    assert [record[0] for record in output_records] == ["viral"] * 5
    row = _summary_row(summary)
    assert row["LYRA_SCORE_RECORDS"] == "3"
    assert row["LYRA_FRAGMENTS"] == "3"
    assert row["LYRA_VIRAL_FRAGMENT_CALLS"] == "1"
    assert row["LYRA_OUTPUT_BAM_RECORDS"] == "5"
    assert int(row["LYRA_OUTPUT_BAM_RECORDS"]) > int(row["LYRA_SCORE_RECORDS"])
    _assert_summary_equations(
        row,
        normalized_rows=len(normalized_rows) - 1,
        actual_bam_records=len(output_records),
    )
