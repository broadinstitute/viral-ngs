from decimal import Decimal
from decimal import localcontext

import pytest
import pysam

from viral_ngs.classify import lyra


def _segment(query_name, *, flag=0x4, query_length=50):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = flag
    record.query_sequence = "A" * query_length
    record.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    record.reference_id = -1
    record.reference_start = -1
    return record


def _write_bam(tmp_path, records, name="reads.bam"):
    path = tmp_path / name
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6"}}) as bam:
        for record in records:
            bam.write(record)
    return path


def _write_score_table(tmp_path, rows, name="scores.tsv"):
    path = tmp_path / name
    path.write_text(
        "read_id\tscore\tcall\n"
        + "".join("{}\t{}\t{}\n".format(*row) for row in rows)
    )
    return path


def _stored_fragments(store):
    return store._connection.execute(
        "SELECT * FROM fragments ORDER BY read_id_key ASC"
    ).fetchall()


def _reconcile_rows(
    tmp_path,
    score_rows,
    bam_records,
    *,
    threshold="0.8",
    name="run",
):
    score_path = _write_score_table(tmp_path, score_rows, name=name + ".tsv")
    bam_path = _write_bam(tmp_path, bam_records, name=name + ".bam")
    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        "sample exact",
        threshold,
        work_dir=tmp_path,
    ) as store:
        rows = [dict(row) for row in _stored_fragments(store)]
        counts = store.counts
    return rows, counts


@pytest.mark.parametrize(
    ("score", "threshold", "expected_call"),
    [
        ("0", "0", "Viral"),
        ("0.799999", "0.8", "Non-viral"),
        ("0.8", "0.8", "Viral"),
        ("1", "1", "Viral"),
    ],
)
def test_fragment_calls_single_end_use_inclusive_threshold(
    tmp_path,
    score,
    threshold,
    expected_call,
):
    rows, _ = _reconcile_rows(
        tmp_path,
        [("single", score, "0")],
        [_segment("single")],
        threshold=threshold,
        name="single-{}-{}".format(score.replace(".", "_"), threshold.replace(".", "_")),
    )

    assert len(rows) == 1
    assert rows[0]["call"] == expected_call
    assert rows[0]["pairing"] == "Single-end"


@pytest.mark.parametrize(
    ("scores", "expected_call"),
    [
        (("0.8", "1"), "Viral"),
        (("0.79", "1"), "Non-viral"),
        (("1", "0.79"), "Non-viral"),
        (("0", "0.79"), "Non-viral"),
    ],
)
def test_fragment_calls_complete_pair_require_both_scores_independent_of_order(
    tmp_path,
    scores,
    expected_call,
):
    calls = []
    for index, current_scores in enumerate((scores, tuple(reversed(scores)))):
        rows, _ = _reconcile_rows(
            tmp_path,
            [("pair", score, "0") for score in current_scores],
            [
                _segment("pair", flag=0x4 | 0x1 | 0x80),
                _segment("pair", flag=0x4 | 0x1 | 0x40),
            ],
            name="pair-{}".format(index),
        )
        calls.append(rows[0]["call"])
        assert rows[0]["pairing"] == "Paired-complete"

    assert calls == [expected_call, expected_call]


def test_fragment_calls_incomplete_pair_are_nonviral_at_zero_threshold(tmp_path):
    rows, counts = _reconcile_rows(
        tmp_path,
        [("incomplete", "1", "1")],
        [_segment("incomplete", flag=0x4 | 0x1 | 0x40)],
        threshold="0",
        name="incomplete",
    )

    assert rows[0]["pairing"] == "Paired-incomplete"
    assert rows[0]["call"] == "Non-viral"
    assert counts.incomplete_pair_fragments == 1
    assert counts.nonviral_fragment_calls == 1


def test_native_call_never_changes_recomputed_fragment_call(tmp_path):
    results = []
    for native_call in ("0", "1"):
        rows, counts = _reconcile_rows(
            tmp_path,
            [("native-disagreement", "0.9", native_call)],
            [_segment("native-disagreement")],
            name="native-call-{}".format(native_call),
        )
        results.append((rows, counts))

    assert results[0] == results[1]
    assert results[0][0][0]["call"] == "Viral"


def test_finalized_fragment_rows_and_counts_have_complete_locked_contract(tmp_path):
    rows, counts = _reconcile_rows(
        tmp_path,
        [
            ("single", "0.75", "1"),
            ("pair", "0.90", "0"),
            ("pair", "0.85", "0"),
            ("incomplete", "1.0", "1"),
        ],
        [
            _segment("pair", flag=0x4 | 0x1 | 0x80),
            _segment("single"),
            _segment("incomplete", flag=0x4 | 0x1 | 0x40),
            _segment("pair", flag=0x4 | 0x1 | 0x40),
        ],
        threshold="0.8",
        name="complete-contract",
    )
    by_read_id = {
        row["read_id_key"].decode("utf-8"): row
        for row in rows
    }

    assert by_read_id["single"] == {
        "read_id_key": b"single",
        "n_scores": 1,
        "pairing": "Single-end",
        "min_score_text": "0.75",
        "max_score_text": "0.75",
        "threshold_text": "0.8",
        "call": "Non-viral",
    }
    assert by_read_id["pair"] == {
        "read_id_key": b"pair",
        "n_scores": 2,
        "pairing": "Paired-complete",
        "min_score_text": "0.85",
        "max_score_text": "0.9",
        "threshold_text": "0.8",
        "call": "Viral",
    }
    assert by_read_id["incomplete"] == {
        "read_id_key": b"incomplete",
        "n_scores": 1,
        "pairing": "Paired-incomplete",
        "min_score_text": "1",
        "max_score_text": "1",
        "threshold_text": "0.8",
        "call": "Non-viral",
    }
    assert counts == lyra.LyraReconciliationCounts(
        input_bam_records=4,
        eligible_bam_records=4,
        score_records=4,
        fragments=3,
        single_end_fragments=1,
        complete_pair_fragments=1,
        incomplete_pair_fragments=1,
        viral_fragment_calls=1,
        nonviral_fragment_calls=2,
    )


def test_equal_score_spellings_finalize_canonically_without_context_rounding(tmp_path):
    exact_text = "0." + ("1234567890" * 8) + "1"
    with localcontext() as context:
        context.prec = 3
        assert lyra._canonical_decimal_text(exact_text + "0000") == exact_text

        results = []
        for index, score_rows in enumerate(
            (
                [("pair", exact_text, "0"), ("pair", exact_text + "0000", "1")],
                [("pair", exact_text + "0000", "1"), ("pair", exact_text, "0")],
            )
        ):
            rows, _ = _reconcile_rows(
                tmp_path,
                score_rows,
                [
                    _segment("pair", flag=0x4 | 0x1 | 0x40),
                    _segment("pair", flag=0x4 | 0x1 | 0x80),
                ],
                threshold="0",
                name="canonical-{}".format(index),
            )
            results.append(rows[0])

    assert results[0] == results[1]
    assert results[0]["min_score_text"] == exact_text
    assert results[0]["max_score_text"] == exact_text
    assert Decimal(results[0]["min_score_text"]).as_tuple() == Decimal(exact_text).as_tuple()


def test_late_malformed_score_fails_before_bam_traversal_or_finalization(
    tmp_path,
    monkeypatch,
):
    score_path = tmp_path / "late-malformed.tsv"
    score_path.write_text(
        "read_id\tscore\tcall\n"
        "valid\t0.8\t1\n"
        "malformed\tnot-a-score\t0\n"
    )
    bam_path = tmp_path / "must-not-open.bam"
    calls = []

    def fail_bam_traversal(*args, **kwargs):
        calls.append("bam")
        raise AssertionError("BAM traversal occurred before score validation completed")

    def fail_finalization(*args, **kwargs):
        calls.append("finalize")
        raise AssertionError("finalization occurred after invalid score input")

    monkeypatch.setattr(lyra.pysam, "AlignmentFile", fail_bam_traversal)
    monkeypatch.setattr(lyra, "_finalize_fragments", fail_finalization)

    with pytest.raises(lyra.LyraInputError) as exc_info:
        with lyra.reconcile_lyra_fragments(
            score_path,
            bam_path,
            "sample",
            "0.8",
            work_dir=tmp_path,
        ):
            pytest.fail("an invalid reconciliation must not yield a store")

    assert exc_info.value.category == "score_syntax"
    assert calls == []


def test_reconciliation_error_prevents_finalization_and_store_exposure(
    tmp_path,
    monkeypatch,
):
    score_path = _write_score_table(
        tmp_path,
        [("score-only", "0.9", "1")],
        name="reconciliation-error.tsv",
    )
    bam_path = _write_bam(tmp_path, [], name="reconciliation-error.bam")
    finalized = []
    original_finalize = lyra._finalize_fragments

    def track_finalization(*args, **kwargs):
        finalized.append(True)
        return original_finalize(*args, **kwargs)

    monkeypatch.setattr(lyra, "_finalize_fragments", track_finalization)

    with pytest.raises(lyra.LyraReconciliationError) as exc_info:
        with lyra.reconcile_lyra_fragments(
            score_path,
            bam_path,
            "sample",
            "0.8",
            work_dir=tmp_path,
        ):
            pytest.fail("an invalid reconciliation must not yield a store")

    assert exc_info.value.category == "score_only"
    assert finalized == []


def test_store_sample_and_exact_fragment_ids_form_unique_output_keys(tmp_path):
    sample_id = "échantillon exact"
    read_ids = ["read-1", "Read-1", "read-1/1", "read-1/2"]
    score_path = _write_score_table(
        tmp_path,
        [(read_id, "0.9", "0") for read_id in reversed(read_ids)],
        name="unique-keys.tsv",
    )
    bam_path = _write_bam(
        tmp_path,
        [_segment(read_id) for read_id in read_ids],
        name="unique-keys.bam",
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        sample_id,
        "0.8",
        work_dir=tmp_path,
    ) as store:
        output_keys = [
            (store.sample_id, row["read_id_key"].decode("utf-8"))
            for row in _stored_fragments(store)
        ]

        assert store.sample_id is sample_id
        assert len(output_keys) == len(set(output_keys)) == len(read_ids)
