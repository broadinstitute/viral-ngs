from dataclasses import FrozenInstanceError
from dataclasses import replace
from decimal import Decimal
from decimal import localcontext
from pathlib import Path

import pytest
import pysam

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


def test_iter_fragments_is_repeatable_ordered_and_immutable(tmp_path):
    read_ids = ["read-2", "read-10", "Read", "read/1"]
    score_path = _write_score_table(
        tmp_path,
        [(read_id, "0.9", "0") for read_id in reversed(read_ids)],
        name="repeatable.tsv",
    )
    bam_path = _write_bam(
        tmp_path,
        [_segment(read_id) for read_id in read_ids],
        name="repeatable.bam",
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        first_pass = list(store.iter_fragments())
        second_pass = list(store.iter_fragments())

        assert first_pass == second_pass
        assert [fragment.read_id for fragment in first_pass] == sorted(
            read_ids,
            key=lambda read_id: read_id.encode("utf-8"),
        )
        assert all(isinstance(fragment, lyra.LyraFragment) for fragment in first_pass)
        with pytest.raises(FrozenInstanceError):
            first_pass[0].call = "Non-viral"

    with pytest.raises(RuntimeError, match="closed"):
        list(store.iter_fragments())


def test_reconciliation_preserves_exact_interior_whitespace_qname(tmp_path):
    read_id = "read  with interior whitespace"
    score_path = _write_score_table(
        tmp_path,
        [(read_id, "0.9", "1")],
        name="interior-whitespace.tsv",
    )
    bam_path = _write_bam(
        tmp_path,
        [_segment(read_id)],
        name="interior-whitespace.bam",
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        fragments = list(store.iter_fragments())

    assert fragments == [
        lyra.LyraFragment(
            read_id=read_id,
            n_scores=1,
            pairing="Single-end",
            min_score=Decimal("0.9"),
            max_score=Decimal("0.9"),
            threshold=Decimal("0.8"),
            call="Viral",
        )
    ]


def test_iter_viral_read_ids_preserves_exact_ids_in_utf8_blob_order(tmp_path):
    read_ids = [
        "read-2",
        "read-10",
        "prefix",
        "prefix!",
        "Read",
        "read",
        "read/1",
        "read/2",
        "001",
        "1",
        "e\u0301",
        "é",
    ]
    nonviral_ids = {"prefix", "read/2"}
    score_path = _write_score_table(
        tmp_path,
        [
            (read_id, "0.1" if read_id in nonviral_ids else "0.9", "1")
            for read_id in reversed(read_ids)
        ],
        name="exact-ids.tsv",
    )
    bam_path = _write_bam(
        tmp_path,
        [_segment(read_id) for read_id in reversed(read_ids)],
        name="exact-ids.bam",
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        fragments = list(store.iter_fragments())
        first_viral_pass = list(store.iter_viral_read_ids())
        second_viral_pass = list(store.iter_viral_read_ids())

        expected_order = sorted(read_ids, key=lambda read_id: read_id.encode("utf-8"))
        assert [fragment.read_id for fragment in fragments] == expected_order
        assert first_viral_pass == second_viral_pass
        assert first_viral_pass == [
            read_id for read_id in expected_order if read_id not in nonviral_ids
        ]

    with pytest.raises(RuntimeError, match="closed"):
        list(store.iter_viral_read_ids())


def test_independently_shuffled_inputs_produce_identical_fragments_and_counts(
    tmp_path,
):
    score_rows = [
        ("single", "0.8", "0"),
        ("pair", "1", "0"),
        ("pair", "0.8", "1"),
        ("incomplete", "1", "1"),
        ("nonviral", "0.1", "1"),
    ]
    bam_specs = [
        ("pair", 0x4 | 0x1 | 0x40),
        ("single", 0x4),
        ("incomplete", 0x4 | 0x1 | 0x80),
        ("pair", 0x4 | 0x1 | 0x80),
        ("nonviral", 0x4),
    ]
    results = []

    for index, (current_scores, current_bam_specs) in enumerate(
        (
            (score_rows, bam_specs),
            (
                [score_rows[i] for i in (3, 1, 4, 0, 2)],
                [bam_specs[i] for i in (4, 2, 0, 3, 1)],
            ),
        )
    ):
        score_path = _write_score_table(
            tmp_path,
            current_scores,
            name="shuffle-{}.tsv".format(index),
        )
        bam_path = _write_bam(
            tmp_path,
            [_segment(read_id, flag=flag) for read_id, flag in current_bam_specs],
            name="shuffle-{}.bam".format(index),
        )
        with lyra.reconcile_lyra_fragments(
            score_path,
            bam_path,
            "sample",
            "0.8",
            work_dir=tmp_path,
        ) as store:
            results.append(
                (
                    list(store.iter_fragments()),
                    store.counts,
                    list(store.iter_viral_read_ids()),
                )
            )

    assert results[0] == results[1]


def test_empty_and_all_ineligible_inputs_have_complete_zero_fragment_counts(tmp_path):
    score_path = _write_score_table(tmp_path, [], name="empty.tsv")
    zero_bam_path = _write_bam(tmp_path, [], name="zero-record.bam")

    with lyra.reconcile_lyra_fragments(
        score_path,
        zero_bam_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        assert list(store.iter_fragments()) == []
        assert list(store.iter_viral_read_ids()) == []
        assert store.counts == lyra.LyraReconciliationCounts(0, 0, 0, 0, 0, 0, 0, 0, 0)

    ineligible_records = [
        _segment("short", query_length=49),
        _segment("sequence-less", sequence_present=False),
        _segment("secondary", flag=0x4 | 0x100),
        _segment("supplementary", flag=0x4 | 0x800),
    ]
    ineligible_bam_path = _write_bam(
        tmp_path,
        ineligible_records,
        name="all-ineligible.bam",
    )
    with lyra.reconcile_lyra_fragments(
        score_path,
        ineligible_bam_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        assert list(store.iter_fragments()) == []
        assert list(store.iter_viral_read_ids()) == []
        assert store.counts == lyra.LyraReconciliationCounts(4, 0, 0, 0, 0, 0, 0, 0, 0)


def test_large_reconciliation_stays_file_backed_and_cleans_up(tmp_path):
    fragment_count = lyra.SQLITE_COMMIT_INTERVAL + 1
    score_path = _write_score_table(
        tmp_path,
        (
            ("large-{:05d}".format(index), "0.9", "0")
            for index in reversed(range(fragment_count))
        ),
        name="large.tsv",
    )
    bam_path = _write_bam(
        tmp_path,
        (_segment("large-{:05d}".format(index)) for index in range(fragment_count)),
        name="large.bam",
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        database_path = Path(store.database_path)
        fragment_iterator = store.iter_fragments()
        viral_id_iterator = store.iter_viral_read_ids()

        assert database_path.is_file()
        assert database_path.stat().st_size > 0
        assert not isinstance(fragment_iterator, (list, set, dict))
        assert not isinstance(viral_id_iterator, (list, set, dict))
        assert sum(1 for _ in fragment_iterator) == fragment_count
        assert sum(1 for _ in viral_id_iterator) == fragment_count
        assert store.counts.fragments == fragment_count
        assert store.counts.viral_fragment_calls == fragment_count
        for name in dir(store):
            if name.startswith("_"):
                continue
            assert not isinstance(getattr(store, name), (list, set, dict))

    assert not database_path.exists()


def test_ignored_companions_change_only_physical_bam_record_count(tmp_path):
    score_rows = [("incomplete", "1", "0")]
    baseline_records = [_segment("incomplete", flag=0x4 | 0x1 | 0x40)]
    companion_records = baseline_records + [
        _segment("incomplete", query_length=49, flag=0x4 | 0x1 | 0x80),
        _segment("incomplete", sequence_present=False, flag=0x4 | 0x1 | 0x80),
        _segment("incomplete", flag=0x4 | 0x100 | 0x1 | 0x80),
        _segment("incomplete", flag=0x4 | 0x800 | 0x1 | 0x80),
    ]
    results = []

    for index, records in enumerate((baseline_records, companion_records)):
        score_path = _write_score_table(
            tmp_path,
            score_rows,
            name="companions-{}.tsv".format(index),
        )
        bam_path = _write_bam(
            tmp_path,
            records,
            name="companions-{}.bam".format(index),
        )
        with lyra.reconcile_lyra_fragments(
            score_path,
            bam_path,
            "sample",
            "0",
            work_dir=tmp_path,
        ) as store:
            results.append((list(store.iter_fragments()), store.counts))

    baseline_fragments, baseline_counts = results[0]
    companion_fragments, companion_counts = results[1]
    assert baseline_fragments == companion_fragments
    assert companion_counts.input_bam_records == 5
    assert replace(
        companion_counts,
        input_bam_records=baseline_counts.input_bam_records,
    ) == baseline_counts
