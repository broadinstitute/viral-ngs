import inspect
import itertools

import pytest
import pysam

from viral_ngs.classify import lyra


_ROLE_LABELS = ("U--", "U1-", "U-2", "U12", "P--", "P1-", "P-2", "P12")


def _role_counts(*role_indexes):
    return tuple(role_indexes.count(index) for index in range(len(_ROLE_LABELS)))


def _insert_evidence(
    store,
    read_id,
    *,
    score_count=0,
    score_line_numbers=(),
    bam_role_counts=None,
):
    bam_role_counts = bam_role_counts or _role_counts()
    eligible_bam_count = sum(bam_role_counts)
    retained_lines = tuple(score_line_numbers[:2]) + (None, None)
    store._connection.execute(
        """
        INSERT INTO evidence (
            read_id_key,
            score_count,
            score_1_text,
            score_1_line,
            score_2_text,
            score_2_line,
            eligible_bam_count,
            unpaired_none_count,
            unpaired_r1_count,
            unpaired_r2_count,
            unpaired_both_count,
            paired_none_count,
            paired_r1_count,
            paired_r2_count,
            paired_both_count
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            lyra._read_id_key(read_id),
            score_count,
            "0.5" if score_count else None,
            retained_lines[0],
            "0.5" if score_count > 1 else None,
            retained_lines[1],
            eligible_bam_count,
            *bam_role_counts,
        ),
    )


def _reconciliation_error(store, score_path="scores.tsv", bam_path="reads.bam"):
    with pytest.raises(lyra.LyraReconciliationError) as exc_info:
        lyra._validate_reconciliation(store, score_path, bam_path)
    return exc_info.value


def _segment(query_name, flag=0x4):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = flag
    record.query_sequence = "A" * 50
    record.query_qualities = pysam.qualitystring_to_array("I" * 50)
    record.reference_id = -1
    record.reference_start = -1
    return record


def _write_score_table(path, read_ids):
    path.write_text(
        "read_id\tscore\tcall\n"
        + "".join("{}\t0.5\t0\n".format(read_id) for read_id in read_ids)
    )


def _write_bam(path, read_ids):
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6"}}) as bam:
        for read_id in read_ids:
            bam.write(_segment(read_id))


def test_pairing_state_exhaustively_classifies_one_and_two_record_role_shapes():
    valid_single_record_shapes = {
        0: lyra.PAIRING_SINGLE_END,
        5: lyra.PAIRING_INCOMPLETE,
        6: lyra.PAIRING_INCOMPLETE,
    }
    for role_index in range(len(_ROLE_LABELS)):
        expected_pairing = valid_single_record_shapes.get(role_index)
        expected_category = None if expected_pairing is not None else "bam_pairing_flags"

        assert lyra._pairing_state(1, _role_counts(role_index)) == (
            expected_pairing,
            expected_category,
        )

    for role_indexes in itertools.combinations_with_replacement(
        range(len(_ROLE_LABELS)),
        2,
    ):
        expected_pairing = (
            lyra.PAIRING_COMPLETE if role_indexes == (5, 6) else None
        )
        expected_category = (
            None if expected_pairing is not None else "bam_pairing_flags"
        )

        assert lyra._pairing_state(2, _role_counts(*role_indexes)) == (
            expected_pairing,
            expected_category,
        )


@pytest.mark.parametrize(
    "bam_role_counts",
    [
        _role_counts(5, 6, 0),
        _role_counts(5, 6, 5),
        _role_counts(0, 0, 0, 0),
        _role_counts(7, 7, 7, 7, 7),
    ],
)
def test_pairing_state_prioritizes_record_limit_for_more_than_two_records(
    bam_role_counts,
):
    assert lyra._pairing_state(sum(bam_role_counts), bam_role_counts) == (
        None,
        "bam_record_limit",
    )


def test_pairing_state_uses_only_count_and_eight_fixed_role_counters():
    assert tuple(inspect.signature(lyra._pairing_state).parameters) == (
        "eligible_bam_count",
        "bam_role_counts",
    )
    role_counts = _role_counts(5, 6)

    assert lyra._pairing_state(2, role_counts) == (lyra.PAIRING_COMPLETE, None)


def test_reconciliation_error_exposes_bounded_structured_intrinsic_context():
    error = lyra.LyraReconciliationError(
        category="bam_record_limit",
        score_path="scores.tsv",
        bam_path="reads.bam",
        read_id="read-1",
        score_count=9,
        eligible_bam_count=8,
        score_line_numbers=(2, 3, 4),
        bam_role_counts=(1, 2, 3, 4, 5, 6, 7, 8, 9),
        reason="too many eligible BAM records",
    )

    assert isinstance(error, lyra.LyraInputError)
    assert error.category == "bam_record_limit"
    assert error.score_path == "scores.tsv"
    assert error.bam_path == "reads.bam"
    assert error.read_id == "read-1"
    assert error.score_count == 9
    assert error.eligible_bam_count == 8
    assert error.score_line_numbers == (2, 3)
    assert error.bam_role_counts == tuple(
        zip(_ROLE_LABELS, range(1, len(_ROLE_LABELS) + 1))
    )
    assert error.reason == "too many eligible BAM records"
    assert len(error.score_line_numbers) <= 2
    assert len(error.bam_role_counts) <= 8


def test_reconciliation_error_message_uses_bounded_rendering_for_text_context():
    control_text = "unsafe\r\n" + ("x" * 10000)
    error = lyra.LyraReconciliationError(
        category="bam_pairing_flags",
        score_path=control_text,
        bam_path=control_text,
        read_id=control_text,
        score_count=2,
        eligible_bam_count=2,
        score_line_numbers=(2, 3),
        bam_role_counts=_role_counts(5, 5),
        reason="incoherent pairing flags",
    )

    assert "\r" not in str(error)
    assert "\n" not in str(error)
    assert "\\r" in str(error)
    assert "\\n" in str(error)
    assert len(str(error)) < 1000


@pytest.mark.parametrize(
    ("score_count", "bam_role_counts", "expected_category"),
    [
        (1, _role_counts(), "score_only"),
        (0, _role_counts(0), "bam_only"),
        (2, _role_counts(0), "record_count_mismatch"),
    ],
)
def test_cross_file_multiplicity_errors_have_stable_categories(
    tmp_path,
    score_count,
    bam_role_counts,
    expected_category,
):
    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        _insert_evidence(
            store,
            "read-1",
            score_count=score_count,
            score_line_numbers=(2, 3),
            bam_role_counts=bam_role_counts,
        )

        error = _reconciliation_error(store)

    assert error.category == expected_category
    assert error.score_count == score_count
    assert error.eligible_bam_count == sum(bam_role_counts)


def test_equal_positive_multiplicities_and_empty_evidence_pass_validation(tmp_path):
    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        lyra._validate_reconciliation(store, "scores.tsv", "reads.bam")
        _insert_evidence(
            store,
            "single",
            score_count=1,
            score_line_numbers=(2,),
            bam_role_counts=_role_counts(0),
        )
        _insert_evidence(
            store,
            "pair",
            score_count=2,
            score_line_numbers=(3, 4),
            bam_role_counts=_role_counts(5, 6),
        )

        assert lyra._validate_reconciliation(
            store,
            "scores.tsv",
            "reads.bam",
        ) is None


@pytest.mark.parametrize(
    ("score_count", "bam_role_counts", "expected_category"),
    [
        (0, _role_counts(5, 6, 0), "bam_record_limit"),
        (0, _role_counts(5, 5), "bam_pairing_flags"),
        (1, _role_counts(5, 5), "bam_pairing_flags"),
    ],
)
def test_intrinsic_bam_errors_precede_cross_file_multiplicity_errors(
    tmp_path,
    score_count,
    bam_role_counts,
    expected_category,
):
    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        _insert_evidence(
            store,
            "read-1",
            score_count=score_count,
            score_line_numbers=(2,),
            bam_role_counts=bam_role_counts,
        )

        error = _reconciliation_error(store)

    assert error.category == expected_category


@pytest.mark.parametrize(
    ("read_ids", "expected_read_id"),
    [
        (("read-2", "read-10"), "read-10"),
        (("prefix!", "prefix"), "prefix"),
        (("a", "A"), "A"),
        (("é", "e\u0301"), "e\u0301"),
        (("read_a", "read-a"), "read-a"),
    ],
)
def test_lowest_failing_read_id_uses_exact_utf8_ordinal_order(
    tmp_path,
    read_ids,
    expected_read_id,
):
    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        for read_id in reversed(read_ids):
            _insert_evidence(store, read_id, score_count=1, score_line_numbers=(2,))

        error = _reconciliation_error(store)

    assert error.category == "score_only"
    assert error.read_id == expected_read_id


def test_shuffled_score_and_bam_streams_preserve_selected_error(tmp_path):
    score_path = tmp_path / "scores.tsv"
    bam_path = tmp_path / "reads.bam"
    score_ids = ["z-score-only", "y-score-only"]
    bam_ids = ["x-bam-only", "A-bam-only"]
    errors = []

    for current_score_ids, current_bam_ids in (
        (score_ids, bam_ids),
        (list(reversed(score_ids)), list(reversed(bam_ids))),
    ):
        _write_score_table(score_path, current_score_ids)
        _write_bam(bam_path, current_bam_ids)
        with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
            lyra._collect_score_evidence(store, score_path, "sample")
            lyra._collect_bam_evidence(store, bam_path)
            errors.append(_reconciliation_error(store, score_path, bam_path))

    attributes = (
        "category",
        "score_path",
        "bam_path",
        "read_id",
        "score_count",
        "eligible_bam_count",
        "score_line_numbers",
        "bam_role_counts",
        "reason",
    )
    assert [getattr(errors[0], name) for name in attributes] == [
        getattr(errors[1], name) for name in attributes
    ]
    assert errors[0].category == "bam_only"
    assert errors[0].read_id == "A-bam-only"


def test_diagnostic_contains_exact_bounded_context_and_converted_paths(tmp_path):
    score_path = tmp_path / ("unsafe\r\n" + ("s" * 1000) + ".tsv")
    bam_path = tmp_path / ("unsafe\r\n" + ("b" * 1000) + ".bam")
    read_id = "unsafe\r\n" + ("q" * 10000)
    role_counts = (1, 1, 1, 1, 1, 1, 1, 1)

    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        _insert_evidence(
            store,
            read_id,
            score_count=9,
            score_line_numbers=(7, 11, 13),
            bam_role_counts=role_counts,
        )

        error = _reconciliation_error(store, score_path, bam_path)

    assert error.category == "bam_record_limit"
    assert error.score_path == str(score_path)
    assert error.bam_path == str(bam_path)
    assert error.read_id == read_id
    assert error.score_count == 9
    assert error.eligible_bam_count == 8
    assert error.score_line_numbers == (7, 11)
    assert error.bam_role_counts == tuple(zip(_ROLE_LABELS, role_counts))
    assert "\r" not in str(error)
    assert "\n" not in str(error)
    assert "\\r" in str(error)
    assert "\\n" in str(error)
    assert len(str(error)) < 1000
