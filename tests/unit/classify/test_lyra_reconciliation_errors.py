import inspect
import itertools

import pytest

from viral_ngs.classify import lyra


_ROLE_LABELS = ("U--", "U1-", "U-2", "U12", "P--", "P1-", "P-2", "P12")


def _role_counts(*role_indexes):
    return tuple(role_indexes.count(index) for index in range(len(_ROLE_LABELS)))


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
