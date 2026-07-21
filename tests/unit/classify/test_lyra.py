from dataclasses import FrozenInstanceError
from decimal import Decimal

import pytest

from viral_ngs.classify import lyra


def test_validate_sample_id_preserves_valid_values():
    sample_ids = [
        "sample A",
        "échantillon-β",
        "é",
        "e\u0301",
        "Sample",
        "sample",
        "x" * 10000,
    ]

    for sample_id in sample_ids:
        assert lyra.validate_sample_id(sample_id) is sample_id


@pytest.mark.parametrize(
    "sample_id",
    [
        None,
        1,
        "",
        " ",
        " sample",
        "sample ",
        "sample\tname",
        "sample\rname",
        "sample\nname",
        "sample\x00name",
        "sample\x85name",
        "sample\u200bname",
        "sample\ud800name",
    ],
)
def test_validate_sample_id_rejects_unsafe_values(sample_id):
    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra.validate_sample_id(sample_id)

    error = exc_info.value
    assert error.category == "sample_id"
    assert error.path is None
    assert error.line_number is None
    assert error.field == "sample_id"
    assert error.reason
    if sample_id is not None:
        assert error.offending_value is not None


def test_input_error_escapes_and_caps_caller_values():
    path = "unsafe\npath-" + ("p" * 10000)
    value = "unsafe\r\nvalue-" + ("v" * 10000)

    error = lyra.LyraInputError(
        category="sample_id",
        path=path,
        line_number=None,
        field="sample_id",
        reason="invalid sample ID",
        offending_value=value,
    )

    assert error.path is path
    assert len(error.offending_value) <= lyra.RENDERED_VALUE_CAP
    assert "\n" not in str(error)
    assert "\r" not in str(error)
    assert "\\n" in str(error)
    assert "\\r" in str(error)
    assert ("p" * 1000) not in str(error)
    assert ("v" * 1000) not in str(error)


def test_score_record_is_immutable_and_retains_exact_decimal_text():
    record = lyra.LyraScoreRecord(
        read_id="read-1",
        score_text="0.800000",
        score=Decimal("0.800000"),
        native_call="1",
        line_number=2,
    )

    assert record.score_text == "0.800000"
    assert record.score == Decimal("0.800000")
    with pytest.raises(FrozenInstanceError):
        record.score_text = "0.8"
