import codecs
from dataclasses import FrozenInstanceError
from decimal import Decimal
import inspect
import re
import zlib

import pytest

from viral_ngs.classify import lyra
import viral_ngs.core.file as util_file


_HEADER = b"read_id\tscore\tcall"
_SCORE_SUFFIXES = (".tsv", ".tsv.gz", ".tsv.zst")


def _write_plain_score_file(tmp_path, contents, name="scores.tsv"):
    path = tmp_path / name
    path.write_bytes(contents)
    return path


def _write_score_file(tmp_path, suffix, contents, stem="scores"):
    path = tmp_path / (stem + suffix)
    with util_file.open_or_gzopen(str(path), "wb") as outf:
        outf.write(contents)
    return path


def _assert_input_error(
    exc_info,
    *,
    category,
    path,
    line_number,
    field,
):
    error = exc_info.value
    assert error.category == category
    assert error.path == str(path)
    assert error.line_number == line_number
    assert error.field == field
    assert error.reason
    return error


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


@pytest.mark.parametrize(
    "contents",
    [
        _HEADER + b"\nread 1\t0.8\t1\n",
        _HEADER + b"\r\nread 1\t0.8\t1\r\n",
        _HEADER + b"\nread 1\t0.8\t1",
        codecs.BOM_UTF8 + _HEADER + b"\nread 1\t0.8\t1\n",
    ],
)
def test_valid_newline_and_bom_forms_preserve_records(tmp_path, contents):
    path = _write_plain_score_file(tmp_path, contents)

    records = list(lyra.iter_lyra_score_records(str(path), "sample A"))

    assert records == [
        lyra.LyraScoreRecord(
            read_id="read 1",
            score_text="0.8",
            score=Decimal("0.8"),
            native_call="1",
            line_number=2,
        )
    ]


def test_valid_scores_retain_each_canonical_spelling(tmp_path):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER
        + b"\nread-0\t0\t0"
        + b"\nread-08\t0.8\t1"
        + b"\nread-0800000\t0.800000\t1"
        + b"\nread-1\t1\t1"
        + b"\nread-10\t1.0\t1\n",
    )

    records = list(lyra.iter_lyra_score_records(str(path), "sample"))

    assert [record.score_text for record in records] == [
        "0",
        "0.8",
        "0.800000",
        "1",
        "1.0",
    ]
    assert [record.score for record in records] == [
        Decimal("0"),
        Decimal("0.8"),
        Decimal("0.800000"),
        Decimal("1"),
        Decimal("1.0"),
    ]


def test_valid_header_only_table_yields_no_records(tmp_path):
    path = _write_plain_score_file(tmp_path, _HEADER + b"\n")

    assert list(lyra.iter_lyra_score_records(str(path), "sample")) == []


@pytest.mark.parametrize(
    ("header", "category", "field"),
    [
        (b"read_id\tscore", "header_width", "call"),
        (b"read_id\tscore\tcall\textra", "header_width", "header"),
        (b"score\tread_id\tcall", "header_schema", "read_id"),
        (b"read\tscore\tcall", "header_schema", "read_id"),
        (b"read_id\tcall\tscore", "header_schema", "score"),
        (b" read_id\tscore\tcall", "header_schema", "read_id"),
        (b"read_id\tscore \tcall", "header_schema", "score"),
        (b"read_id\tscore\tcall ", "header_schema", "call"),
        (codecs.BOM_UTF8 + codecs.BOM_UTF8 + _HEADER, "header_schema", "read_id"),
    ],
)
def test_header_contract_is_exact(tmp_path, header, category, field):
    path = _write_plain_score_file(tmp_path, header + b"\n")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category=category,
        path=path,
        line_number=1,
        field=field,
    )


@pytest.mark.parametrize(
    ("contents", "line_number"),
    [
        (b"\n" + _HEADER + b"\n", 1),
        (_HEADER + b"\nread-1\t0.8\t1\n\nread-2\t0.9\t1\n", 3),
        (_HEADER + b"\nread-1\t0.8\t1\n\n", 3),
    ],
)
def test_blank_rows_fail_on_their_physical_line(tmp_path, contents, line_number):
    path = _write_plain_score_file(tmp_path, contents)

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="blank_row",
        path=path,
        line_number=line_number,
        field="row",
    )


@pytest.mark.parametrize(
    ("contents", "line_number"),
    [
        (_HEADER + b"\rread-1\t0.8\t1\r", 1),
        (_HEADER + b"\nread-1\t0.8\r1\n", 2),
        (_HEADER + b"\nread-1\t0.8\t1\r\r\n", 2),
    ],
)
def test_newline_contract_rejects_every_bare_cr(tmp_path, contents, line_number):
    path = _write_plain_score_file(tmp_path, contents)

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="line_ending",
        path=path,
        line_number=line_number,
        field="row",
    )


@pytest.mark.parametrize(
    ("row", "field"),
    [
        (b"read-\xff\t0.8\t1", "read_id"),
        (b"read-1\t0.\xff\t1", "score"),
        (b"read-1\t0.8\t\xff", "call"),
    ],
)
def test_utf8_errors_identify_the_field(tmp_path, row, field):
    path = _write_plain_score_file(tmp_path, _HEADER + b"\n" + row + b"\n")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    error = _assert_input_error(
        exc_info,
        category="utf8",
        path=path,
        line_number=2,
        field=field,
    )
    assert isinstance(error.__cause__, UnicodeDecodeError)


@pytest.mark.parametrize(
    "row",
    [
        b"read-1\t0.8",
        b"read-1\t0.8\t1\textra",
    ],
)
def test_row_width_is_exact(tmp_path, row):
    path = _write_plain_score_file(tmp_path, _HEADER + b"\n" + row + b"\n")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="row_width",
        path=path,
        line_number=2,
        field="row",
    )


@pytest.mark.parametrize(
    ("row", "field"),
    [
        (b" read-1\t0.8\t1", "read_id"),
        (b"read-1 \t0.8\t1", "read_id"),
        (b"read-1\t 0.8\t1", "score"),
        (b"read-1\t0.8 \t1", "score"),
        (b"read-1\t0.8\t 1", "call"),
        (b"read-1\t0.8\t1 ", "call"),
    ],
)
def test_row_boundary_whitespace_identifies_the_field(tmp_path, row, field):
    path = _write_plain_score_file(tmp_path, _HEADER + b"\n" + row + b"\n")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="field_whitespace",
        path=path,
        line_number=2,
        field=field,
    )


def test_empty_read_id_is_rejected_after_whitespace_checks(tmp_path):
    path = _write_plain_score_file(tmp_path, _HEADER + b"\n\tbad-score\tbad-call\n")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="empty_read_id",
        path=path,
        line_number=2,
        field="read_id",
    )


@pytest.mark.parametrize(
    "score",
    [
        b"",
        b"00",
        b"00.80",
        b"01.0",
        b".8",
        b"1.",
        b"8e-1",
        b"+0.8",
        b"-0.1",
        b"NaN",
        b"Inf",
        b"0_8",
        b"0,8",
        b"True",
        "٠.٨".encode("utf-8"),
    ],
)
def test_score_syntax_rejects_noncanonical_spellings(tmp_path, score):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\nread-1\t" + score + b"\tbad-call\n",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="score_syntax",
        path=path,
        line_number=2,
        field="score",
    )


@pytest.mark.parametrize("score", [b"2", b"1.0000001", b"999999999999"])
def test_score_range_is_exact_and_inclusive(tmp_path, score):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\nread-1\t" + score + b"\t1\n",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="score_range",
        path=path,
        line_number=2,
        field="score",
    )


@pytest.mark.parametrize(
    "native_call",
    [b"", b"2", b"-1", b"0.0", b"1.0", b"True", b"False"],
)
def test_call_accepts_only_exact_zero_or_one(tmp_path, native_call):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\nread-1\t0.8\t" + native_call + b"\n",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="call",
        path=path,
        line_number=2,
        field="call",
    )


def test_ordering_row_width_precedes_field_validation(tmp_path):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\n read-1\tbad-score\tbad-call\textra\n",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    assert exc_info.value.category == "row_width"


def test_ordering_score_precedes_call(tmp_path):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\nread-1\tbad-score\tbad-call\n",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    assert exc_info.value.category == "score_syntax"


def test_ordering_invalid_sample_precedes_extension_and_file_access(tmp_path):
    missing_path = tmp_path / "missing.bad"

    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra.iter_lyra_score_records(str(missing_path), " invalid")

    assert exc_info.value.category == "sample_id"


def test_extension_allowlist_is_case_sensitive(tmp_path):
    path = tmp_path / "scores.TSV"

    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra.iter_lyra_score_records(str(path), "sample")

    _assert_input_error(
        exc_info,
        category="extension",
        path=path,
        line_number=None,
        field="file",
    )


def test_empty_file_has_context_without_an_invented_line_number(tmp_path):
    path = _write_plain_score_file(tmp_path, b"")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="empty_file",
        path=path,
        line_number=None,
        field="file",
    )


def test_missing_file_remains_distinguishable(tmp_path):
    path = tmp_path / "missing.tsv"

    with pytest.raises(FileNotFoundError):
        list(lyra.iter_lyra_score_records(str(path), "sample"))


def test_valid_iterator_is_lazy_after_each_record(tmp_path):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\nread-1\t0.8\t1\nread-2\tbad-score\t1\n",
    )

    records = lyra.iter_lyra_score_records(str(path), "sample")

    assert not inspect.isgeneratorfunction(lyra.iter_lyra_score_records)
    assert next(records).read_id == "read-1"
    with pytest.raises(lyra.LyraInputError) as exc_info:
        next(records)
    assert exc_info.value.category == "score_syntax"
    assert exc_info.value.line_number == 3


def test_compression_parity_produces_equal_exact_records(tmp_path):
    contents = (
        codecs.BOM_UTF8
        + _HEADER
        + b"\r\nread-1\t0.800000\t1\r\nread-2\t0\t0"
    )
    record_sets = []
    for suffix in _SCORE_SUFFIXES:
        path = _write_score_file(tmp_path, suffix, contents)
        record_sets.append(
            list(lyra.iter_lyra_score_records(str(path), "sample A"))
        )

    assert record_sets[0] == record_sets[1] == record_sets[2]
    assert record_sets[0] == [
        lyra.LyraScoreRecord(
            read_id="read-1",
            score_text="0.800000",
            score=Decimal("0.800000"),
            native_call="1",
            line_number=2,
        ),
        lyra.LyraScoreRecord(
            read_id="read-2",
            score_text="0",
            score=Decimal("0"),
            native_call="0",
            line_number=3,
        ),
    ]


@pytest.mark.parametrize("suffix", _SCORE_SUFFIXES)
def test_compression_parity_accepts_header_only_table(tmp_path, suffix):
    path = _write_score_file(tmp_path, suffix, _HEADER + b"\n")

    assert list(lyra.iter_lyra_score_records(str(path), "sample")) == []


@pytest.mark.parametrize("suffix", _SCORE_SUFFIXES)
def test_compression_parity_rejects_logically_empty_file(tmp_path, suffix):
    path = _write_score_file(tmp_path, suffix, b"")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category="empty_file",
        path=path,
        line_number=None,
        field="file",
    )


@pytest.mark.parametrize("suffix", _SCORE_SUFFIXES)
@pytest.mark.parametrize(
    ("contents", "category", "line_number", "field"),
    [
        (
            codecs.BOM_UTF8 + codecs.BOM_UTF8 + _HEADER + b"\n",
            "header_schema",
            1,
            "read_id",
        ),
        (_HEADER + b"\nread-1\t0.8\r1\n", "line_ending", 2, "row"),
        (_HEADER + b"\n\n", "blank_row", 2, "row"),
        (_HEADER + b"\nread-\xff\t0.8\t1\n", "utf8", 2, "read_id"),
        (_HEADER + b"\nread-1\t0.\xff\t1\n", "utf8", 2, "score"),
        (_HEADER + b"\nread-1\t0.8\t\xff\n", "utf8", 2, "call"),
        (b"score\tread_id\tcall\n", "header_schema", 1, "read_id"),
        (_HEADER + b"\nread-1\t0.8\n", "row_width", 2, "row"),
        (
            _HEADER + b"\nread-1\t 0.8\tbad-call\n",
            "field_whitespace",
            2,
            "score",
        ),
        (_HEADER + b"\n\t0.8\t1\n", "empty_read_id", 2, "read_id"),
        (_HEADER + b"\nread-1\t00\t1\n", "score_syntax", 2, "score"),
        (_HEADER + b"\nread-1\t00.80\t1\n", "score_syntax", 2, "score"),
        (_HEADER + b"\nread-1\t01.0\t1\n", "score_syntax", 2, "score"),
        (_HEADER + b"\nread-1\t2\t1\n", "score_range", 2, "score"),
        (_HEADER + b"\nread-1\t0.8\tTrue\n", "call", 2, "call"),
    ],
)
def test_compression_parity_preserves_malformed_input_context(
    tmp_path,
    suffix,
    contents,
    category,
    line_number,
    field,
):
    path = _write_score_file(tmp_path, suffix, contents)

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    _assert_input_error(
        exc_info,
        category=category,
        path=path,
        line_number=line_number,
        field=field,
    )


@pytest.mark.parametrize("suffix", [".tsv.gz", ".tsv.zst"])
@pytest.mark.parametrize("raw_bytes", [b"not compressed", b"\x00\xffcorrupt"])
def test_compression_errors_wrap_mislabeled_or_corrupt_files(
    tmp_path,
    suffix,
    raw_bytes,
):
    path = tmp_path / ("scores" + suffix)
    path.write_bytes(raw_bytes)

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    error = _assert_input_error(
        exc_info,
        category="compression",
        path=path,
        line_number=None,
        field="file",
    )
    assert error.__cause__ is not None


def test_compression_errors_wrap_invalid_gzip_deflate_block(tmp_path):
    path = tmp_path / "scores.tsv.gz"
    path.write_bytes(bytes.fromhex("1f8b08000000000000ff07"))

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    error = _assert_input_error(
        exc_info,
        category="compression",
        path=path,
        line_number=None,
        field="file",
    )
    assert isinstance(error.__cause__, zlib.error)


@pytest.mark.parametrize(
    "filesystem_error",
    [PermissionError("permission denied"), OSError("filesystem failure")],
)
def test_compressed_filesystem_errors_remain_distinguishable(
    tmp_path,
    monkeypatch,
    filesystem_error,
):
    path = tmp_path / "scores.tsv.gz"

    def raise_filesystem_error(*args, **kwargs):
        raise filesystem_error

    monkeypatch.setattr(util_file, "open_or_gzopen", raise_filesystem_error)

    with pytest.raises(type(filesystem_error)) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    assert exc_info.value is filesystem_error


@pytest.mark.parametrize("suffix", [".tsv.gz", ".tsv.zst"])
def test_compression_errors_wrap_truncated_frames(tmp_path, suffix):
    path = _write_score_file(
        tmp_path,
        suffix,
        _HEADER + b"\nread-1\t0.800000\t1\n" + (b"x" * 10000),
    )
    compressed = path.read_bytes()
    path.write_bytes(compressed[: max(1, len(compressed) // 2)])

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    error = _assert_input_error(
        exc_info,
        category="compression",
        path=path,
        line_number=None,
        field="file",
    )
    assert error.__cause__ is not None


def test_diagnostic_escapes_control_path_and_caps_huge_field(tmp_path):
    path = _write_plain_score_file(
        tmp_path,
        _HEADER + b"\nread-1\t" + (b"9" * 1000000) + b"x\t1\n",
        name="unsafe\npath.tsv",
    )

    with pytest.raises(lyra.LyraInputError) as exc_info:
        list(lyra.iter_lyra_score_records(str(path), "sample"))

    error = exc_info.value
    assert error.category == "score_syntax"
    assert len(error.offending_value) <= lyra.RENDERED_VALUE_CAP
    assert "\n" not in str(error)
    assert "\\n" in str(error)
    assert len(str(error)) < 600


def test_score_regex_is_one_fixed_ascii_linear_pattern():
    assert lyra._SCORE_PATTERN.pattern == r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?"
    assert lyra._SCORE_PATTERN.flags & re.ASCII


@pytest.mark.parametrize("suffix", _SCORE_SUFFIXES)
def test_compression_iterator_does_not_read_later_bad_row_early(tmp_path, suffix):
    path = _write_score_file(
        tmp_path,
        suffix,
        _HEADER + b"\nread-1\t0.8\t1\nread-2\tbad-score\t1\n",
    )

    records = lyra.iter_lyra_score_records(str(path), "sample")

    assert next(records).read_id == "read-1"
    with pytest.raises(lyra.LyraInputError) as exc_info:
        next(records)
    assert exc_info.value.category == "score_syntax"
    assert exc_info.value.line_number == 3
