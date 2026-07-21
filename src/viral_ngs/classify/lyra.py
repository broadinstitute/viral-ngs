"""Strict validation helpers for native Lyra score tables."""

import codecs
import gzip
import io
import os
import re
import zlib
from dataclasses import dataclass
from decimal import Decimal

import zstandard as zstd

import viral_ngs.core.file as util_file


RENDERED_VALUE_CAP = 160
EXPECTED_SCORE_HEADER = ("read_id", "score", "call")

_SUPPORTED_SCORE_SUFFIXES = (".tsv", ".tsv.gz", ".tsv.zst")
_SCORE_PATTERN = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", re.ASCII)
_MIN_SCORE = Decimal("0")
_MAX_SCORE = Decimal("1")
_COMPRESSION_EXCEPTIONS = (
    EOFError,
    gzip.BadGzipFile,
    zlib.error,
    zstd.ZstdError,
)


def _bounded_repr(value):
    """Render a caller value without emitting controls or unbounded text."""
    if isinstance(value, os.PathLike):
        value = os.fspath(value)

    if isinstance(value, str):
        prefix = value[:RENDERED_VALUE_CAP]
        was_truncated = len(value) > len(prefix)
    elif isinstance(value, bytes):
        prefix = value[:RENDERED_VALUE_CAP]
        was_truncated = len(value) > len(prefix)
    else:
        prefix = value
        was_truncated = False

    rendered = repr(prefix)
    if was_truncated or len(rendered) > RENDERED_VALUE_CAP:
        return rendered[:RENDERED_VALUE_CAP - 3] + "..."
    return rendered


class LyraInputError(ValueError):
    """A native Lyra input failure with stable machine-readable context."""

    def __init__(
        self,
        category,
        path=None,
        line_number=None,
        field=None,
        reason=None,
        offending_value=None,
    ):
        self.category = category
        self.path = path
        self.line_number = line_number
        self.field = field
        self.reason = reason
        self.offending_value = (
            None if offending_value is None else _bounded_repr(offending_value)
        )

        details = ["category={}".format(_bounded_repr(category))]
        if path is not None:
            details.append("path={}".format(_bounded_repr(path)))
        if line_number is not None:
            details.append("line={}".format(line_number))
        if field is not None:
            details.append("field={}".format(_bounded_repr(field)))
        if reason is not None:
            details.append("reason={}".format(_bounded_repr(reason)))
        if self.offending_value is not None:
            details.append("offending_value={}".format(self.offending_value))
        super().__init__("Lyra input error: " + "; ".join(details))


@dataclass(frozen=True)
class LyraScoreRecord:
    """One validated native Lyra score row."""

    read_id: str
    score_text: str
    score: Decimal
    native_call: str
    line_number: int


def validate_sample_id(sample_id):
    """Return a valid sample ID unchanged or raise ``LyraInputError``."""
    if not isinstance(sample_id, str):
        raise LyraInputError(
            category="sample_id",
            field="sample_id",
            reason="sample ID must be a string",
            offending_value=sample_id,
        )
    if sample_id == "":
        raise LyraInputError(
            category="sample_id",
            field="sample_id",
            reason="sample ID must not be empty",
            offending_value=sample_id,
        )
    if sample_id.isspace():
        raise LyraInputError(
            category="sample_id",
            field="sample_id",
            reason="sample ID must not contain only whitespace",
            offending_value=sample_id,
        )
    if sample_id != sample_id.strip():
        raise LyraInputError(
            category="sample_id",
            field="sample_id",
            reason="sample ID must not have boundary whitespace",
            offending_value=sample_id,
        )
    if not sample_id.isprintable():
        raise LyraInputError(
            category="sample_id",
            field="sample_id",
            reason="sample ID must contain only printable characters",
            offending_value=sample_id,
        )
    try:
        sample_id.encode("utf-8", errors="strict")
    except UnicodeEncodeError as exc:
        raise LyraInputError(
            category="sample_id",
            field="sample_id",
            reason="sample ID must be UTF-8 encodable",
            offending_value=sample_id,
        ) from exc
    return sample_id


def iter_lyra_score_records(score_path, sample_id):
    """Return a lazy iterator over a strictly validated native score table."""
    validate_sample_id(sample_id)
    path = os.fspath(score_path)
    if not path.endswith(_SUPPORTED_SCORE_SUFFIXES):
        raise LyraInputError(
            category="extension",
            path=score_path,
            field="file",
            reason="score path must end with .tsv, .tsv.gz, or .tsv.zst",
            offending_value=path,
        )
    return _iter_lyra_score_records(path, score_path)


def _iter_lyra_score_records(path, caller_path):
    rows = _iter_physical_rows(path, caller_path)
    line_number, header_text = next(rows)
    _validate_header(header_text, caller_path, line_number)

    for line_number, row_text in rows:
        yield _parse_score_row(row_text, caller_path, line_number)


def _iter_physical_rows(path, caller_path):
    try:
        with util_file.open_or_gzopen(path, "rb") as raw_stream:
            if isinstance(raw_stream, io.BufferedIOBase):
                yield from _iter_buffered_rows(raw_stream, caller_path)
            else:
                with io.BufferedReader(raw_stream) as buffered_stream:
                    yield from _iter_buffered_rows(buffered_stream, caller_path)
        if path.endswith(".tsv.zst"):
            _validate_complete_zstd_frames(path)
    except LyraInputError as exc:
        if path.endswith(".tsv.zst") and exc.category == "empty_file":
            try:
                _validate_complete_zstd_frames(path)
            except _COMPRESSION_EXCEPTIONS as compression_exc:
                raise _compression_error(caller_path) from compression_exc
        raise
    except FileNotFoundError:
        raise
    except _COMPRESSION_EXCEPTIONS as exc:
        if not path.endswith((".tsv.gz", ".tsv.zst")):
            raise
        raise _compression_error(caller_path) from exc


def _validate_complete_zstd_frames(path):
    decompressor = None
    saw_complete_frame = False
    with open(path, "rb") as compressed_stream:
        for chunk in iter(lambda: compressed_stream.read(64 * 1024), b""):
            pending = chunk
            while pending:
                if decompressor is None:
                    decompressor = zstd.ZstdDecompressor().decompressobj()
                decompressor.decompress(pending)
                pending = decompressor.unused_data
                if decompressor.eof:
                    saw_complete_frame = True
                    decompressor = None
                elif pending:
                    raise zstd.ZstdError("invalid trailing Zstandard frame data")
                else:
                    break

    if decompressor is not None or not saw_complete_frame:
        raise EOFError("truncated Zstandard frame")


def _compression_error(caller_path):
    return LyraInputError(
        category="compression",
        path=caller_path,
        field="file",
        reason="compressed score table could not be decoded",
    )


def _iter_buffered_rows(stream, caller_path):
    saw_line = False
    for line_number, raw_line in enumerate(stream, start=1):
        saw_line = True
        content = raw_line
        if content.endswith(b"\n"):
            content = content[:-1]
            if content.endswith(b"\r"):
                content = content[:-1]

        if b"\r" in content:
            raise LyraInputError(
                category="line_ending",
                path=caller_path,
                line_number=line_number,
                field="row",
                reason="bare carriage return is not allowed",
                offending_value=content,
            )

        if line_number == 1 and content.startswith(codecs.BOM_UTF8):
            content = content[len(codecs.BOM_UTF8):]

        try:
            text = content.decode("utf-8", errors="strict")
        except UnicodeDecodeError as exc:
            raise LyraInputError(
                category="utf8",
                path=caller_path,
                line_number=line_number,
                field=_field_at_byte_offset(content, exc.start),
                reason="row is not valid UTF-8",
                offending_value=content,
            ) from exc

        if text == "":
            raise LyraInputError(
                category="blank_row",
                path=caller_path,
                line_number=line_number,
                field="row",
                reason="blank physical rows are not allowed",
                offending_value=text,
            )
        yield line_number, text

    if not saw_line:
        raise LyraInputError(
            category="empty_file",
            path=caller_path,
            field="file",
            reason="score table has no header",
        )


def _field_at_byte_offset(content, offset):
    field_index = min(content[:offset].count(b"\t"), len(EXPECTED_SCORE_HEADER) - 1)
    return EXPECTED_SCORE_HEADER[field_index]


def _validate_header(header_text, caller_path, line_number):
    fields = tuple(header_text.split("\t"))
    if len(fields) != len(EXPECTED_SCORE_HEADER):
        field = (
            EXPECTED_SCORE_HEADER[len(fields)]
            if len(fields) < len(EXPECTED_SCORE_HEADER)
            else "header"
        )
        raise LyraInputError(
            category="header_width",
            path=caller_path,
            line_number=line_number,
            field=field,
            reason="header must contain exactly three fields",
            offending_value=header_text,
        )

    for expected, actual in zip(EXPECTED_SCORE_HEADER, fields):
        if actual != expected:
            raise LyraInputError(
                category="header_schema",
                path=caller_path,
                line_number=line_number,
                field=expected,
                reason="header field does not match the native Lyra schema",
                offending_value=actual,
            )


def _parse_score_row(row_text, caller_path, line_number):
    fields = row_text.split("\t")
    if len(fields) != len(EXPECTED_SCORE_HEADER):
        raise LyraInputError(
            category="row_width",
            path=caller_path,
            line_number=line_number,
            field="row",
            reason="score row must contain exactly three fields",
            offending_value=row_text,
        )

    read_id, score_text, native_call = fields
    for field, value in zip(EXPECTED_SCORE_HEADER, fields):
        if value != value.strip():
            raise LyraInputError(
                category="field_whitespace",
                path=caller_path,
                line_number=line_number,
                field=field,
                reason="field must not have leading or trailing whitespace",
                offending_value=value,
            )

    if read_id == "":
        raise LyraInputError(
            category="empty_read_id",
            path=caller_path,
            line_number=line_number,
            field="read_id",
            reason="read ID must not be empty",
            offending_value=read_id,
        )

    if _SCORE_PATTERN.fullmatch(score_text) is None:
        raise LyraInputError(
            category="score_syntax",
            path=caller_path,
            line_number=line_number,
            field="score",
            reason="score does not use the canonical ASCII decimal syntax",
            offending_value=score_text,
        )

    score = Decimal(score_text)
    if not score.is_finite() or not _MIN_SCORE <= score <= _MAX_SCORE:
        raise LyraInputError(
            category="score_range",
            path=caller_path,
            line_number=line_number,
            field="score",
            reason="score must be finite and between zero and one inclusive",
            offending_value=score_text,
        )

    if native_call not in {"0", "1"}:
        raise LyraInputError(
            category="call",
            path=caller_path,
            line_number=line_number,
            field="call",
            reason="native call must be exactly 0 or 1",
            offending_value=native_call,
        )

    return LyraScoreRecord(
        read_id=read_id,
        score_text=score_text,
        score=score,
        native_call=native_call,
        line_number=line_number,
    )
