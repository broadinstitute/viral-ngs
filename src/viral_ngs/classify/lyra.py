"""Strict validation helpers for native Lyra score tables."""

import os
from dataclasses import dataclass
from decimal import Decimal


RENDERED_VALUE_CAP = 160


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
