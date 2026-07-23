"""Strict validation helpers for native Lyra score tables."""

import codecs
import copy
import gzip
import io
import os
import re
import secrets
import sqlite3
import stat
import subprocess
import tempfile
import zlib
from contextlib import contextmanager
from dataclasses import dataclass
from decimal import Decimal
from decimal import InvalidOperation

import pysam
import zstandard as zstd

import viral_ngs.core.file as util_file
import viral_ngs.core.misc as util_misc


RENDERED_VALUE_CAP = 160
EXPECTED_SCORE_HEADER = ("read_id", "score", "call")
SQLITE_COMMIT_INTERVAL = 10000
PAIRING_SINGLE_END = "Single-end"
PAIRING_COMPLETE = "Paired-complete"
PAIRING_INCOMPLETE = "Paired-incomplete"
CALL_VIRAL = "Viral"
CALL_NONVIRAL = "Non-viral"
NORMALIZED_HEADER = (
    "SAMPLE_ID",
    "READ_ID",
    "LYRA_N_SCORES",
    "LYRA_PAIRING",
    "LYRA_MIN_SCORE",
    "LYRA_MAX_SCORE",
    "LYRA_THRESHOLD",
    "LYRA_CALL",
)
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

_SUPPORTED_SCORE_SUFFIXES = (".tsv", ".tsv.gz", ".tsv.zst")
_MAX_CANONICAL_THRESHOLD_LENGTH = 160
_MAX_PUBLICATION_CLEANUP_FAILURES = 16
_MAX_PUBLICATION_COMMAND_FIELDS = 8
_SCORE_PATTERN = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", re.ASCII)
_MIN_SCORE = Decimal("0")
_MAX_SCORE = Decimal("1")
_COMPRESSION_EXCEPTIONS = (
    EOFError,
    gzip.BadGzipFile,
    zlib.error,
    zstd.ZstdError,
)
_SCORE_EVIDENCE_UPSERT = """
    INSERT INTO evidence (
        read_id_key,
        score_count,
        score_1_text,
        score_1_line
    ) VALUES (?, ?, ?, ?)
    ON CONFLICT(read_id_key) DO UPDATE SET
        score_1_text = CASE
            WHEN evidence.score_count = 0 THEN excluded.score_1_text
            ELSE evidence.score_1_text
        END,
        score_1_line = CASE
            WHEN evidence.score_count = 0 THEN excluded.score_1_line
            ELSE evidence.score_1_line
        END,
        score_2_text = CASE
            WHEN evidence.score_count = 1 THEN excluded.score_1_text
            ELSE evidence.score_2_text
        END,
        score_2_line = CASE
            WHEN evidence.score_count = 1 THEN excluded.score_1_line
            ELSE evidence.score_2_line
        END,
        score_count = evidence.score_count + excluded.score_count
"""
_BAM_EVIDENCE_UPSERT = """
    INSERT INTO evidence (
        read_id_key,
        eligible_bam_count,
        unpaired_none_count,
        unpaired_r1_count,
        unpaired_r2_count,
        unpaired_both_count,
        paired_none_count,
        paired_r1_count,
        paired_r2_count,
        paired_both_count
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(read_id_key) DO UPDATE SET
        eligible_bam_count = (
            evidence.eligible_bam_count + excluded.eligible_bam_count
        ),
        unpaired_none_count = (
            evidence.unpaired_none_count + excluded.unpaired_none_count
        ),
        unpaired_r1_count = (
            evidence.unpaired_r1_count + excluded.unpaired_r1_count
        ),
        unpaired_r2_count = (
            evidence.unpaired_r2_count + excluded.unpaired_r2_count
        ),
        unpaired_both_count = (
            evidence.unpaired_both_count + excluded.unpaired_both_count
        ),
        paired_none_count = (
            evidence.paired_none_count + excluded.paired_none_count
        ),
        paired_r1_count = (
            evidence.paired_r1_count + excluded.paired_r1_count
        ),
        paired_r2_count = (
            evidence.paired_r2_count + excluded.paired_r2_count
        ),
        paired_both_count = (
            evidence.paired_both_count + excluded.paired_both_count
        )
"""
_BAM_ROLE_INDEX = {
    (False, False, False): 0,
    (False, True, False): 1,
    (False, False, True): 2,
    (False, True, True): 3,
    (True, False, False): 4,
    (True, True, False): 5,
    (True, False, True): 6,
    (True, True, True): 7,
}
_BAM_ROLE_LABELS = (
    "U--",
    "U1-",
    "U-2",
    "U12",
    "P--",
    "P1-",
    "P-2",
    "P12",
)
_BAM_ROLE_COLUMNS = (
    "unpaired_none_count",
    "unpaired_r1_count",
    "unpaired_r2_count",
    "unpaired_both_count",
    "paired_none_count",
    "paired_r1_count",
    "paired_r2_count",
    "paired_both_count",
)
_SINGLE_END_ROLE_COUNTS = (1, 0, 0, 0, 0, 0, 0, 0)
_PAIRED_R1_ROLE_COUNTS = (0, 0, 0, 0, 0, 1, 0, 0)
_PAIRED_R2_ROLE_COUNTS = (0, 0, 0, 0, 0, 0, 1, 0)
_COMPLETE_PAIR_ROLE_COUNTS = (0, 0, 0, 0, 0, 1, 1, 0)
_RECONCILIATION_REASONS = {
    "bam_record_limit": "more than two eligible BAM records share this read ID",
    "bam_pairing_flags": "eligible BAM records have incoherent pairing flags",
    "score_only": "score records have no matching eligible BAM records",
    "bam_only": "eligible BAM records have no matching score records",
    "record_count_mismatch": "score and eligible BAM record counts differ",
}


def _bounded_repr(value):
    """Render a caller value without emitting controls or unbounded text."""
    if type(value) is str:
        prefix = value[:RENDERED_VALUE_CAP]
        was_truncated = len(value) > len(prefix)
    elif type(value) is bytes:
        prefix = value[:RENDERED_VALUE_CAP]
        was_truncated = len(value) > len(prefix)
    else:
        return "<non-string value>"

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


class LyraPathError(LyraInputError):
    """A pre-processing path conflict with stable bounded context."""

    def __init__(
        self,
        category,
        role,
        path,
        reason,
        conflicting_role=None,
        conflicting_path=None,
        stage=None,
    ):
        self.category = category
        self.role = role
        self.path = path
        self.reason = reason
        self.conflicting_role = conflicting_role
        self.conflicting_path = conflicting_path
        self.stage = stage

        details = [
            "category={}".format(_bounded_repr(category)),
            "role={}".format(_bounded_repr(role)),
            "path={}".format(_bounded_repr(path)),
            "reason={}".format(_bounded_repr(reason)),
        ]
        if conflicting_role is not None:
            details.append(
                "conflicting_role={}".format(_bounded_repr(conflicting_role))
            )
        if conflicting_path is not None:
            details.append(
                "conflicting_path={}".format(_bounded_repr(conflicting_path))
            )
        if stage is not None:
            details.append("stage={}".format(_bounded_repr(stage)))
        ValueError.__init__(self, "Lyra path error: " + "; ".join(details))


class LyraReconciliationError(LyraInputError):
    """A whole-store reconciliation failure with bounded evidence context."""

    def __init__(
        self,
        category,
        score_path,
        bam_path,
        read_id,
        score_count,
        eligible_bam_count,
        score_line_numbers,
        bam_role_counts,
        reason,
    ):
        self.category = category
        self.score_path = score_path
        self.bam_path = bam_path
        self.read_id = read_id
        self.score_count = score_count
        self.eligible_bam_count = eligible_bam_count
        self.score_line_numbers = tuple(score_line_numbers)[:2]
        self.bam_role_counts = tuple(
            (label, count)
            for label, count in zip(_BAM_ROLE_LABELS, bam_role_counts)
            if count
        )
        self.reason = reason

        details = [
            "category={}".format(_bounded_repr(category)),
            "score_path={}".format(_bounded_repr(score_path)),
            "bam_path={}".format(_bounded_repr(bam_path)),
            "read_id={}".format(_bounded_repr(read_id)),
            "score_count={}".format(score_count),
            "eligible_bam_count={}".format(eligible_bam_count),
            "score_line_numbers={}".format(self.score_line_numbers),
            "bam_role_counts={}".format(self.bam_role_counts),
            "reason={}".format(_bounded_repr(reason)),
        ]
        ValueError.__init__(self, "Lyra reconciliation error: " + "; ".join(details))


class LyraArtifactConsistencyError(RuntimeError):
    """An impossible finalized-artifact mismatch with bounded scalar context."""

    def __init__(
        self,
        category,
        field=None,
        expected=None,
        actual=None,
        read_id=None,
    ):
        self.category = category
        self.field = field
        self.expected = expected
        self.actual = actual
        self.read_id = read_id

        details = ["category={}".format(_bounded_repr(category))]
        for name, value in (
            ("field", field),
            ("expected", expected),
            ("actual", actual),
            ("read_id", read_id),
        ):
            if value is not None:
                details.append(
                    "{}={}".format(name, _bounded_repr(str(value)))
                )
        super().__init__("Lyra artifact consistency error: " + "; ".join(details))


@dataclass(frozen=True)
class FileIdentity:
    """Stable scalar identity facts for one opened file."""

    device: int
    inode: int
    size: int
    mtime_ns: int


@dataclass(frozen=True)
class ArtifactDestination:
    """One resolved, immutable artifact destination."""

    role: str
    caller_path: str
    final_path: str
    parent_path: str
    destination_key: tuple
    suffix: str


@dataclass(frozen=True)
class ParentDirectoryAnchor:
    """One planned output parent retained by descriptor and identity."""

    parent_path: str
    identity: tuple
    descriptor: int


@dataclass(frozen=True)
class ArtifactStage:
    """One hidden same-parent artifact under transaction ownership."""

    role: str
    basename: str
    display_path: str
    descriptor: int
    object_identity: tuple
    destination: ArtifactDestination
    parent: ParentDirectoryAnchor

    @property
    def stage_path(self):
        """Retain the former diagnostic-path spelling for direct callers."""
        return self.display_path


@dataclass(frozen=True)
class PublishedArtifact:
    """One final directory entry verified as owned by this transaction."""

    role: str
    final_path: str
    identity: FileIdentity


@dataclass(frozen=True)
class CleanupFailure:
    """One immutable bounded publication-cleanup diagnostic fact."""

    operation: str
    role: str = None
    path: str = None
    error_type: str = None
    errno: int = None
    category: str = None


@dataclass(frozen=True)
class CleanupOutcome:
    """One bounded rollback or stage-cleanup result."""

    operation: str
    role: str
    path: str
    status: str
    expected: FileIdentity = None
    actual: FileIdentity = None
    error_type: str = None
    errno: int = None


@dataclass(frozen=True)
class LyraArtifactPathPlan:
    """Converted inputs and resolved destinations for one artifact set."""

    score_path: str
    bam_path: str
    normalized: ArtifactDestination
    summary: ArtifactDestination
    viral_bam: ArtifactDestination


def _file_identity(stat_result):
    return FileIdentity(
        device=stat_result.st_dev,
        inode=stat_result.st_ino,
        size=stat_result.st_size,
        mtime_ns=stat_result.st_mtime_ns,
    )


def _render_file_identity(identity):
    if identity is None:
        return "none"
    return "device={},inode={},size={},mtime_ns={}".format(
        identity.device,
        identity.inode,
        identity.size,
        identity.mtime_ns,
    )


class LyraSourceIdentityError(RuntimeError):
    """A source BAM replacement or mutation with bounded identity facts."""

    category = "source_bam_identity"

    def __init__(
        self,
        stage,
        path,
        expected,
        actual=None,
        actual_status=None,
    ):
        self.stage = stage
        self.path = path
        self.expected = expected
        self.actual = actual
        self.actual_status = actual_status

        details = [
            "category={}".format(_bounded_repr(self.category)),
            "stage={}".format(_bounded_repr(stage)),
            "path={}".format(_bounded_repr(path)),
            "expected={}".format(_render_file_identity(expected)),
        ]
        if actual is not None:
            details.append("actual={}".format(_render_file_identity(actual)))
        if actual_status is not None:
            details.append(
                "actual_status={}".format(_bounded_repr(actual_status))
            )
        super().__init__("Lyra source identity error: " + "; ".join(details))


def _known_primary_category(primary):
    if isinstance(primary, LyraSourceIdentityError):
        return LyraSourceIdentityError.category
    if isinstance(
        primary,
        (
            LyraInputError,
            LyraArtifactConsistencyError,
        ),
    ):
        category = vars(primary).get("category")
        if type(category) is str:
            return category
    return None


def _known_primary_path(primary):
    path = None
    if isinstance(
        primary,
        (
            LyraInputError,
            LyraArtifactConsistencyError,
            LyraSourceIdentityError,
        ),
    ):
        path = vars(primary).get("path")
    elif isinstance(primary, OSError):
        path = vars(primary).get("filename")
    if type(path) in (str, bytes):
        return _bounded_repr(path)
    return None


def _known_primary_command(primary):
    if not isinstance(primary, subprocess.CalledProcessError):
        return (), False
    command = vars(primary).get("cmd")
    if type(command) in (str, bytes):
        return (_bounded_repr(command),), False
    if type(command) not in (list, tuple):
        return (), False
    rendered = tuple(
        _bounded_repr(value)
        for value in command[:_MAX_PUBLICATION_COMMAND_FIELDS]
    )
    return rendered, len(command) > len(rendered)


def _bam_filter_cleanup_failures(primary):
    """Normalize trusted ReadIdStore cleanup facts for publication output."""
    attached = vars(primary).get("cleanup_failures", ())
    if type(attached) is not tuple:
        return ()

    failures = []
    allowed_operations = {
        "writer_stdin_close",
        "reader_stdout_close",
        "writer_wait",
        "reader_wait",
    }
    for attached_failure in attached:
        if type(attached_failure) is not util_misc.ReadIdCleanupFailure:
            continue
        values = vars(attached_failure)
        operation = values.get("operation")
        error_type = values.get("error_type")
        error_number = values.get("errno")
        if operation not in allowed_operations or type(error_type) is not str:
            continue
        failures.append(
            CleanupFailure(
                operation=operation,
                error_type=error_type,
                errno=error_number if type(error_number) is int else None,
                category="bam_filter_cleanup",
            )
        )
    return tuple(failures)


class LyraPublicationError(RuntimeError):
    """One structured publication failure retaining its original cause."""

    category = "publication"

    def __init__(self, stage, primary, cleanup_failures=()):
        self.stage = stage
        self.primary_type = type(primary).__name__
        self.primary_category = _known_primary_category(primary)

        primary_errno = None
        if isinstance(primary, OSError):
            candidate = getattr(primary, "errno", None)
            if type(candidate) is int:
                primary_errno = candidate
        self.primary_errno = primary_errno

        primary_returncode = None
        if isinstance(primary, subprocess.CalledProcessError):
            candidate = vars(primary).get("returncode")
            if type(candidate) is int:
                primary_returncode = candidate
        self.primary_returncode = primary_returncode
        self.primary_path = _known_primary_path(primary)
        (
            self.primary_command,
            self.primary_command_truncated,
        ) = _known_primary_command(primary)

        failures = tuple(cleanup_failures)
        self.cleanup_failures = failures[:_MAX_PUBLICATION_CLEANUP_FAILURES]
        self.cleanup_failures_truncated = len(failures) > len(
            self.cleanup_failures
        )

        details = [
            "category={}".format(_bounded_repr(self.category)),
            "stage={}".format(_bounded_repr(stage)),
            "primary_type={}".format(_bounded_repr(self.primary_type)),
        ]
        if self.primary_category is not None:
            details.append(
                "primary_category={}".format(
                    _bounded_repr(self.primary_category)
                )
            )
        if self.primary_errno is not None:
            details.append("primary_errno={}".format(self.primary_errno))
        if self.primary_returncode is not None:
            details.append(
                "primary_returncode={}".format(self.primary_returncode)
            )
        if self.primary_path is not None:
            details.append("primary_path={}".format(self.primary_path))
        if self.primary_command:
            details.append("primary_command={}".format(self.primary_command))
        details.append(
            "cleanup_failures={}".format(len(self.cleanup_failures))
        )
        if self.cleanup_failures_truncated:
            details.append("cleanup_failures_truncated=true")
        for index, failure in enumerate(self.cleanup_failures):
            failure_details = [
                "operation={}".format(_bounded_repr(failure.operation)),
                "role={}".format(_bounded_repr(failure.role)),
                "path={}".format(_bounded_repr(failure.path)),
                "error_type={}".format(_bounded_repr(failure.error_type)),
                "category={}".format(_bounded_repr(failure.category)),
            ]
            if failure.errno is not None:
                failure_details.append("errno={}".format(failure.errno))
            details.append(
                "cleanup[{}]=({})".format(
                    index,
                    ",".join(failure_details),
                )
            )
        super().__init__("Lyra publication error: " + "; ".join(details))


def _pairing_state(eligible_bam_count, bam_role_counts):
    """Return the locked pairing state and intrinsic BAM error category."""
    if eligible_bam_count > 2:
        return None, "bam_record_limit"
    if eligible_bam_count == 0:
        return None, None
    if eligible_bam_count == 1:
        if bam_role_counts == _SINGLE_END_ROLE_COUNTS:
            return PAIRING_SINGLE_END, None
        if bam_role_counts in (_PAIRED_R1_ROLE_COUNTS, _PAIRED_R2_ROLE_COUNTS):
            return PAIRING_INCOMPLETE, None
        return None, "bam_pairing_flags"
    if eligible_bam_count == 2:
        if bam_role_counts == _COMPLETE_PAIR_ROLE_COUNTS:
            return PAIRING_COMPLETE, None
        return None, "bam_pairing_flags"
    return None, "bam_pairing_flags"


@dataclass(frozen=True)
class LyraScoreRecord:
    """One validated native Lyra score row."""

    read_id: str
    score_text: str
    score: Decimal
    native_call: str
    line_number: int


@dataclass(frozen=True)
class LyraFragment:
    """One finalized Lyra fragment classification."""

    read_id: str
    n_scores: int
    pairing: str
    min_score: Decimal
    max_score: Decimal
    threshold: Decimal
    call: str


@dataclass(frozen=True)
class LyraReconciliationCounts:
    """Fixed aggregate counters for one Lyra reconciliation run."""

    input_bam_records: int
    eligible_bam_records: int
    score_records: int
    fragments: int
    single_end_fragments: int
    complete_pair_fragments: int
    incomplete_pair_fragments: int
    viral_fragment_calls: int
    nonviral_fragment_calls: int


def _canonical_output_decimal_length(value):
    """Predict canonical ordinary-decimal length without rendering it."""
    if value.is_zero():
        return 1

    value_tuple = value.as_tuple()
    digit_count = len(value_tuple.digits)
    exponent = value_tuple.exponent
    for digit in reversed(value_tuple.digits):
        if digit != 0:
            break
        digit_count -= 1
        exponent += 1

    if exponent >= 0:
        return digit_count + exponent
    if digit_count + exponent > 0:
        return digit_count + 1
    return 2 - exponent


def validate_lyra_threshold(threshold):
    """Return a finite inclusive ``[0, 1]`` threshold as a Decimal."""
    if type(threshold) is bool:
        raise LyraInputError(
            category="threshold",
            field="threshold",
            reason="threshold must not be a boolean",
            offending_value=threshold,
        )
    try:
        value = threshold if isinstance(threshold, Decimal) else Decimal(str(threshold))
    except (InvalidOperation, TypeError, ValueError) as exc:
        raise LyraInputError(
            category="threshold",
            field="threshold",
            reason="threshold must be a decimal value",
            offending_value=threshold,
        ) from exc
    if not value.is_finite() or not _MIN_SCORE <= value <= _MAX_SCORE:
        raise LyraInputError(
            category="threshold",
            field="threshold",
            reason="threshold must be finite and between zero and one inclusive",
            offending_value=threshold,
        )
    if _canonical_output_decimal_length(value) > _MAX_CANONICAL_THRESHOLD_LENGTH:
        raise LyraInputError(
            category="threshold",
            field="threshold",
            reason="canonical threshold text exceeds 160 characters",
            offending_value=threshold,
        )
    return value


def _read_id_key(read_id):
    return read_id.encode("utf-8", errors="strict")


def _canonical_decimal_text(text):
    """Remove redundant fractional zeroes from validated score text."""
    if "." not in text:
        return text
    whole, fractional = text.split(".", 1)
    fractional = fractional.rstrip("0")
    if fractional:
        return whole + "." + fractional
    return whole or "0"


def _canonical_output_decimal(value):
    """Render a Decimal as canonical ordinary text without arithmetic."""
    if value.is_zero():
        return "0"
    return _canonical_decimal_text(format(value, "f"))


def _write_tsv_row(binary_stream, values):
    """Write one strict UTF-8 TSV row with an LF terminator."""
    row = "\t".join(str(value) for value in values) + "\n"
    binary_stream.write(row.encode("utf-8", errors="strict"))


@contextmanager
def _duplicate_stage_raw_stream(stage, mode, truncate=False):
    """Yield a stream around a duplicate while retaining the owner handle."""
    _assert_stage_handle(stage, require_named_entry=False)
    if truncate:
        os.ftruncate(stage.descriptor, 0)
    os.lseek(stage.descriptor, 0, os.SEEK_SET)
    duplicate = os.dup(stage.descriptor)
    try:
        with os.fdopen(duplicate, mode) as stream:
            duplicate = None
            yield stream
    finally:
        if duplicate is not None:
            os.close(duplicate)
    _assert_stage_handle(stage, require_named_entry=False)


@contextmanager
def _stage_output_stream(stage):
    with _duplicate_stage_raw_stream(stage, "wb", truncate=True) as raw_stream:
        suffix = stage.destination.suffix
        if suffix == ".tsv.gz":
            with gzip.GzipFile(fileobj=raw_stream, mode="wb") as output_stream:
                yield output_stream
        elif suffix == ".tsv.zst":
            compressor = zstd.ZstdCompressor(
                level=10,
                threads=util_misc.sanitize_thread_count(None),
            )
            with compressor.stream_writer(raw_stream) as output_stream:
                yield output_stream
        else:
            yield raw_stream


@contextmanager
def _stage_input_stream(stage):
    with _duplicate_stage_raw_stream(stage, "rb") as raw_stream:
        suffix = stage.destination.suffix
        if suffix == ".tsv.gz":
            with gzip.GzipFile(fileobj=raw_stream, mode="rb") as input_stream:
                yield input_stream
        elif suffix == ".tsv.zst":
            decompressor = zstd.ZstdDecompressor()
            with decompressor.stream_reader(raw_stream) as input_stream:
                yield input_stream
        else:
            yield raw_stream


@contextmanager
def _normalized_output_stream(normalized_output):
    if isinstance(normalized_output, ArtifactStage):
        with _stage_output_stream(normalized_output) as output_stream:
            yield output_stream
        return
    path = os.fspath(normalized_output)
    with util_file.open_or_gzopen(path, "wb") as output_stream:
        yield output_stream


def _write_normalized(store, normalized_output):
    """Stream one exact normalized row per finalized fragment."""
    threshold = store.threshold
    threshold_text = _canonical_output_decimal(threshold)
    pairing_values = {
        PAIRING_SINGLE_END,
        PAIRING_COMPLETE,
        PAIRING_INCOMPLETE,
    }
    call_values = {CALL_VIRAL, CALL_NONVIRAL}
    row_count = 0
    with _normalized_output_stream(normalized_output) as output_stream:
        _write_tsv_row(output_stream, NORMALIZED_HEADER)
        for fragment in store.iter_fragments():
            if fragment.pairing not in pairing_values:
                raise LyraArtifactConsistencyError(
                    category="pairing",
                    field="LYRA_PAIRING",
                    expected="locked pairing vocabulary",
                    actual=fragment.pairing,
                    read_id=fragment.read_id,
                )
            if fragment.call not in call_values:
                raise LyraArtifactConsistencyError(
                    category="call",
                    field="LYRA_CALL",
                    expected="locked call vocabulary",
                    actual=fragment.call,
                    read_id=fragment.read_id,
                )
            if fragment.threshold != threshold:
                raise LyraArtifactConsistencyError(
                    category="threshold",
                    field="LYRA_THRESHOLD",
                    expected=threshold,
                    actual=fragment.threshold,
                    read_id=fragment.read_id,
                )
            _write_tsv_row(
                output_stream,
                (
                    store.sample_id,
                    fragment.read_id,
                    fragment.n_scores,
                    fragment.pairing,
                    _canonical_output_decimal(fragment.min_score),
                    _canonical_output_decimal(fragment.max_score),
                    threshold_text,
                    fragment.call,
                ),
            )
            row_count += 1
    return row_count


def _validate_artifact_counts(counts, normalized_rows, output_bam_records):
    """Require finalized and in-pass artifact counts to agree exactly."""
    if normalized_rows != counts.fragments:
        raise LyraArtifactConsistencyError(
            category="normalized_row_count",
            field="LYRA_FRAGMENTS",
            expected=counts.fragments,
            actual=normalized_rows,
        )
    if counts.eligible_bam_records != counts.score_records:
        raise LyraArtifactConsistencyError(
            category="eligible_score_record_count",
            field="LYRA_ELIGIBLE_BAM_RECORDS",
            expected=counts.score_records,
            actual=counts.eligible_bam_records,
        )

    pairing_score_records = (
        counts.single_end_fragments
        + counts.incomplete_pair_fragments
        + 2 * counts.complete_pair_fragments
    )
    if counts.score_records != pairing_score_records:
        raise LyraArtifactConsistencyError(
            category="score_pairing_record_count",
            field="LYRA_SCORE_RECORDS",
            expected=pairing_score_records,
            actual=counts.score_records,
        )

    pairing_fragments = (
        counts.single_end_fragments
        + counts.complete_pair_fragments
        + counts.incomplete_pair_fragments
    )
    if counts.fragments != pairing_fragments:
        raise LyraArtifactConsistencyError(
            category="fragment_pairing_count",
            field="LYRA_FRAGMENTS",
            expected=pairing_fragments,
            actual=counts.fragments,
        )

    call_fragments = (
        counts.viral_fragment_calls + counts.nonviral_fragment_calls
    )
    if counts.fragments != call_fragments:
        raise LyraArtifactConsistencyError(
            category="fragment_call_count",
            field="LYRA_FRAGMENTS",
            expected=call_fragments,
            actual=counts.fragments,
        )
    if not 0 <= output_bam_records <= counts.input_bam_records:
        raise LyraArtifactConsistencyError(
            category="output_bam_record_range",
            field="LYRA_OUTPUT_BAM_RECORDS",
            expected="0..{}".format(counts.input_bam_records),
            actual=output_bam_records,
        )
    if (output_bam_records == 0) != (counts.viral_fragment_calls == 0):
        raise LyraArtifactConsistencyError(
            category="output_bam_viral_equivalence",
            field="LYRA_OUTPUT_BAM_RECORDS",
            expected=(counts.viral_fragment_calls == 0),
            actual=(output_bam_records == 0),
        )


def _write_summary(store, summary_output, output_bam_records):
    """Write one exact summary row from finalized and checked counts."""
    counts = store.counts
    if isinstance(summary_output, ArtifactStage):
        output_context = _stage_output_stream(summary_output)
    else:
        output_context = open(os.fspath(summary_output), "wb")
    with output_context as output_stream:
        _write_tsv_row(output_stream, SUMMARY_HEADER)
        _write_tsv_row(
            output_stream,
            (
                store.sample_id,
                _canonical_output_decimal(store.threshold),
                counts.input_bam_records,
                counts.eligible_bam_records,
                counts.score_records,
                counts.fragments,
                counts.single_end_fragments,
                counts.complete_pair_fragments,
                counts.incomplete_pair_fragments,
                counts.viral_fragment_calls,
                counts.nonviral_fragment_calls,
                output_bam_records,
            ),
        )


def _write_viral_bam(
    store,
    viral_bam_output,
    work_dir=None,
):
    """Stream exact Viral QNAMEs into a faithful source-order BAM."""
    expected_id_count = store.counts.viral_fragment_calls
    output_stage = (
        viral_bam_output
        if isinstance(viral_bam_output, ArtifactStage)
        else None
    )
    if output_stage is None:
        output_path = os.fspath(viral_bam_output)
        output_descriptor = None
    else:
        _assert_stage_handle(output_stage, require_named_entry=False)
        os.ftruncate(output_stage.descriptor, 0)
        os.lseek(output_stage.descriptor, 0, os.SEEK_SET)
        output_path = output_stage.display_path
        output_descriptor = output_stage.descriptor
    with tempfile.TemporaryDirectory(
        prefix="lyra_viral_bam_",
        dir=work_dir,
    ) as temporary_directory:
        database_path = os.path.join(
            temporary_directory,
            "viral_read_ids.sqlite3",
        )
        with util_misc.ReadIdStore(database_path) as read_ids:
            actual_id_count = read_ids.extend(store.iter_viral_read_ids())
            if actual_id_count != expected_id_count:
                raise LyraArtifactConsistencyError(
                    category="viral_read_id_count",
                    field="LYRA_VIRAL_FRAGMENT_CALLS",
                    expected=expected_id_count,
                    actual=actual_id_count,
                )
            output_bam_records = read_ids.filter_bam_by_ids(
                store.source_bam_display_path,
                output_path,
                include=True,
                in_bam_fd=store.source_bam_fd,
                out_bam_fd=output_descriptor,
            )
            _assert_source_bam_identity(store, "after_bam_filter")
            if output_stage is not None:
                _assert_stage_handle(
                    output_stage,
                    require_named_entry=False,
                )
            return output_bam_records


def _assert_parent_anchor(anchor, stage):
    """Require a retained parent descriptor and current path to stay planned."""
    expected = anchor.identity
    try:
        retained_stat = os.fstat(anchor.descriptor)
    except OSError as exc:
        raise LyraPathError(
            category="output_parent_identity",
            role="output_parent",
            path=anchor.parent_path,
            reason="retained output parent descriptor is unreadable",
            stage=stage,
        ) from exc
    retained_identity = _path_identity(retained_stat)
    if retained_identity != expected or not stat.S_ISDIR(retained_stat.st_mode):
        raise LyraPathError(
            category="output_parent_identity",
            role="output_parent",
            path=anchor.parent_path,
            reason="retained output parent descriptor does not match planned identity",
            stage=stage,
        )

    try:
        current_stat = os.stat(anchor.parent_path)
    except OSError as exc:
        raise LyraPathError(
            category="output_parent_identity",
            role="output_parent",
            path=anchor.parent_path,
            reason="planned output parent path is missing or unreadable",
            stage=stage,
        ) from exc
    if (
        _path_identity(current_stat) != expected
        or not stat.S_ISDIR(current_stat.st_mode)
    ):
        raise LyraPathError(
            category="output_parent_identity",
            role="output_parent",
            path=anchor.parent_path,
            reason="planned output parent path no longer matches retained identity",
            stage=stage,
        )


def _assert_stage_handle(stage, require_named_entry=True):
    """Require a stage descriptor and optional directory entry to stay owned."""
    expected = stage.object_identity
    try:
        retained_stat = os.fstat(stage.descriptor)
    except OSError as exc:
        raise LyraArtifactConsistencyError(
            category="stage_identity",
            field=stage.role,
            expected=expected,
            actual="descriptor_unreadable",
        ) from exc
    retained_identity = _path_identity(retained_stat)
    if retained_identity != expected or not stat.S_ISREG(retained_stat.st_mode):
        raise LyraArtifactConsistencyError(
            category="stage_identity",
            field=stage.role,
            expected=expected,
            actual=retained_identity,
        )
    if not require_named_entry:
        return

    try:
        entry_stat = os.stat(
            stage.basename,
            dir_fd=stage.parent.descriptor,
            follow_symlinks=False,
        )
    except FileNotFoundError as exc:
        raise LyraArtifactConsistencyError(
            category="stage_identity",
            field=stage.role,
            expected=expected,
            actual="missing",
        ) from exc
    except OSError as exc:
        raise LyraArtifactConsistencyError(
            category="stage_identity",
            field=stage.role,
            expected=expected,
            actual="entry_unreadable",
        ) from exc
    entry_identity = _path_identity(entry_stat)
    if entry_identity != expected or not stat.S_ISREG(entry_stat.st_mode):
        raise LyraArtifactConsistencyError(
            category="stage_identity",
            field=stage.role,
            expected=expected,
            actual=entry_identity,
        )


def _open_parent_anchor(destination):
    """Open and verify one immutable planned output parent."""
    expected = destination.destination_key[:2]
    descriptor = os.open(
        destination.parent_path,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
    )
    anchor = ParentDirectoryAnchor(
        parent_path=destination.parent_path,
        identity=expected,
        descriptor=descriptor,
    )
    try:
        _assert_parent_anchor(anchor, "create_stages")
    except BaseException:
        os.close(descriptor)
        raise
    return anchor


def _create_artifact_stage(destination, parent):
    """Create and retain one exclusive stage relative to its parent handle."""
    flags = (
        os.O_RDWR
        | os.O_CREAT
        | os.O_EXCL
        | os.O_NOFOLLOW
        | os.O_CLOEXEC
    )
    last_collision = None
    for _ in range(100):
        basename = ".lyra-{}-{}{}".format(
            destination.role,
            secrets.token_hex(16),
            destination.suffix,
        )
        try:
            descriptor = os.open(
                basename,
                flags,
                0o600,
                dir_fd=parent.descriptor,
            )
        except FileExistsError as exc:
            last_collision = exc
            continue

        try:
            stage_stat = os.fstat(descriptor)
            artifact_stage = ArtifactStage(
                role=destination.role,
                basename=basename,
                display_path=os.path.join(parent.parent_path, basename),
                descriptor=descriptor,
                object_identity=_path_identity(stage_stat),
                destination=destination,
                parent=parent,
            )
            _assert_stage_handle(artifact_stage)
            return artifact_stage
        except BaseException:
            try:
                os.unlink(basename, dir_fd=parent.descriptor)
            finally:
                os.close(descriptor)
            raise
    raise last_collision


def _readback_error(category, field=None, expected=None, actual=None, read_id=None):
    return LyraArtifactConsistencyError(
        category=category,
        field=field,
        expected=expected,
        actual=actual,
        read_id=read_id,
    )


def _stage_suffix(stage_or_path):
    if isinstance(stage_or_path, ArtifactStage):
        return stage_or_path.destination.suffix
    return os.fspath(stage_or_path)


@contextmanager
def _staged_input_stream(stage_or_path):
    if isinstance(stage_or_path, ArtifactStage):
        with _stage_input_stream(stage_or_path) as input_stream:
            yield input_stream
        return
    path = os.fspath(stage_or_path)
    with util_file.open_or_gzopen(path, "rb") as input_stream:
        yield input_stream


def _iter_strict_staged_rows(stage_or_path, category):
    """Stream strict UTF-8, LF-terminated physical rows from one stage."""
    try:
        with _staged_input_stream(stage_or_path) as raw_stream:
            if isinstance(raw_stream, io.BufferedIOBase):
                buffered_stream = raw_stream
            else:
                buffered_stream = io.BufferedReader(raw_stream)
            try:
                for line_number, raw_line in enumerate(buffered_stream, start=1):
                    if not raw_line.endswith(b"\n"):
                        raise _readback_error(
                            category,
                            field="line_ending",
                            expected="LF-terminated row",
                            actual="missing final LF at line {}".format(line_number),
                        )
                    content = raw_line[:-1]
                    if b"\r" in content:
                        raise _readback_error(
                            category,
                            field="line_ending",
                            expected="LF only",
                            actual="carriage return at line {}".format(line_number),
                        )
                    try:
                        text = content.decode("utf-8", errors="strict")
                    except UnicodeDecodeError as exc:
                        raise _readback_error(
                            category,
                            field="utf8",
                            expected="strict UTF-8",
                            actual="invalid UTF-8 at line {}".format(line_number),
                        ) from exc
                    yield line_number, tuple(text.split("\t"))
            finally:
                if buffered_stream is not raw_stream:
                    buffered_stream.close()
        if _stage_suffix(stage_or_path).endswith(".tsv.zst"):
            if isinstance(stage_or_path, ArtifactStage):
                with _duplicate_stage_raw_stream(
                    stage_or_path,
                    "rb",
                ) as compressed_stream:
                    _validate_complete_zstd_frames_stream(compressed_stream)
            else:
                _validate_complete_zstd_frames(os.fspath(stage_or_path))
    except LyraArtifactConsistencyError:
        raise
    except _COMPRESSION_EXCEPTIONS as exc:
        raise _readback_error(
            category,
            field="compression",
            expected="complete compressed stream",
            actual="decoder rejected staged output",
        ) from exc


def _validate_staged_normalized(store, stage_or_path):
    """Stream-compare a closed normalized stage with finalized fragments."""
    rows = _iter_strict_staged_rows(stage_or_path, "normalized_readback")
    try:
        _, header = next(rows)
    except StopIteration as exc:
        raise _readback_error(
            "normalized_readback",
            field="header",
            expected=NORMALIZED_HEADER,
            actual="missing",
        ) from exc
    if header != NORMALIZED_HEADER:
        raise _readback_error(
            "normalized_readback",
            field="header",
            expected=NORMALIZED_HEADER,
            actual=header,
        )

    pairing_values = {
        PAIRING_SINGLE_END,
        PAIRING_COMPLETE,
        PAIRING_INCOMPLETE,
    }
    call_values = {CALL_VIRAL, CALL_NONVIRAL}
    threshold_text = _canonical_output_decimal(store.threshold)
    fragments = store.iter_fragments()
    row_count = 0
    try:
        for line_number, fields in rows:
            if len(fields) != len(NORMALIZED_HEADER):
                raise _readback_error(
                    "normalized_readback",
                    field="row_width",
                    expected=len(NORMALIZED_HEADER),
                    actual=len(fields),
                )
            try:
                fragment = next(fragments)
            except StopIteration as exc:
                raise _readback_error(
                    "normalized_readback",
                    field="row_count",
                    expected=store.counts.fragments,
                    actual="extra row at line {}".format(line_number),
                ) from exc
            if fields[3] not in pairing_values:
                raise _readback_error(
                    "normalized_readback",
                    field="LYRA_PAIRING",
                    expected="locked pairing vocabulary",
                    actual=fields[3],
                    read_id=fragment.read_id,
                )
            if fields[7] not in call_values:
                raise _readback_error(
                    "normalized_readback",
                    field="LYRA_CALL",
                    expected="locked call vocabulary",
                    actual=fields[7],
                    read_id=fragment.read_id,
                )
            expected_fields = (
                store.sample_id,
                fragment.read_id,
                str(fragment.n_scores),
                fragment.pairing,
                _canonical_output_decimal(fragment.min_score),
                _canonical_output_decimal(fragment.max_score),
                threshold_text,
                fragment.call,
            )
            for field_name, expected, actual in zip(
                NORMALIZED_HEADER,
                expected_fields,
                fields,
            ):
                if actual != expected:
                    raise _readback_error(
                        "normalized_readback",
                        field=field_name,
                        expected=expected,
                        actual=actual,
                        read_id=fragment.read_id,
                    )
            row_count += 1

        try:
            missing_fragment = next(fragments)
        except StopIteration:
            missing_fragment = None
        if missing_fragment is not None:
            raise _readback_error(
                "normalized_readback",
                field="row_count",
                expected=store.counts.fragments,
                actual=row_count,
                read_id=missing_fragment.read_id,
            )
    finally:
        close = getattr(fragments, "close", None)
        if close is not None:
            close()

    if row_count != store.counts.fragments:
        raise _readback_error(
            "normalized_readback",
            field="row_count",
            expected=store.counts.fragments,
            actual=row_count,
        )
    return row_count


def _validate_staged_bam(store, stage_or_path):
    """Validate the staged BAM header and stream its actual record count."""
    try:
        if isinstance(stage_or_path, ArtifactStage):
            bam_context = _duplicate_stage_raw_stream(stage_or_path, "rb")
        else:
            bam_context = open(os.fspath(stage_or_path), "rb")
        with bam_context as bam_stream:
            with pysam.AlignmentFile(
                bam_stream,
                "rb",
                check_sq=False,
            ) as bam:
                actual_header = bam.header.to_dict()
                if actual_header != store.source_bam_header:
                    raise _readback_error(
                        "bam_readback",
                        field="header",
                        expected=store.source_bam_header,
                        actual=actual_header,
                    )
                record_count = 0
                for _ in bam.fetch(until_eof=True):
                    record_count += 1
    except LyraArtifactConsistencyError:
        raise
    except (OSError, ValueError) as exc:
        raise _readback_error(
            "bam_readback",
            field="bam",
            expected="readable complete BAM",
            actual=type(exc).__name__,
        ) from exc
    return record_count


def _validate_staged_summary(
    store,
    stage_or_path,
    normalized_rows,
    bam_records,
):
    """Validate the exact one-row summary and rerun all count equations."""
    rows = _iter_strict_staged_rows(stage_or_path, "summary_readback")
    try:
        _, header = next(rows)
        _, fields = next(rows)
    except StopIteration as exc:
        raise _readback_error(
            "summary_readback",
            field="rows",
            expected="header and one data row",
            actual="missing row",
        ) from exc
    if header != SUMMARY_HEADER:
        raise _readback_error(
            "summary_readback",
            field="header",
            expected=SUMMARY_HEADER,
            actual=header,
        )
    if len(fields) != len(SUMMARY_HEADER):
        raise _readback_error(
            "summary_readback",
            field="row_width",
            expected=len(SUMMARY_HEADER),
            actual=len(fields),
        )
    try:
        _, extra = next(rows)
    except StopIteration:
        extra = None
    if extra is not None:
        raise _readback_error(
            "summary_readback",
            field="rows",
            expected="exactly one data row",
            actual="extra row",
        )

    counts = store.counts
    expected_fields = (
        store.sample_id,
        _canonical_output_decimal(store.threshold),
        str(counts.input_bam_records),
        str(counts.eligible_bam_records),
        str(counts.score_records),
        str(counts.fragments),
        str(counts.single_end_fragments),
        str(counts.complete_pair_fragments),
        str(counts.incomplete_pair_fragments),
        str(counts.viral_fragment_calls),
        str(counts.nonviral_fragment_calls),
        str(bam_records),
    )
    for field_name, expected, actual in zip(
        SUMMARY_HEADER,
        expected_fields,
        fields,
    ):
        if actual != expected:
            raise _readback_error(
                "summary_readback",
                field=field_name,
                expected=expected,
                actual=actual,
            )
    try:
        _validate_artifact_counts(counts, normalized_rows, bam_records)
    except LyraArtifactConsistencyError as exc:
        raise _readback_error(
            "summary_readback",
            field="equations",
            expected="all artifact count equations",
            actual=exc.category,
        ) from exc


def _validate_staged_artifacts(
    store,
    normalized_stage,
    summary_stage,
    bam_stage,
    producer_normalized_rows,
    producer_bam_records,
    stage_callback=None,
    validated_stage_identities=None,
):
    if stage_callback is not None:
        stage_callback("validate_normalized_readback")
    normalized_rows = _validate_staged_normalized(
        store,
        normalized_stage,
    )
    if stage_callback is not None:
        stage_callback("validate_bam_readback")
    bam_records = _validate_staged_bam(store, bam_stage)
    if normalized_rows != producer_normalized_rows:
        raise _readback_error(
            "normalized_readback",
            field="producer_row_count",
            expected=producer_normalized_rows,
            actual=normalized_rows,
        )
    if bam_records != producer_bam_records:
        raise _readback_error(
            "bam_readback",
            field="producer_record_count",
            expected=producer_bam_records,
            actual=bam_records,
        )
    if stage_callback is not None:
        stage_callback("validate_summary_readback")
    _validate_staged_summary(
        store,
        summary_stage,
        normalized_rows,
        bam_records,
    )
    if stage_callback is not None:
        stage_callback("validate_stage_identities")
    identities = {}
    for artifact_stage in (
        normalized_stage,
        summary_stage,
        bam_stage,
    ):
        _assert_parent_anchor(artifact_stage.parent, "validate_stage_identities")
        _assert_stage_handle(artifact_stage)
        identities[artifact_stage.role] = _file_identity(
            os.fstat(artifact_stage.descriptor)
        )
    if validated_stage_identities is not None:
        validated_stage_identities.clear()
        validated_stage_identities.update(identities)
    return normalized_rows, bam_records


def _validate_artifact_output_suffixes(
    normalized_output,
    summary_output,
    viral_bam_output,
):
    """Return output paths after exact case-sensitive suffix validation."""
    caller_paths = (
        normalized_output,
        summary_output,
        viral_bam_output,
    )
    converted_paths = tuple(os.fspath(path) for path in caller_paths)
    checks = (
        (
            "normalized_output",
            converted_paths[0],
            _SUPPORTED_SCORE_SUFFIXES,
            "normalized output must end with .tsv, .tsv.gz, or .tsv.zst",
        ),
        (
            "summary_output",
            converted_paths[1],
            (".tsv",),
            "summary output must end with .tsv",
        ),
        (
            "viral_bam_output",
            converted_paths[2],
            (".bam",),
            "viral BAM output must end with .bam",
        ),
    )
    for index, (field, path, suffixes, reason) in enumerate(checks):
        if not isinstance(path, str) or not path.endswith(suffixes):
            raise LyraInputError(
                category="output_extension",
                path=caller_paths[index],
                field=field,
                reason=reason,
                offending_value=path,
            )
    return converted_paths


def _artifact_suffix(role, path, suffixes, reason):
    for suffix in sorted(suffixes, key=len, reverse=True):
        if isinstance(path, str) and path.endswith(suffix):
            return suffix
    raise LyraInputError(
        category="output_extension",
        path=path,
        field=role + "_output",
        reason=reason,
        offending_value=path,
    )


def _resolved_artifact_destination(role, caller_path, suffixes, reason):
    suffix = _artifact_suffix(role, caller_path, suffixes, reason)
    absolute_path = caller_path
    if not os.path.isabs(absolute_path):
        absolute_path = os.path.join(os.getcwd(), absolute_path)
    unresolved_parent, basename = os.path.split(absolute_path)
    try:
        parent_path = os.path.realpath(unresolved_parent, strict=True)
        parent_stat = os.stat(parent_path)
    except OSError as exc:
        raise LyraPathError(
            category="output_parent",
            role=role,
            path=caller_path,
            reason="output parent must be an existing directory",
        ) from exc
    if not stat.S_ISDIR(parent_stat.st_mode):
        raise LyraPathError(
            category="output_parent",
            role=role,
            path=caller_path,
            reason="output parent must be an existing directory",
        )
    final_path = os.path.join(parent_path, basename)
    return ArtifactDestination(
        role=role,
        caller_path=caller_path,
        final_path=final_path,
        parent_path=parent_path,
        destination_key=(parent_stat.st_dev, parent_stat.st_ino, basename),
        suffix=suffix,
    )


def _path_identity(stat_result):
    return (stat_result.st_dev, stat_result.st_ino)


def _build_artifact_path_plan(
    score_path,
    bam_path,
    normalized_output,
    summary_output,
    viral_bam_output,
):
    """Convert paths once and reject every deterministic destination hazard."""
    converted_paths = tuple(
        os.fspath(path)
        for path in (
            score_path,
            bam_path,
            normalized_output,
            summary_output,
            viral_bam_output,
        )
    )
    converted_score, converted_bam = converted_paths[:2]
    if not isinstance(converted_score, str) or not converted_score.endswith(
        _SUPPORTED_SCORE_SUFFIXES
    ):
        raise LyraInputError(
            category="extension",
            path=converted_score,
            field="file",
            reason="score path must end with .tsv, .tsv.gz, or .tsv.zst",
            offending_value=converted_score,
        )

    destinations = (
        _resolved_artifact_destination(
            "normalized",
            converted_paths[2],
            _SUPPORTED_SCORE_SUFFIXES,
            "normalized output must end with .tsv, .tsv.gz, or .tsv.zst",
        ),
        _resolved_artifact_destination(
            "summary",
            converted_paths[3],
            (".tsv",),
            "summary output must end with .tsv",
        ),
        _resolved_artifact_destination(
            "viral_bam",
            converted_paths[4],
            (".bam",),
            "viral BAM output must end with .bam",
        ),
    )

    destinations_by_key = {}
    for destination in destinations:
        conflict = destinations_by_key.get(destination.destination_key)
        if conflict is not None:
            raise LyraPathError(
                category="output_output_alias",
                role=destination.role,
                path=destination.caller_path,
                conflicting_role=conflict.role,
                conflicting_path=conflict.caller_path,
                reason="artifact outputs must name distinct directory entries",
            )
        destinations_by_key[destination.destination_key] = destination

    input_paths = (
        ("score", converted_score),
        ("bam", converted_bam),
    )
    input_identities = tuple(
        (role, path, _path_identity(os.stat(path)))
        for role, path in input_paths
    )

    existing_entries = []
    for destination in destinations:
        try:
            entry_stat = os.lstat(destination.final_path)
        except FileNotFoundError:
            continue
        except OSError as exc:
            raise LyraPathError(
                category="output_exists",
                role=destination.role,
                path=destination.caller_path,
                reason="artifact output directory entry is not safely absent",
            ) from exc

        identity = None
        try:
            identity = _path_identity(os.stat(destination.final_path))
        except OSError:
            pass
        existing_entries.append((destination, entry_stat, identity))

    for destination, _, identity in existing_entries:
        if identity is None:
            continue
        for input_role, input_path, input_identity in input_identities:
            if identity == input_identity:
                raise LyraPathError(
                    category="input_output_alias",
                    role=destination.role,
                    path=destination.caller_path,
                    conflicting_role=input_role,
                    conflicting_path=input_path,
                    reason="artifact output must not identify an input file",
                )

    for index, (destination, _, identity) in enumerate(existing_entries):
        if identity is None:
            continue
        for conflict, _, conflict_identity in existing_entries[:index]:
            if conflict_identity is not None and identity == conflict_identity:
                raise LyraPathError(
                    category="output_output_alias",
                    role=destination.role,
                    path=destination.caller_path,
                    conflicting_role=conflict.role,
                    conflicting_path=conflict.caller_path,
                    reason="artifact outputs must not identify the same entry",
                )

    if existing_entries:
        destination = existing_entries[0][0]
        raise LyraPathError(
            category="output_exists",
            role=destination.role,
            path=destination.caller_path,
            reason="artifact output directory entry already exists",
        )

    return LyraArtifactPathPlan(
        score_path=converted_score,
        bam_path=converted_bam,
        normalized=destinations[0],
        summary=destinations[1],
        viral_bam=destinations[2],
    )


def postprocess_lyra(
    score_path,
    bam_path,
    sample_id,
    threshold,
    normalized_output,
    summary_output,
    viral_bam_output,
    work_dir=None,
):
    """Preflight paths before reconciling and generating Lyra artifacts."""
    validated_sample_id = validate_sample_id(sample_id)
    validated_threshold = validate_lyra_threshold(threshold)
    path_plan = _build_artifact_path_plan(
        score_path,
        bam_path,
        normalized_output,
        summary_output,
        viral_bam_output,
    )
    with reconcile_lyra_fragments(
        path_plan.score_path,
        path_plan.bam_path,
        validated_sample_id,
        validated_threshold,
        work_dir=work_dir,
    ) as store:
        write_lyra_artifacts(
            store,
            path_plan,
            work_dir=work_dir,
        )


def _assert_path_plan_available(path_plan, stage):
    """Recheck resolved final names without rebuilding their immutable plan."""
    if not isinstance(path_plan, LyraArtifactPathPlan):
        raise TypeError("path_plan must be a LyraArtifactPathPlan")
    for destination in (
        path_plan.normalized,
        path_plan.summary,
        path_plan.viral_bam,
    ):
        try:
            os.lstat(destination.final_path)
        except FileNotFoundError:
            continue
        except OSError as exc:
            raise LyraPathError(
                category="output_exists",
                role=destination.role,
                path=destination.caller_path,
                reason="artifact output directory entry is not safely absent",
                stage=stage,
            ) from exc
        raise LyraPathError(
            category="output_exists",
            role=destination.role,
            path=destination.caller_path,
            reason="artifact output directory entry already exists",
            stage=stage,
        )


def _link_no_clobber(stage_path, final_path):
    os.link(stage_path, final_path, follow_symlinks=False)


def _unlink_path(path):
    os.unlink(path)


def _fsync_file(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _cleanup_failure(operation, role, path, status, error):
    return CleanupOutcome(
        operation=operation,
        role=role,
        path=path,
        status=status,
        error_type=type(error).__name__,
        errno=getattr(error, "errno", None),
    )


class LyraArtifactTransaction:
    """Own staged generation, durable publication, and bounded rollback."""

    def __init__(self, store, path_plan, work_dir=None):
        if not isinstance(path_plan, LyraArtifactPathPlan):
            raise TypeError("path_plan must be a LyraArtifactPathPlan")
        self.store = store
        self.path_plan = path_plan
        self.work_dir = work_dir
        self.stage = "initialized"
        self._parent_anchors = []
        self._stages = []
        self._published = []
        self._pending_publication = None
        self._cleanup_outcomes = ()
        self._owned_descriptors_closed = False
        self._validated_stage_identities = {}

    @property
    def stages(self):
        return tuple(self._stages)

    @property
    def published(self):
        return tuple(self._published)

    @property
    def cleanup_outcomes(self):
        return self._cleanup_outcomes

    def _set_stage(self, stage):
        self.stage = stage

    def _create_stages(self):
        self.stage = "create_stages"
        destinations = (
            self.path_plan.normalized,
            self.path_plan.summary,
            self.path_plan.viral_bam,
        )
        anchors_by_identity = {}
        try:
            for destination in destinations:
                parent_identity = destination.destination_key[:2]
                anchor = anchors_by_identity.get(parent_identity)
                if anchor is None:
                    anchor = _open_parent_anchor(destination)
                    self._parent_anchors.append(anchor)
                    anchors_by_identity[parent_identity] = anchor
                else:
                    _assert_parent_anchor(anchor, self.stage)

            for destination in destinations:
                parent = anchors_by_identity[destination.destination_key[:2]]
                self._stages.append(
                    _create_artifact_stage(destination, parent)
                )
        except BaseException:
            self._cleanup_partial_acquisition()
            raise

    def _cleanup_partial_acquisition(self):
        for artifact_stage in reversed(self._stages):
            try:
                entry_stat = os.stat(
                    artifact_stage.basename,
                    dir_fd=artifact_stage.parent.descriptor,
                    follow_symlinks=False,
                )
                if _path_identity(entry_stat) == artifact_stage.object_identity:
                    os.unlink(
                        artifact_stage.basename,
                        dir_fd=artifact_stage.parent.descriptor,
                    )
            except BaseException:
                pass
        self._close_owned_descriptors()

    def _close_owned_descriptors(self):
        if self._owned_descriptors_closed:
            return ()
        self._owned_descriptors_closed = True
        failures = []
        for artifact_stage in reversed(self._stages):
            try:
                os.close(artifact_stage.descriptor)
            except BaseException as error:
                failures.append(
                    _cleanup_failure(
                        "close_stage_descriptor",
                        artifact_stage.role,
                        artifact_stage.display_path,
                        "close_failed",
                        error,
                    )
                )
        for anchor in reversed(self._parent_anchors):
            try:
                os.close(anchor.descriptor)
            except BaseException as error:
                failures.append(
                    _cleanup_failure(
                        "close_parent_descriptor",
                        "output_parent",
                        anchor.parent_path,
                        "close_failed",
                        error,
                    )
                )
        return tuple(failures)

    def _stages_in_publication_order(self):
        by_role = {stage.role: stage for stage in self._stages}
        return (
            by_role["normalized"],
            by_role["viral_bam"],
            by_role["summary"],
        )

    def _generate_and_validate(self):
        normalized_stage, bam_stage, summary_stage = (
            self._stages_in_publication_order()
        )
        self.stage = "generate_normalized"
        normalized_rows = _write_normalized(
            self.store,
            normalized_stage,
        )
        self.stage = "generate_viral_bam"
        output_bam_records = _write_viral_bam(
            self.store,
            bam_stage,
            work_dir=self.work_dir,
        )
        self.stage = "validate_generation_counts"
        _validate_artifact_counts(
            self.store.counts,
            normalized_rows,
            output_bam_records,
        )
        self.stage = "generate_summary"
        _write_summary(
            self.store,
            summary_stage,
            output_bam_records,
        )
        _validate_staged_artifacts(
            self.store,
            normalized_stage,
            summary_stage,
            bam_stage,
            normalized_rows,
            output_bam_records,
            stage_callback=self._set_stage,
            validated_stage_identities=self._validated_stage_identities,
        )

    def _sync_stages(self):
        publication_stages = self._stages_in_publication_order()
        for stage in publication_stages:
            self.stage = "fsync_{}_stage".format(stage.role)
            _fsync_file(stage.stage_path)
        for stage in publication_stages:
            self.stage = "fsync_{}_stage_parent".format(stage.role)
            _fsync_directory(stage.destination.parent_path)

    def _publish_stage(self, artifact_stage):
        stage_identity = _file_identity(os.lstat(artifact_stage.stage_path))
        pending = PublishedArtifact(
            role=artifact_stage.role,
            final_path=artifact_stage.destination.final_path,
            identity=stage_identity,
        )
        try:
            _link_no_clobber(
                artifact_stage.stage_path,
                artifact_stage.destination.final_path,
            )
        except BaseException:
            try:
                final_identity = _file_identity(os.lstat(pending.final_path))
            except BaseException:
                pass
            else:
                if final_identity == stage_identity:
                    self._pending_publication = pending
            raise

        self._pending_publication = pending
        final_identity = _file_identity(os.lstat(pending.final_path))
        if final_identity != stage_identity:
            raise LyraArtifactConsistencyError(
                category="publication_identity",
                field=artifact_stage.role,
                expected=stage_identity,
                actual=final_identity,
            )
        self._published.append(pending)
        self._pending_publication = None

    def _remove_stage_and_sync(self, artifact_stage):
        _unlink_path(artifact_stage.stage_path)
        _fsync_directory(artifact_stage.destination.parent_path)

    def _publish(self):
        normalized_stage, bam_stage, summary_stage = (
            self._stages_in_publication_order()
        )

        self.stage = "publish_normalized"
        self._publish_stage(normalized_stage)
        self.stage = "sync_normalized_final"
        _fsync_directory(normalized_stage.destination.parent_path)

        self.stage = "publish_viral_bam"
        self._publish_stage(bam_stage)
        self.stage = "sync_viral_bam_final"
        _fsync_directory(bam_stage.destination.parent_path)

        self.stage = "remove_normalized_stage"
        self._remove_stage_and_sync(normalized_stage)
        self.stage = "remove_viral_bam_stage"
        self._remove_stage_and_sync(bam_stage)

        self.stage = "publish_summary"
        self._publish_stage(summary_stage)
        self.stage = "remove_summary_stage"
        _unlink_path(summary_stage.stage_path)
        self.stage = "sync_summary_final"
        _fsync_directory(summary_stage.destination.parent_path)
        self.stage = "complete"

    def _rollback_published(self, artifact):
        try:
            observed = _file_identity(os.lstat(artifact.final_path))
        except FileNotFoundError:
            return CleanupOutcome(
                operation="rollback",
                role=artifact.role,
                path=artifact.final_path,
                status="absent",
                expected=artifact.identity,
            )
        except BaseException as error:
            return _cleanup_failure(
                "rollback",
                artifact.role,
                artifact.final_path,
                "lstat_failed",
                error,
            )

        if observed != artifact.identity:
            return CleanupOutcome(
                operation="rollback",
                role=artifact.role,
                path=artifact.final_path,
                status="identity_mismatch",
                expected=artifact.identity,
                actual=observed,
            )

        try:
            _unlink_path(artifact.final_path)
        except FileNotFoundError:
            return CleanupOutcome(
                operation="rollback",
                role=artifact.role,
                path=artifact.final_path,
                status="absent",
                expected=artifact.identity,
                actual=observed,
            )
        except BaseException as error:
            return _cleanup_failure(
                "rollback",
                artifact.role,
                artifact.final_path,
                "unlink_failed",
                error,
            )

        try:
            _fsync_directory(os.path.dirname(artifact.final_path))
        except BaseException as error:
            return _cleanup_failure(
                "rollback",
                artifact.role,
                artifact.final_path,
                "removed_sync_failed",
                error,
            )
        return CleanupOutcome(
            operation="rollback",
            role=artifact.role,
            path=artifact.final_path,
            status="removed",
            expected=artifact.identity,
            actual=observed,
        )

    def _cleanup_stage(self, artifact_stage):
        try:
            _unlink_path(artifact_stage.stage_path)
        except FileNotFoundError:
            return CleanupOutcome(
                operation="cleanup",
                role=artifact_stage.role,
                path=artifact_stage.stage_path,
                status="absent",
            )
        except BaseException as error:
            return _cleanup_failure(
                "cleanup",
                artifact_stage.role,
                artifact_stage.stage_path,
                "unlink_failed",
                error,
            )

        try:
            _fsync_directory(artifact_stage.destination.parent_path)
        except BaseException as error:
            return _cleanup_failure(
                "cleanup",
                artifact_stage.role,
                artifact_stage.stage_path,
                "removed_sync_failed",
                error,
            )
        return CleanupOutcome(
            operation="cleanup",
            role=artifact_stage.role,
            path=artifact_stage.stage_path,
            status="removed",
        )

    def rollback_and_cleanup(self):
        """Rollback owned finals and remove all known stages without raising."""
        outcomes = []
        rollback_artifacts = []
        if self._pending_publication is not None:
            rollback_artifacts.append(self._pending_publication)
        rollback_artifacts.extend(reversed(self._published))
        for artifact in rollback_artifacts:
            outcomes.append(self._rollback_published(artifact))
        for artifact_stage in self._stages:
            outcomes.append(self._cleanup_stage(artifact_stage))
        outcomes.extend(self._close_owned_descriptors())
        self._cleanup_outcomes = tuple(outcomes)
        return self._cleanup_outcomes

    def _cleanup_failures(self):
        failures = []
        for outcome in self._cleanup_outcomes:
            if outcome.status in {"absent", "removed"}:
                continue
            if outcome.status == "identity_mismatch":
                failures.append(
                    CleanupFailure(
                        operation=outcome.operation,
                        role=outcome.role,
                        path=outcome.path,
                        error_type="FileIdentityMismatch",
                        category="rollback_identity_mismatch",
                    )
                )
                continue

            prefix = (
                "stage_cleanup"
                if outcome.operation == "cleanup"
                else "rollback"
            )
            failures.append(
                CleanupFailure(
                    operation=outcome.operation,
                    role=outcome.role,
                    path=outcome.path,
                    error_type=outcome.error_type,
                    errno=(
                        outcome.errno
                        if type(outcome.errno) is int
                        else None
                    ),
                    category="{}_{}".format(prefix, outcome.status),
                )
            )
        return tuple(failures)

    def _failure_stage(self, primary):
        if isinstance(primary, (LyraSourceIdentityError, LyraPathError)):
            primary_stage = vars(primary).get("stage")
            if type(primary_stage) is str:
                return primary_stage
        return self.stage

    def generate_validate_and_publish(self):
        try:
            self.stage = "pre_generation"
            _assert_path_plan_available(self.path_plan, self.stage)
            _assert_source_bam_identity(self.store, self.stage)
            self._create_stages()
            self._generate_and_validate()
            self._sync_stages()
            self._publish()
        except BaseException as primary:
            failure_stage = self._failure_stage(primary)
            self.rollback_and_cleanup()
            if not isinstance(primary, Exception):
                raise
            raise LyraPublicationError(
                stage=failure_stage,
                primary=primary,
                cleanup_failures=(
                    self._cleanup_failures()
                    + _bam_filter_cleanup_failures(primary)
                ),
            ) from primary
        else:
            close_outcomes = self._close_owned_descriptors()
            if close_outcomes:
                self._cleanup_outcomes = close_outcomes
                primary = OSError("failed to close artifact transaction descriptors")
                failure_stage = "close_transaction_descriptors"
                self.rollback_and_cleanup()
                raise LyraPublicationError(
                    stage=failure_stage,
                    primary=primary,
                    cleanup_failures=self._cleanup_failures(),
                ) from primary


def write_lyra_artifacts(
    store,
    path_plan,
    work_dir=None,
):
    """Generate, read back, and publish one coordinated artifact set."""
    transaction = LyraArtifactTransaction(
        store,
        path_plan,
        work_dir=work_dir,
    )
    transaction.generate_validate_and_publish()


class LyraFragmentStore:
    """Context-bound SQLite state for bounded Lyra reconciliation."""

    def __init__(self, sample_id, work_dir=None):
        self._sample_id = validate_sample_id(sample_id)
        self._score_records = 0
        self._input_bam_records = 0
        self._eligible_bam_records = 0
        self._counts = None
        self._threshold = None
        self._source_bam_fd = None
        self._source_bam_path = None
        self._source_bam_display_path = None
        self._source_bam_identity = None
        self._source_bam_header = None
        self._closed = False
        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="lyra_reconciliation_",
            dir=work_dir,
        )
        self._database_path = os.path.join(
            self._temporary_directory.name,
            "reconciliation.sqlite3",
        )
        self._connection = None
        try:
            self._connection = sqlite3.connect(self._database_path)
            self._connection.row_factory = sqlite3.Row
            self._configure_database()
            self._create_schema()
        except BaseException:
            self.close()
            raise

    @property
    def sample_id(self):
        return self._sample_id

    @property
    def database_path(self):
        return self._database_path

    @property
    def score_records(self):
        return self._score_records

    @property
    def input_bam_records(self):
        return self._input_bam_records

    @property
    def eligible_bam_records(self):
        return self._eligible_bam_records

    @property
    def counts(self):
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        if self._counts is None:
            raise RuntimeError("Lyra fragment store is not finalized")
        return self._counts

    @property
    def threshold(self):
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        if self._threshold is None:
            raise RuntimeError("Lyra fragment store is not finalized")
        return self._threshold

    def _require_source_metadata(self):
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        if self._source_bam_fd is None or self._source_bam_identity is None:
            raise RuntimeError("Lyra fragment store is not finalized")

    @property
    def source_bam_fd(self):
        self._require_source_metadata()
        return self._source_bam_fd

    @property
    def source_bam_path(self):
        self._require_source_metadata()
        return self._source_bam_path

    @property
    def source_bam_display_path(self):
        self._require_source_metadata()
        return self._source_bam_display_path

    @property
    def source_bam_identity(self):
        self._require_source_metadata()
        return self._source_bam_identity

    @property
    def source_bam_header(self):
        self._require_source_metadata()
        return copy.deepcopy(self._source_bam_header)

    def _install_source_metadata(
        self,
        source_bam_fd,
        source_bam_path,
        source_bam_display_path,
        source_bam_identity,
        source_bam_header,
    ):
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        if self._counts is None or self._threshold is None:
            raise RuntimeError("Lyra fragment store is not finalized")
        if self._source_bam_fd is not None or self._source_bam_identity is not None:
            raise RuntimeError("Lyra source metadata is already installed")
        if type(source_bam_fd) is not int or source_bam_fd < 0:
            raise ValueError("source_bam_fd must be a non-negative built-in int")
        os.fstat(source_bam_fd)

        retained_header = copy.deepcopy(source_bam_header)
        (
            self._source_bam_fd,
            self._source_bam_path,
            self._source_bam_display_path,
            self._source_bam_identity,
            self._source_bam_header,
        ) = (
            source_bam_fd,
            source_bam_path,
            source_bam_display_path,
            source_bam_identity,
            retained_header,
        )

    def _configure_database(self):
        self._connection.execute("PRAGMA journal_mode = OFF")
        self._connection.execute("PRAGMA synchronous = OFF")
        self._connection.execute("PRAGMA temp_store = FILE")
        self._connection.execute("PRAGMA cache_size = -64000")

    def _create_schema(self):
        self._connection.executescript(
            """
            CREATE TABLE evidence (
                read_id_key BLOB PRIMARY KEY,
                score_count INTEGER NOT NULL DEFAULT 0,
                score_1_text TEXT,
                score_1_line INTEGER,
                score_2_text TEXT,
                score_2_line INTEGER,
                eligible_bam_count INTEGER NOT NULL DEFAULT 0,
                unpaired_none_count INTEGER NOT NULL DEFAULT 0,
                unpaired_r1_count INTEGER NOT NULL DEFAULT 0,
                unpaired_r2_count INTEGER NOT NULL DEFAULT 0,
                unpaired_both_count INTEGER NOT NULL DEFAULT 0,
                paired_none_count INTEGER NOT NULL DEFAULT 0,
                paired_r1_count INTEGER NOT NULL DEFAULT 0,
                paired_r2_count INTEGER NOT NULL DEFAULT 0,
                paired_both_count INTEGER NOT NULL DEFAULT 0
            ) WITHOUT ROWID;

            CREATE TABLE fragments (
                read_id_key BLOB PRIMARY KEY,
                n_scores INTEGER NOT NULL,
                pairing TEXT NOT NULL,
                min_score_text TEXT NOT NULL,
                max_score_text TEXT NOT NULL,
                threshold_text TEXT NOT NULL,
                call TEXT NOT NULL
            ) WITHOUT ROWID;
            """
        )
        self._connection.commit()

    def _ordered_evidence_cursor(self):
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        return self._connection.execute(
            "SELECT * FROM evidence ORDER BY read_id_key ASC"
        )

    def iter_fragments(self):
        """Return a fresh streaming cursor of fragments in exact key order."""
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        cursor = self._connection.execute(
            """
            SELECT
                read_id_key,
                n_scores,
                pairing,
                min_score_text,
                max_score_text,
                threshold_text,
                call
            FROM fragments
            ORDER BY read_id_key ASC
            """
        )

        def _iterate():
            try:
                for row in cursor:
                    yield LyraFragment(
                        read_id=row["read_id_key"].decode(
                            "utf-8",
                            errors="strict",
                        ),
                        n_scores=row["n_scores"],
                        pairing=row["pairing"],
                        min_score=Decimal(row["min_score_text"]),
                        max_score=Decimal(row["max_score_text"]),
                        threshold=Decimal(row["threshold_text"]),
                        call=row["call"],
                    )
            finally:
                cursor.close()

        return _iterate()

    def iter_viral_read_ids(self):
        """Return a fresh streaming cursor of exact Viral read IDs."""
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        cursor = self._connection.execute(
            """
            SELECT read_id_key
            FROM fragments
            WHERE call = ?
            ORDER BY read_id_key ASC
            """,
            (CALL_VIRAL,),
        )

        def _iterate():
            try:
                for row in cursor:
                    yield row["read_id_key"].decode("utf-8", errors="strict")
            finally:
                cursor.close()

        return _iterate()

    def close(self):
        if self._closed:
            return
        self._closed = True
        cleanup_error = None
        cleanup_traceback = None

        source_bam_fd = self._source_bam_fd
        self._source_bam_fd = None
        connection = self._connection
        self._connection = None
        cleanup_operations = []
        if source_bam_fd is not None:
            cleanup_operations.append(lambda: os.close(source_bam_fd))
        if connection is not None:
            cleanup_operations.append(connection.close)
        cleanup_operations.append(self._temporary_directory.cleanup)

        for cleanup in cleanup_operations:
            try:
                cleanup()
            except BaseException as exc:
                if cleanup_error is None:
                    cleanup_error = exc
                    cleanup_traceback = exc.__traceback__

        if cleanup_error is not None:
            raise cleanup_error.with_traceback(cleanup_traceback)

    def __enter__(self):
        if self._closed:
            raise RuntimeError("Lyra fragment store is closed")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False


def _collect_score_evidence(store, score_path, sample_id):
    """Stream validated native score evidence into ``store``."""
    for record in iter_lyra_score_records(score_path, sample_id):
        store._connection.execute(
            _SCORE_EVIDENCE_UPSERT,
            (
                _read_id_key(record.read_id),
                1,
                record.score_text,
                record.line_number,
            ),
        )
        store._score_records += 1
        if store._score_records % SQLITE_COMMIT_INTERVAL == 0:
            store._connection.commit()
    store._connection.commit()


def _collect_bam_evidence(store, bam, display_path=None):
    """Stream exactly Lyra-eligible BAM role evidence into ``store``."""
    if display_path is None:
        display_path = bam
        path = os.fspath(bam)
        with pysam.AlignmentFile(path, "rb", check_sq=False) as opened_bam:
            return _collect_bam_evidence(store, opened_bam, display_path)

    for read in bam.fetch(until_eof=True):
        store._input_bam_records += 1
        if read.is_secondary or read.is_supplementary:
            continue
        sequence = read.query_sequence
        if sequence is None or len(sequence) < 50:
            continue
        if read.query_name is None:
            raise LyraInputError(
                category="bam_read_id",
                path=display_path,
                field="read_id",
                reason="eligible BAM record must have a query name",
            )

        role_index = _BAM_ROLE_INDEX[
            (read.is_paired, read.is_read1, read.is_read2)
        ]
        role_counts = tuple(
            1 if index == role_index else 0
            for index in range(len(_BAM_ROLE_INDEX))
        )
        store._connection.execute(
            _BAM_EVIDENCE_UPSERT,
            (_read_id_key(read.query_name), 1, *role_counts),
        )
        store._eligible_bam_records += 1
        if store._eligible_bam_records % SQLITE_COMMIT_INTERVAL == 0:
            store._connection.commit()
    store._connection.commit()


def _source_identity_status(error):
    if isinstance(error, FileNotFoundError):
        return "missing"
    return "unreadable"


def _assert_source_bam_identity(store, stage):
    """Require the retained descriptor and path to match reconciled identity."""
    expected = store.source_bam_identity
    display_path = store.source_bam_display_path
    try:
        descriptor_actual = _file_identity(os.fstat(store.source_bam_fd))
    except OSError as exc:
        raise LyraSourceIdentityError(
            stage=stage,
            path=display_path,
            expected=expected,
            actual_status=_source_identity_status(exc),
        ) from exc

    if descriptor_actual != expected:
        raise LyraSourceIdentityError(
            stage=stage,
            path=display_path,
            expected=expected,
            actual=descriptor_actual,
        )

    try:
        path_actual = _file_identity(os.stat(store.source_bam_path))
    except OSError as exc:
        raise LyraSourceIdentityError(
            stage=stage,
            path=display_path,
            expected=expected,
            actual_status=_source_identity_status(exc),
        ) from exc

    if path_actual != expected:
        raise LyraSourceIdentityError(
            stage=stage,
            path=display_path,
            expected=expected,
            actual=path_actual,
        )


def _validate_reconciliation(store, score_path, bam_path):
    """Raise one deterministic error unless all stored evidence is coherent."""
    converted_score_path = os.fspath(score_path)
    converted_bam_path = os.fspath(bam_path)
    for row in store._ordered_evidence_cursor():
        score_count = row["score_count"]
        eligible_bam_count = row["eligible_bam_count"]
        bam_role_counts = tuple(row[column] for column in _BAM_ROLE_COLUMNS)
        _, category = _pairing_state(eligible_bam_count, bam_role_counts)

        if category is None:
            if score_count > 0 and eligible_bam_count == 0:
                category = "score_only"
            elif eligible_bam_count > 0 and score_count == 0:
                category = "bam_only"
            elif score_count > 0 and eligible_bam_count > 0:
                if score_count != eligible_bam_count:
                    category = "record_count_mismatch"

        if category is not None:
            score_line_numbers = tuple(
                line_number
                for line_number in (row["score_1_line"], row["score_2_line"])
                if line_number is not None
            )
            raise LyraReconciliationError(
                category=category,
                score_path=converted_score_path,
                bam_path=converted_bam_path,
                read_id=row["read_id_key"].decode("utf-8", errors="strict"),
                score_count=score_count,
                eligible_bam_count=eligible_bam_count,
                score_line_numbers=score_line_numbers,
                bam_role_counts=bam_role_counts,
                reason=_RECONCILIATION_REASONS[category],
            )


def _finalize_fragments(store, threshold):
    """Finalize one exact Decimal classification for every evidence key."""
    fragments = 0
    single_end_fragments = 0
    complete_pair_fragments = 0
    incomplete_pair_fragments = 0
    viral_fragment_calls = 0
    nonviral_fragment_calls = 0
    threshold_text = str(threshold)

    for row in store._ordered_evidence_cursor():
        bam_role_counts = tuple(row[column] for column in _BAM_ROLE_COLUMNS)
        pairing, intrinsic_category = _pairing_state(
            row["eligible_bam_count"],
            bam_role_counts,
        )
        assert intrinsic_category is None

        n_scores = row["score_count"]
        score_texts = [row["score_1_text"]]
        if n_scores == 2:
            score_texts.append(row["score_2_text"])
        scores = [
            (Decimal(text), _canonical_decimal_text(text))
            for text in score_texts
        ]
        minimum = min(value for value, _ in scores)
        maximum = max(value for value, _ in scores)
        min_score_text = min(
            text for value, text in scores if value == minimum
        )
        max_score_text = min(
            text for value, text in scores if value == maximum
        )

        if pairing == PAIRING_SINGLE_END:
            is_viral = scores[0][0] >= threshold
            single_end_fragments += 1
        elif pairing == PAIRING_COMPLETE:
            is_viral = all(value >= threshold for value, _ in scores)
            complete_pair_fragments += 1
        else:
            assert pairing == PAIRING_INCOMPLETE
            is_viral = False
            incomplete_pair_fragments += 1
        call = CALL_VIRAL if is_viral else CALL_NONVIRAL

        store._connection.execute(
            """
            INSERT INTO fragments (
                read_id_key,
                n_scores,
                pairing,
                min_score_text,
                max_score_text,
                threshold_text,
                call
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                row["read_id_key"],
                n_scores,
                pairing,
                min_score_text,
                max_score_text,
                threshold_text,
                call,
            ),
        )
        fragments += 1
        if is_viral:
            viral_fragment_calls += 1
        else:
            nonviral_fragment_calls += 1
        if fragments % SQLITE_COMMIT_INTERVAL == 0:
            store._connection.commit()

    store._connection.commit()
    store._counts = LyraReconciliationCounts(
        input_bam_records=store.input_bam_records,
        eligible_bam_records=store.eligible_bam_records,
        score_records=store.score_records,
        fragments=fragments,
        single_end_fragments=single_end_fragments,
        complete_pair_fragments=complete_pair_fragments,
        incomplete_pair_fragments=incomplete_pair_fragments,
        viral_fragment_calls=viral_fragment_calls,
        nonviral_fragment_calls=nonviral_fragment_calls,
    )
    store._threshold = threshold


@contextmanager
def reconcile_lyra_fragments(
    score_path,
    bam_path,
    sample_id,
    threshold,
    work_dir=None,
):
    """Yield one fully validated and finalized file-backed fragment store."""
    iter_lyra_score_records(score_path, sample_id)
    validated_threshold = validate_lyra_threshold(threshold)
    source_bam_display_path = os.fspath(bam_path)
    source_bam_path = os.path.abspath(source_bam_display_path)
    with LyraFragmentStore(sample_id=sample_id, work_dir=work_dir) as store:
        _collect_score_evidence(store, score_path, sample_id)
        retained_source_fd = None
        try:
            with open(source_bam_path, "rb") as source_stream:
                source_identity = _file_identity(os.fstat(source_stream.fileno()))
                with pysam.AlignmentFile(
                    source_stream,
                    "rb",
                    check_sq=False,
                ) as bam:
                    source_header = bam.header.to_dict()
                    _collect_bam_evidence(
                        store,
                        bam,
                        source_bam_display_path,
                    )
                descriptor_identity = _file_identity(
                    os.fstat(source_stream.fileno())
                )
                retained_source_fd = os.dup(source_stream.fileno())
                retained_descriptor_identity = _file_identity(
                    os.fstat(retained_source_fd)
                )

            try:
                path_identity = _file_identity(os.stat(source_bam_path))
            except OSError as exc:
                raise LyraSourceIdentityError(
                    stage="reconciliation",
                    path=source_bam_display_path,
                    expected=source_identity,
                    actual_status=_source_identity_status(exc),
                ) from exc

            actual_identity = None
            for candidate_identity in (
                descriptor_identity,
                retained_descriptor_identity,
                path_identity,
            ):
                if candidate_identity != source_identity:
                    actual_identity = candidate_identity
                    break
            if actual_identity is not None:
                raise LyraSourceIdentityError(
                    stage="reconciliation",
                    path=source_bam_display_path,
                    expected=source_identity,
                    actual=actual_identity,
                )

            _validate_reconciliation(store, score_path, source_bam_display_path)
            _finalize_fragments(store, validated_threshold)
            store._install_source_metadata(
                source_bam_fd=retained_source_fd,
                source_bam_path=source_bam_path,
                source_bam_display_path=source_bam_display_path,
                source_bam_identity=source_identity,
                source_bam_header=source_header,
            )
            retained_source_fd = None
        finally:
            if retained_source_fd is not None:
                os.close(retained_source_fd)
        yield store


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
    with open(path, "rb") as compressed_stream:
        _validate_complete_zstd_frames_stream(compressed_stream)


def _validate_complete_zstd_frames_stream(compressed_stream):
    """Require one or more complete Zstandard frames from a binary stream."""
    decompressor = None
    saw_complete_frame = False
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
