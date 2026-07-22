"""Strict validation helpers for native Lyra score tables."""

import codecs
import gzip
import io
import os
import re
import sqlite3
import tempfile
import zlib
from dataclasses import dataclass
from decimal import Decimal
from decimal import InvalidOperation

import pysam
import zstandard as zstd

import viral_ngs.core.file as util_file


RENDERED_VALUE_CAP = 160
EXPECTED_SCORE_HEADER = ("read_id", "score", "call")
SQLITE_COMMIT_INTERVAL = 10000

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
    return value


def _read_id_key(read_id):
    return read_id.encode("utf-8", errors="strict")


class LyraFragmentStore:
    """Context-bound SQLite state for bounded Lyra reconciliation."""

    def __init__(self, sample_id, work_dir=None):
        self._sample_id = validate_sample_id(sample_id)
        self._score_records = 0
        self._input_bam_records = 0
        self._eligible_bam_records = 0
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

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            if self._connection is not None:
                self._connection.close()
                self._connection = None
        finally:
            self._temporary_directory.cleanup()

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


def _collect_bam_evidence(store, bam_path):
    """Stream exactly Lyra-eligible BAM role evidence into ``store``."""
    path = os.fspath(bam_path)
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
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
                    path=bam_path,
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
