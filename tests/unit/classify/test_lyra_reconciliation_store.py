from dataclasses import fields
from dataclasses import FrozenInstanceError
from decimal import Decimal
from pathlib import Path
import sqlite3

import pytest

from viral_ngs.classify import lyra


_EVIDENCE_COLUMNS = [
    "read_id_key",
    "score_count",
    "score_1_text",
    "score_1_line",
    "score_2_text",
    "score_2_line",
    "eligible_bam_count",
    "unpaired_none_count",
    "unpaired_r1_count",
    "unpaired_r2_count",
    "unpaired_both_count",
    "paired_none_count",
    "paired_r1_count",
    "paired_r2_count",
    "paired_both_count",
]

_FRAGMENT_COLUMNS = [
    "read_id_key",
    "n_scores",
    "pairing",
    "min_score_text",
    "max_score_text",
    "threshold_text",
    "call",
]


def _insert_empty_evidence(store, read_id):
    store._connection.execute(
        "INSERT INTO evidence (read_id_key) VALUES (?)",
        (lyra._read_id_key(read_id),),
    )


def test_contract_dataclasses_are_frozen_with_exact_fields():
    assert [field.name for field in fields(lyra.LyraFragment)] == [
        "read_id",
        "n_scores",
        "pairing",
        "min_score",
        "max_score",
        "threshold",
        "call",
    ]
    assert [field.name for field in fields(lyra.LyraReconciliationCounts)] == [
        "input_bam_records",
        "eligible_bam_records",
        "score_records",
        "fragments",
        "single_end_fragments",
        "complete_pair_fragments",
        "incomplete_pair_fragments",
        "viral_fragment_calls",
        "nonviral_fragment_calls",
    ]

    fragment = lyra.LyraFragment(
        read_id="read-1",
        n_scores=1,
        pairing="Single-end",
        min_score=Decimal("0.8"),
        max_score=Decimal("0.8"),
        threshold=Decimal("0.8"),
        call="Viral",
    )
    counts = lyra.LyraReconciliationCounts(0, 0, 0, 0, 0, 0, 0, 0, 0)

    with pytest.raises(FrozenInstanceError):
        fragment.call = "Non-viral"
    with pytest.raises(FrozenInstanceError):
        counts.fragments = 1


@pytest.mark.parametrize(
    ("threshold", "expected"),
    [
        (0, Decimal("0")),
        (1, Decimal("1")),
        (0.8, Decimal("0.8")),
        ("0.800000", Decimal("0.800000")),
        (Decimal("0.25"), Decimal("0.25")),
    ],
)
def test_threshold_accepts_finite_inclusive_decimal_values(threshold, expected):
    value = lyra.validate_lyra_threshold(threshold)

    assert value == expected
    assert isinstance(value, Decimal)


@pytest.mark.parametrize(
    "threshold",
    [
        True,
        False,
        None,
        "",
        "not-a-number",
        -0.0001,
        1.0001,
        float("nan"),
        float("inf"),
        float("-inf"),
        Decimal("NaN"),
        Decimal("Infinity"),
        Decimal("-Infinity"),
    ],
)
def test_threshold_rejects_boolean_malformed_nonfinite_and_out_of_range_values(
    threshold,
):
    with pytest.raises(lyra.LyraInputError) as exc_info:
        lyra.validate_lyra_threshold(threshold)

    error = exc_info.value
    assert error.category == "threshold"
    assert error.path is None
    assert error.line_number is None
    assert error.field == "threshold"
    assert error.reason


def test_read_id_key_is_exact_strict_utf8_and_has_explicit_order(tmp_path):
    read_ids = [
        "read-2",
        "é",
        "prefix!",
        "a",
        "prefix",
        "A",
        "e\u0301",
        "read-10",
    ]

    assert [lyra._read_id_key(read_id) for read_id in read_ids] == [
        read_id.encode("utf-8") for read_id in read_ids
    ]
    assert len(set(map(lyra._read_id_key, read_ids))) == len(read_ids)

    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        for read_id in reversed(read_ids):
            _insert_empty_evidence(store, read_id)
        ordered_ids = [
            row["read_id_key"].decode("utf-8")
            for row in store._ordered_evidence_cursor()
        ]

    assert ordered_ids == sorted(read_ids, key=lambda value: value.encode("utf-8"))


def test_sqlite_store_uses_exact_bounded_schemas_and_file_pragmas(tmp_path):
    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        database_path = Path(store.database_path)
        assert database_path.is_file()
        assert database_path.parent.parent == tmp_path

        evidence_info = store._connection.execute(
            "PRAGMA table_info(evidence)"
        ).fetchall()
        fragment_info = store._connection.execute(
            "PRAGMA table_info(fragments)"
        ).fetchall()
        table_sql = {
            row["name"]: row["sql"]
            for row in store._connection.execute(
                "SELECT name, sql FROM sqlite_master WHERE type = ?",
                ("table",),
            )
        }

        assert [row["name"] for row in evidence_info] == _EVIDENCE_COLUMNS
        assert [row["name"] for row in fragment_info] == _FRAGMENT_COLUMNS
        assert evidence_info[0]["type"] == "BLOB"
        assert evidence_info[0]["pk"] == 1
        assert fragment_info[0]["type"] == "BLOB"
        assert fragment_info[0]["pk"] == 1
        assert "WITHOUT ROWID" in table_sql["evidence"].upper()
        assert "WITHOUT ROWID" in table_sql["fragments"].upper()
        assert store._connection.execute("PRAGMA journal_mode").fetchone()[0] == "off"
        assert store._connection.execute("PRAGMA synchronous").fetchone()[0] == 0
        assert store._connection.execute("PRAGMA temp_store").fetchone()[0] == 1
        assert store._connection.execute("PRAGMA cache_size").fetchone()[0] < 0
        assert lyra.SQLITE_COMMIT_INTERVAL == 10000


def test_sqlite_store_cleanup_is_idempotent_and_rejects_iteration_after_close(
    tmp_path,
):
    store = lyra.LyraFragmentStore("sample", work_dir=tmp_path)
    connection = store._connection
    database_path = Path(store.database_path)
    temporary_directory = database_path.parent

    assert database_path.exists()
    store.close()
    store.close()

    assert not temporary_directory.exists()
    with pytest.raises(sqlite3.ProgrammingError):
        connection.execute("SELECT 1")
    with pytest.raises(RuntimeError, match="closed"):
        store._ordered_evidence_cursor()


def test_sqlite_store_cleans_up_after_normal_and_exceptional_context_exit(tmp_path):
    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as normal_store:
        normal_directory = Path(normal_store.database_path).parent
    assert not normal_directory.exists()

    with pytest.raises(RuntimeError, match="expected failure"):
        with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as failed_store:
            failed_directory = Path(failed_store.database_path).parent
            raise RuntimeError("expected failure")
    assert not failed_directory.exists()


def test_sql_looking_and_control_bearing_read_ids_remain_bound_values(tmp_path):
    read_ids = [
        "read'); DROP TABLE evidence; --",
        "read\x00control\x01value",
    ]

    with lyra.LyraFragmentStore("sample", work_dir=tmp_path) as store:
        for read_id in read_ids:
            _insert_empty_evidence(store, read_id)
        stored_keys = {
            row["read_id_key"] for row in store._ordered_evidence_cursor()
        }
        evidence_count = store._connection.execute(
            "SELECT COUNT(*) FROM evidence"
        ).fetchone()[0]

    assert stored_keys == {read_id.encode("utf-8") for read_id in read_ids}
    assert evidence_count == 2


def test_store_preserves_one_immutable_sample_scalar_for_unique_fragment_keys(tmp_path):
    sample_id = "échantillon exact"
    read_ids = ["read-1", "read-2", "Read-2"]

    with lyra.LyraFragmentStore(sample_id, work_dir=tmp_path) as store:
        assert store.sample_id is sample_id
        with pytest.raises(AttributeError):
            store.sample_id = "replacement"

        for read_id in read_ids:
            _insert_empty_evidence(store, read_id)
        fragment_keys = [
            (store.sample_id, row["read_id_key"].decode("utf-8"))
            for row in store._ordered_evidence_cursor()
        ]
        evidence_columns = {
            row["name"]
            for row in store._connection.execute("PRAGMA table_info(evidence)")
        }

    assert len(fragment_keys) == len(set(fragment_keys)) == len(read_ids)
    assert "sample_id" not in evidence_columns
