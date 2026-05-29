"""Unit tests for the shared DuckDB execution helpers."""

import pytest

from viral_ngs.core import duckdb_utils, misc


def test_import_duckdb_returns_module():
    duckdb = pytest.importorskip("duckdb")
    assert duckdb_utils.import_duckdb("unit test") is duckdb


def test_connect_none_consults_default(tmp_path, monkeypatch):
    pytest.importorskip("duckdb")
    calls = []
    monkeypatch.setattr(
        duckdb_utils, "default_memory_limit", lambda: calls.append(1) or "256MB"
    )
    con = duckdb_utils.connect(str(tmp_path))
    try:
        assert calls == [1]
    finally:
        con.close()


def test_connect_empty_string_opts_out(tmp_path, monkeypatch):
    pytest.importorskip("duckdb")
    # An empty string must neither consult the cgroup default nor set a limit.
    monkeypatch.setattr(
        duckdb_utils,
        "default_memory_limit",
        lambda: pytest.fail("default_memory_limit must not be called for ''"),
    )
    con = duckdb_utils.connect(str(tmp_path), memory_limit="")
    con.close()


def test_connect_explicit_limit_is_applied(tmp_path):
    pytest.importorskip("duckdb")
    # DuckDB normalizes the displayed unit (e.g. MB->MiB), so assert distinct
    # caller limits produce distinct settings rather than matching a literal.
    small = duckdb_utils.connect(str(tmp_path), memory_limit="256MB")
    large = duckdb_utils.connect(str(tmp_path), memory_limit="2GB")
    try:
        small_v = small.execute("SELECT current_setting('memory_limit')").fetchone()[0]
        large_v = large.execute("SELECT current_setting('memory_limit')").fetchone()[0]
        assert small_v != large_v
    finally:
        small.close()
        large.close()


def test_default_memory_limit_none_when_uncapped(monkeypatch):
    monkeypatch.setattr(misc, "available_mem_bytes", lambda fraction=1.0: None)
    assert duckdb_utils.default_memory_limit() is None


def test_default_memory_limit_formats_mib(monkeypatch):
    monkeypatch.setattr(
        misc, "available_mem_bytes", lambda fraction=1.0: 4 * 1024 * 1024 * 1024
    )
    assert duckdb_utils.default_memory_limit() == "{}MB".format(4 * 1024)


def test_available_mem_bytes_none_when_no_cgroup(monkeypatch):
    monkeypatch.setattr(misc, "_cgroup_mem_limit_bytes", lambda: None)
    assert misc.available_mem_bytes() is None
