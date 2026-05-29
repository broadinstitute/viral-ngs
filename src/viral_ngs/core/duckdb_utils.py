"""Shared DuckDB helpers.

DuckDB is a classify-image dependency, so this module imports it lazily and is
not imported from ``viral_ngs.core.__init__``.
"""

import logging

from viral_ngs.core import misc

log = logging.getLogger(__name__)

_DUCKDB_MEM_FRACTION = 0.75


def import_duckdb(context="this command"):
    try:
        import duckdb
    except ImportError as exc:
        raise ImportError(
            "DuckDB is required for {}. Use a viral-ngs image that includes "
            "the classify dependencies.".format(context)
        ) from exc
    return duckdb


def default_memory_limit():
    """Return a DuckDB memory-limit string derived from the cgroup, or None."""
    mem_bytes = misc.available_mem_bytes(fraction=_DUCKDB_MEM_FRACTION)
    if mem_bytes is None:
        return None

    mem_mib = mem_bytes // (1024 * 1024)
    if mem_mib <= 0:
        return None
    return "{}MB".format(mem_mib)


def connect(temp_directory, memory_limit=None, duckdb_module=None):
    """Open an in-memory DuckDB connection with shared resource defaults.

    Caller-supplied ``memory_limit`` values always win. Pass an empty string to
    opt out of setting a DuckDB memory limit.
    """
    duckdb = duckdb_module or import_duckdb("DuckDB-backed viral-ngs helpers")
    config = {"temp_directory": temp_directory}

    if memory_limit is None:
        memory_limit = default_memory_limit()
    if memory_limit:
        config["memory_limit"] = str(memory_limit)
        log.debug("DuckDB memory_limit set to %s", config["memory_limit"])

    return duckdb.connect(database=":memory:", config=config)
