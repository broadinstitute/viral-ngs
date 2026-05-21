from unittest.mock import MagicMock

import pandas as pd
import pytest

from viral_ngs.classify import virnucpro
from viral_ngs.core import file as util_file


def test_classify_contigs_writes_sorted_contig_calls(tmp_path):
    highestscore_tsv = tmp_path / "highestscore.tsv"
    output_tsv = tmp_path / "contigs.tsv"
    highestscore_tsv.write_text(
        "Modified_ID\tIs_Virus\tmax_score_0\tmax_score_1\n"
        "NODE_10_length_100_chunk_1\tTrue\t0.10\t0.95\n"
        "NODE_10_length_100_chunk_2\tTrue\t0.20\t0.90\n"
        "NODE_2_length_100_chunk_1\tFalse\t0.95\t0.10\n"
        "NODE_2_length_100_chunk_2\tFalse\t0.90\t0.20\n"
        "NODE_3_length_100_chunk_1\tFalse\t0.75\t0.75\n"
    )

    virnucpro.classify_contigs(str(highestscore_tsv), str(output_tsv), min_chunks=1)

    result = pd.read_csv(output_tsv, sep="\t")
    assert list(result["ID"]) == [
        "NODE_2_length_100",
        "NODE_3_length_100",
        "NODE_10_length_100",
    ]
    calls = dict(zip(result["ID"], result["call"]))
    tiers = dict(zip(result["ID"], result["tier"]))
    assert calls["NODE_10_length_100"] == "Viral"
    assert tiers["NODE_10_length_100"] == "high_confidence"
    assert calls["NODE_2_length_100"] == "Non-viral"
    assert tiers["NODE_2_length_100"] == "high_confidence"
    assert calls["NODE_3_length_100"] == "Ambiguous"
    assert tiers["NODE_3_length_100"] == "review"


def test_classify_contigs_rejects_missing_required_columns(tmp_path):
    highestscore_tsv = tmp_path / "highestscore.tsv"
    highestscore_tsv.write_text("Modified_ID\tmax_score_0\nNODE_1_chunk_1\t0.1\n")

    with pytest.raises(ValueError, match="Missing required columns"):
        virnucpro.classify_contigs(
            str(highestscore_tsv),
            str(tmp_path / "contigs.tsv"),
        )


def test_classify_contigs_rejects_unmatched_ids(tmp_path):
    highestscore_tsv = tmp_path / "highestscore.tsv"
    highestscore_tsv.write_text(
        "Modified_ID\tIs_Virus\tmax_score_0\tmax_score_1\n"
        "contigA_chunk_1\tTrue\t0.1\t0.9\n"
    )

    with pytest.raises(ValueError, match="No valid IDs extracted"):
        virnucpro.classify_contigs(
            str(highestscore_tsv),
            str(tmp_path / "contigs.tsv"),
        )


def test_classify_reads_by_contig_writes_read_level_calls(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    paf_file = tmp_path / "reads.paf"
    output_tsv = tmp_path / "reads_classified.tsv.zst"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
        "NODE_1\tViral\thigh_confidence\t0.9\t5\t5\t0\t0\t1.0\t0.0\n"
        "NODE_2\tNon-viral\thigh_confidence\t-0.9\t5\t0\t5\t0\t0.0\t1.0\n"
    )
    paf_file.write_text(
        # read1 maps well to a viral contig.
        "read1\t100\t0\t100\t+\tNODE_1\t1000\t0\t100\t95\t100\t60\t"
        "tp:A:P\tcg:Z:100M\t95.0\t100.0\n"
        # Secondary alignment is ignored.
        "read1\t100\t0\t100\t+\tNODE_2\t1000\t0\t100\t95\t100\t60\t"
        "tp:A:S\t95.0\t100.0\n"
        # read2 maps below the default mapq threshold to a non-viral contig.
        "read2\t100\t0\t100\t+\tNODE_2\t1000\t0\t100\t95\t100\t3\t"
        "tp:A:P\t95.0\t100.0\n"
        # read3 maps to contigs with distinct calls and is flagged Multi-mapped.
        "read3\t100\t0\t100\t+\tNODE_1\t1000\t0\t100\t90\t100\t20\t"
        "tp:A:P\t90.0\t100.0\n"
        "read3\t100\t0\t100\t-\tNODE_2\t1000\t0\t100\t98\t100\t30\t"
        "tp:A:P\t98.0\t100.0\n"
        # Missing contig classification is filled as Unclassified.
        "read4\t100\t0\t100\t+\tNODE_3\t1000\t0\t100\t91\t100\t30\t"
        "tp:A:P\t91.0\t100.0\n"
    )

    virnucpro.classify_reads_by_contig(
        str(paf_file),
        str(contig_classifications),
        str(output_tsv),
    )

    with util_file.open_or_gzopen(str(output_tsv), "rt") as inf:
        result = pd.read_csv(inf, sep="\t")

    assert list(result.columns) == virnucpro.READ_CLASSIFICATION_COLUMNS
    calls = dict(zip(result["read_id"], result["call"]))
    tiers = dict(zip(result["read_id"], result["tier"]))
    mapped_well = dict(zip(result["read_id"], result["mapped_well"]))
    contigs = dict(zip(result["read_id"], result["contig_id"]))

    assert calls == {
        "read1": "Viral",
        "read2": "Non-viral",
        "read3": "Multi-mapped",
        "read4": "Unclassified",
    }
    assert tiers["read3"] == "review"
    assert mapped_well == {
        "read1": True,
        "read2": False,
        "read3": True,
        "read4": True,
    }
    assert contigs["read3"] == "NODE_2"


def test_classify_reads_by_contig_detects_fractional_identity_scale(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    paf_file = tmp_path / "reads.paf"
    output_tsv = tmp_path / "reads_classified.tsv"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
        "NODE_1\tViral\thigh_confidence\t0.9\t5\t5\t0\t0\t1.0\t0.0\n"
    )
    paf_file.write_text(
        "read1\t100\t0\t100\t+\tNODE_1\t1000\t0\t100\t95\t100\t60\t"
        "tp:A:P\t0.95\t0.90\n"
    )

    virnucpro.classify_reads_by_contig(
        str(paf_file),
        str(contig_classifications),
        str(output_tsv),
        min_identity=90.0,
        min_query_cov=80.0,
    )

    result = pd.read_csv(output_tsv, sep="\t")

    assert result.loc[0, "mapped_well"]


def test_classify_reads_by_contig_writes_header_for_empty_paf(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    paf_file = tmp_path / "reads.paf"
    output_tsv = tmp_path / "reads_classified.tsv.zst"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
    )
    paf_file.write_text("")

    virnucpro.classify_reads_by_contig(
        str(paf_file),
        str(contig_classifications),
        str(output_tsv),
    )

    with util_file.open_or_gzopen(str(output_tsv), "rt") as inf:
        result = pd.read_csv(inf, sep="\t")

    assert list(result.columns) == virnucpro.READ_CLASSIFICATION_COLUMNS
    assert result.empty


def test_classify_reads_by_contig_rejects_missing_classification_columns(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    paf_file = tmp_path / "reads.paf"

    contig_classifications.write_text("ID\tcall\nNODE_1\tViral\n")
    paf_file.write_text(
        "read1\t100\t0\t100\t+\tNODE_1\t1000\t0\t100\t95\t100\t60\t"
        "tp:A:P\t95.0\t100.0\n"
    )

    with pytest.raises(ValueError, match="Missing required columns"):
        virnucpro.classify_reads_by_contig(
            str(paf_file),
            str(contig_classifications),
            str(tmp_path / "reads_classified.tsv"),
        )


# -- _default_memory_limit cgroup auto-detect -------------------------------


def _patch_cgroup_paths(monkeypatch, v2_path=None, v1_path=None):
    """Point the helper at tmp paths instead of /sys/fs/cgroup."""
    monkeypatch.setattr(
        virnucpro, "_CGROUP_V2_PATH", str(v2_path) if v2_path else "/nonexistent/v2",
    )
    monkeypatch.setattr(
        virnucpro, "_CGROUP_V1_PATH", str(v1_path) if v1_path else "/nonexistent/v1",
    )


def test_default_memory_limit_reads_cgroup_v2_byte_count(tmp_path, monkeypatch):
    v2 = tmp_path / "memory.max"
    v2.write_text("8589934592\n")  # 8 GiB
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    # 8 GiB * 0.75 = 6 GiB = 6144 MiB
    assert virnucpro._default_memory_limit() == "6144MB"


def test_default_memory_limit_returns_none_when_cgroup_v2_unlimited(tmp_path, monkeypatch):
    v2 = tmp_path / "memory.max"
    v2.write_text("max\n")
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    assert virnucpro._default_memory_limit() is None


def test_default_memory_limit_falls_back_to_cgroup_v1(tmp_path, monkeypatch):
    v1 = tmp_path / "memory.limit_in_bytes"
    v1.write_text("4294967296\n")  # 4 GiB
    _patch_cgroup_paths(monkeypatch, v1_path=v1)

    # 4 GiB * 0.75 = 3 GiB = 3072 MiB
    assert virnucpro._default_memory_limit() == "3072MB"


def test_default_memory_limit_returns_none_for_cgroup_v1_unlimited_sentinel(tmp_path, monkeypatch):
    v1 = tmp_path / "memory.limit_in_bytes"
    v1.write_text("{}\n".format(virnucpro._CGROUP_V1_UNLIMITED))
    _patch_cgroup_paths(monkeypatch, v1_path=v1)

    assert virnucpro._default_memory_limit() is None


def test_default_memory_limit_returns_none_when_no_cgroup_files(monkeypatch):
    _patch_cgroup_paths(monkeypatch)
    assert virnucpro._default_memory_limit() is None


def test_default_memory_limit_returns_none_on_malformed_content(tmp_path, monkeypatch):
    v2 = tmp_path / "memory.max"
    v2.write_text("not-a-number\n")
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    assert virnucpro._default_memory_limit() is None


def test_default_memory_limit_returns_none_for_zero_byte_limit(tmp_path, monkeypatch):
    v2 = tmp_path / "memory.max"
    v2.write_text("0\n")
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    assert virnucpro._default_memory_limit() is None


def test_connect_duckdb_caller_value_wins_over_autodetect(tmp_path, monkeypatch):
    """Explicit caller-supplied limit must not be overridden by cgroup auto-detect."""
    v2 = tmp_path / "memory.max"
    v2.write_text("8589934592\n")  # would auto-detect to 6144MB
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    fake_conn = MagicMock()
    fake_duckdb = MagicMock()
    fake_duckdb.connect.return_value = fake_conn

    virnucpro._connect_duckdb(
        fake_duckdb,
        duckdb_temp_dir=str(tmp_path),
        duckdb_memory_limit="2GB",
    )

    config = fake_duckdb.connect.call_args.kwargs["config"]
    assert config["memory_limit"] == "2GB"


def test_connect_duckdb_empty_string_opts_out_of_limit(tmp_path, monkeypatch):
    """Empty-string limit must skip both the override and the auto-detect."""
    v2 = tmp_path / "memory.max"
    v2.write_text("8589934592\n")
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    fake_duckdb = MagicMock()
    virnucpro._connect_duckdb(
        fake_duckdb,
        duckdb_temp_dir=str(tmp_path),
        duckdb_memory_limit="",
    )

    config = fake_duckdb.connect.call_args.kwargs["config"]
    assert "memory_limit" not in config


def test_connect_duckdb_autodetect_when_caller_passes_none(tmp_path, monkeypatch):
    """When caller passes None, helper should call _default_memory_limit() and use it."""
    v2 = tmp_path / "memory.max"
    v2.write_text("8589934592\n")
    _patch_cgroup_paths(monkeypatch, v2_path=v2)

    fake_duckdb = MagicMock()
    virnucpro._connect_duckdb(
        fake_duckdb,
        duckdb_temp_dir=str(tmp_path),
        duckdb_memory_limit=None,
    )

    config = fake_duckdb.connect.call_args.kwargs["config"]
    assert config["memory_limit"] == "6144MB"
