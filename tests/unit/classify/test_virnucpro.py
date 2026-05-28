from unittest.mock import MagicMock

import pandas as pd
import pytest
import pysam

from viral_ngs.classify import virnucpro
from viral_ngs.core import file as util_file


_REF_IDS = {
    "NODE_1": 0,
    "NODE_2": 1,
    "NODE_3": 2,
}

_BAM_HEADER = {
    "HD": {"VN": "1.6"},
    "SQ": [
        {"SN": "NODE_1", "LN": 1000},
        {"SN": "NODE_2", "LN": 1000},
        {"SN": "NODE_3", "LN": 1000},
    ],
}


def _write_bam(path, records):
    with pysam.AlignmentFile(str(path), "wb", header=_BAM_HEADER) as outf:
        for record in records:
            outf.write(record)


def _segment(
    query_name,
    reference_name="NODE_1",
    reference_start=0,
    cigar=None,
    mapq=60,
    nm=0,
    flag=0,
    query_length=None,
):
    cigar = cigar or [(pysam.CMATCH, 100)]
    if query_length is None:
        query_length = sum(
            length
            for op, length in cigar
            if op in (
                pysam.CMATCH,
                pysam.CINS,
                pysam.CSOFT_CLIP,
                pysam.CEQUAL,
                pysam.CDIFF,
            )
        )
    if query_length == 0:
        query_length = 100

    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.query_sequence = "A" * query_length
    record.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    record.flag = flag

    if flag & 0x4:
        record.reference_id = -1
        record.reference_start = -1
        record.mapping_quality = 0
    else:
        record.reference_id = _REF_IDS[reference_name]
        record.reference_start = reference_start
        record.mapping_quality = mapq
        record.cigartuples = cigar
        if nm is not None:
            record.set_tag("NM", nm)

    return record


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


def test_classify_contigs_uses_custom_score_thresholds(tmp_path):
    highestscore_tsv = tmp_path / "highestscore.tsv"
    output_tsv = tmp_path / "contigs.tsv"
    highestscore_tsv.write_text(
        "Modified_ID\tIs_Virus\tmax_score_0\tmax_score_1\n"
        "NODE_1_chunk_1\tTrue\t0.25\t0.75\n"
        "NODE_2_chunk_1\tFalse\t0.65\t0.65\n"
        "NODE_3_chunk_1\tTrue\t0.45\t0.72\n"
    )

    virnucpro.classify_contigs(
        str(highestscore_tsv),
        str(output_tsv),
        min_chunks=1,
        min_confident_score=0.7,
        max_opposing_score=0.3,
        min_ambiguous_score=0.6,
        min_weighted_delta=0.25,
        high_confidence_delta=0.4,
    )

    result = pd.read_csv(output_tsv, sep="\t")
    result_by_id = result.set_index("ID")
    assert result_by_id.loc["NODE_1", "call"] == "Viral"
    assert result_by_id.loc["NODE_1", "tier"] == "high_confidence"
    assert result_by_id.loc["NODE_2", "n_ambiguous"] == 1
    assert result_by_id.loc["NODE_3", "call"] == "Viral"


@pytest.mark.parametrize(
    "contents",
    [
        "",
        "Modified_ID\tIs_Virus\tmax_score_0\tmax_score_1\n",
    ],
)
def test_classify_contigs_writes_header_for_empty_input(tmp_path, contents):
    highestscore_tsv = tmp_path / "highestscore.tsv"
    output_tsv = tmp_path / "contigs.tsv"
    highestscore_tsv.write_text(contents)

    virnucpro.classify_contigs(str(highestscore_tsv), str(output_tsv))

    result = pd.read_csv(output_tsv, sep="\t")
    assert list(result.columns) == virnucpro.CONTIG_CLASSIFICATION_COLUMNS
    assert result.empty


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
    aligned_bam = tmp_path / "reads.bam"
    output_tsv = tmp_path / "reads_classified.tsv.zst"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
        "NODE_1\tViral\thigh_confidence\t0.9\t5\t5\t0\t0\t1.0\t0.0\n"
        "NODE_2\tNon-viral\thigh_confidence\t-0.9\t5\t0\t5\t0\t0.0\t1.0\n"
    )
    _write_bam(
        aligned_bam,
        [
            # read1 maps well to a viral contig.
            _segment("read1", "NODE_1", mapq=60, nm=5),
            # Secondary alignment is ignored.
            _segment("read1", "NODE_2", mapq=60, nm=5, flag=0x100),
            # Supplementary alignment is ignored before NM validation.
            _segment("read1", "NODE_2", mapq=60, nm=None, flag=0x800),
            # read2 maps below the default mapq threshold to a non-viral contig.
            _segment("read2", "NODE_2", mapq=3, nm=5),
            # read3 maps to contigs with distinct calls and is flagged Multi-mapped.
            _segment("read3", "NODE_1", mapq=20, nm=10),
            _segment("read3", "NODE_2", mapq=30, nm=2, flag=0x10),
            # Missing contig classification is filled as Unclassified.
            _segment("read4", "NODE_3", mapq=30, nm=9),
            # Unmapped reads are not read-classification candidates.
            _segment("read5", flag=0x4),
        ],
    )

    virnucpro.classify_reads_by_contig(
        str(aligned_bam),
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


def test_classify_reads_by_contig_computes_bam_identity_and_query_cov(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"
    output_tsv = tmp_path / "reads_classified.tsv"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
        "NODE_1\tViral\thigh_confidence\t0.9\t5\t5\t0\t0\t1.0\t0.0\n"
    )
    _write_bam(
        aligned_bam,
        [
            _segment(
                "read1",
                "NODE_1",
                cigar=[
                    (pysam.CSOFT_CLIP, 10),
                    (pysam.CMATCH, 80),
                    (pysam.CSOFT_CLIP, 10),
                ],
                nm=8,
            ),
        ],
    )

    virnucpro.classify_reads_by_contig(
        str(aligned_bam),
        str(contig_classifications),
        str(output_tsv),
        min_identity=90.0,
        min_query_cov=80.0,
    )

    result = pd.read_csv(output_tsv, sep="\t")

    assert result.loc[0, "mapped_well"]
    assert result.loc[0, "pct_identity"] == pytest.approx(90.0)
    assert result.loc[0, "pct_query_cov"] == pytest.approx(80.0)


def test_classify_reads_by_contig_counts_indels_in_identity_denominator(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"
    output_tsv = tmp_path / "reads_classified.tsv"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
        "NODE_1\tViral\thigh_confidence\t0.9\t5\t5\t0\t0\t1.0\t0.0\n"
    )
    _write_bam(
        aligned_bam,
        [
            _segment(
                "read1",
                "NODE_1",
                cigar=[
                    (pysam.CMATCH, 50),
                    (pysam.CDEL, 5),
                    (pysam.CMATCH, 50),
                ],
                nm=5,
            ),
        ],
    )

    virnucpro.classify_reads_by_contig(
        str(aligned_bam),
        str(contig_classifications),
        str(output_tsv),
        min_identity=95.0,
        min_query_cov=100.0,
    )

    result = pd.read_csv(output_tsv, sep="\t")

    assert result.loc[0, "mapped_well"]
    assert result.loc[0, "pct_identity"] == pytest.approx(10000.0 / 105.0)
    assert result.loc[0, "pct_query_cov"] == pytest.approx(100.0)


def test_classify_reads_by_contig_writes_header_for_empty_bam(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"
    output_tsv = tmp_path / "reads_classified.tsv.zst"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
    )
    _write_bam(aligned_bam, [])

    virnucpro.classify_reads_by_contig(
        str(aligned_bam),
        str(contig_classifications),
        str(output_tsv),
    )

    with util_file.open_or_gzopen(str(output_tsv), "rt") as inf:
        result = pd.read_csv(inf, sep="\t")

    assert list(result.columns) == virnucpro.READ_CLASSIFICATION_COLUMNS
    assert result.empty


def test_classify_reads_by_contig_writes_header_for_all_unmapped_bam(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"
    output_tsv = tmp_path / "reads_classified.tsv"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
    )
    _write_bam(
        aligned_bam,
        [
            _segment("read1", flag=0x4),
            _segment("read2", flag=0x4),
        ],
    )

    virnucpro.classify_reads_by_contig(
        str(aligned_bam),
        str(contig_classifications),
        str(output_tsv),
    )

    result = pd.read_csv(output_tsv, sep="\t")

    assert list(result.columns) == virnucpro.READ_CLASSIFICATION_COLUMNS
    assert result.empty


def test_classify_reads_by_contig_rejects_missing_classification_columns(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"

    contig_classifications.write_text("ID\tcall\nNODE_1\tViral\n")
    _write_bam(aligned_bam, [_segment("read1", "NODE_1", mapq=60, nm=5)])

    with pytest.raises(ValueError, match="Missing required columns"):
        virnucpro.classify_reads_by_contig(
            str(aligned_bam),
            str(contig_classifications),
            str(tmp_path / "reads_classified.tsv"),
        )


def test_classify_reads_by_contig_rejects_mapped_bam_without_nm(tmp_path):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
        "NODE_1\tViral\thigh_confidence\t0.9\t5\t5\t0\t0\t1.0\t0.0\n"
    )
    _write_bam(aligned_bam, [_segment("read1", "NODE_1", nm=None)])

    with pytest.raises(ValueError, match="missing required NM tag"):
        virnucpro.classify_reads_by_contig(
            str(aligned_bam),
            str(contig_classifications),
            str(tmp_path / "reads_classified.tsv"),
        )


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"min_identity": 0.9}, "min_identity must be specified as a percent"),
        ({"min_query_cov": 0.8}, "min_query_cov must be specified as a percent"),
        ({"min_identity": -1}, "min_identity must be between 0 and 100"),
        ({"min_query_cov": 101}, "min_query_cov must be between 0 and 100"),
    ],
)
def test_classify_reads_by_contig_rejects_invalid_percent_thresholds(
        tmp_path, kwargs, match):
    contig_classifications = tmp_path / "contigs.tsv"
    aligned_bam = tmp_path / "reads.bam"

    contig_classifications.write_text(
        "ID\tcall\ttier\tweighted_delta\tn_chunks\tn_confident_viral\t"
        "n_confident_nonviral\tn_ambiguous\tviral_proportion\t"
        "nonviral_proportion\n"
    )
    _write_bam(aligned_bam, [])

    with pytest.raises(ValueError, match=match):
        virnucpro.classify_reads_by_contig(
            str(aligned_bam),
            str(contig_classifications),
            str(tmp_path / "reads_classified.tsv"),
            **kwargs
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
