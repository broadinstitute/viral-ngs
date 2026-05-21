import pandas as pd
import pytest

from viral_ngs.classify import virnucpro


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
