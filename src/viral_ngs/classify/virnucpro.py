"""VirNucPro post-processing helpers."""

import logging
import os
import re

import pandas as pd

from viral_ngs.core import file

log = logging.getLogger(__name__)


def _classify_contig_group(
    group,
    min_viral_proportion=0.1,
    min_nonviral_proportion=0.1,
    min_chunk_count=5,
):
    """Classify a contig from VirNucPro chunk-level highest-score rows."""
    group = group.copy()
    group["delta"] = group["max_score_1"] - group["max_score_0"]
    group["confidence"] = group["delta"].abs().pow(0.5)

    if group["confidence"].sum() > 0:
        weighted_delta = (
            group["delta"] * group["confidence"]
        ).sum() / group["confidence"].sum()
    else:
        weighted_delta = group["delta"].mean()

    confident_viral = (group["max_score_1"] > 0.8) & (group["max_score_0"] < 0.3)
    confident_nonviral = (group["max_score_0"] > 0.8) & (group["max_score_1"] < 0.3)
    ambiguous = (group["max_score_1"] > 0.7) & (group["max_score_0"] > 0.7)

    n_chunks = len(group)
    n_confident_viral = int(confident_viral.sum())
    n_confident_nonviral = int(confident_nonviral.sum())
    n_ambiguous = int(ambiguous.sum())
    n_effective = n_chunks - n_ambiguous
    viral_proportion = n_confident_viral / n_effective if n_effective > 0 else 0
    nonviral_proportion = n_confident_nonviral / n_effective if n_effective > 0 else 0

    if weighted_delta > 0.3:
        call = "Viral"
        if n_confident_viral >= 1 and viral_proportion >= min_viral_proportion:
            tier = "high_confidence" if weighted_delta > 0.6 else "moderate_confidence"
        else:
            tier = "low_confidence"
    elif weighted_delta < -0.3:
        call = "Non-viral"
        if (
            n_confident_nonviral >= 1
            and nonviral_proportion >= min_nonviral_proportion
        ):
            tier = "high_confidence" if weighted_delta < -0.6 else "moderate_confidence"
        else:
            tier = "low_confidence"
    else:
        call = "Ambiguous"
        tier = "review"

    if n_chunks < min_chunk_count:
        if tier in ("high_confidence", "moderate_confidence"):
            tier = "low_confidence"
        elif tier == "low_confidence":
            tier = "review"

    return {
        "call": call,
        "tier": tier,
        "weighted_delta": round(weighted_delta, 3),
        "n_chunks": n_chunks,
        "n_confident_viral": n_confident_viral,
        "n_confident_nonviral": n_confident_nonviral,
        "n_ambiguous": n_ambiguous,
        "viral_proportion": round(viral_proportion, 3),
        "nonviral_proportion": round(nonviral_proportion, 3),
    }


def classify_contigs(
    highestscore_tsv,
    output_tsv,
    min_viral_prop=0.1,
    min_nonviral_prop=0.1,
    min_chunks=5,
    id_col="Modified_ID",
    id_pattern=r"(NODE\_\d+)",
):
    """Classify contigs from a VirNucPro highest-score TSV."""
    df = pd.read_csv(highestscore_tsv, sep="\t")
    required_cols = [id_col, "max_score_0", "max_score_1"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError("Missing required columns: {}".format(", ".join(missing)))

    ids = df[id_col].astype(str)
    df["ID"] = ids.str.replace(r"_chunk_\d+$", "", regex=True)
    df["_group"] = ids.str.extract(id_pattern, expand=False)

    n_unmatched = int(df["_group"].isna().sum())
    if n_unmatched == len(df):
        raise ValueError("No valid IDs extracted. Check id_pattern.")
    if n_unmatched > 0:
        log.warning(
            "%s of %s rows did not match id_pattern and were excluded",
            n_unmatched,
            len(df),
        )

    df = df.dropna(subset=["_group"])
    rows = []
    for _, group in df.groupby("_group", sort=False):
        row = _classify_contig_group(
            group,
            min_viral_proportion=min_viral_prop,
            min_nonviral_proportion=min_nonviral_prop,
            min_chunk_count=min_chunks,
        )
        row["ID"] = group["ID"].iloc[0]
        rows.append(row)

    results = pd.DataFrame(rows)
    columns = ["ID"] + [col for col in results.columns if col != "ID"]
    results = results[columns]
    results = results.sort_values("ID", key=lambda col: col.map(_natural_sort_key))

    _ensure_parent_dir(output_tsv)
    results.to_csv(output_tsv, sep="\t", index=False)


def _ensure_parent_dir(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        file.mkdir_p(parent)


def _natural_sort_key(value):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]
