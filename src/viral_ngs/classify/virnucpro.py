"""VirNucPro post-processing helpers."""

import csv
import logging
import os
import re
import shutil
import tempfile

import pandas as pd

import viral_ngs.core.duckdb_utils as duckdb_utils
from viral_ngs.core import file
from . import contig_alignments

log = logging.getLogger(__name__)

READ_CLASSIFICATION_COLUMNS = [
    "read_id",
    "read_length",
    "contig_id",
    "contig_length",
    "strand",
    "mapping_quality",
    "pct_identity",
    "pct_query_cov",
    "mapped_well",
    "call",
    "tier",
    "weighted_delta",
    "n_chunks",
    "n_confident_viral",
    "n_confident_nonviral",
    "n_ambiguous",
    "viral_proportion",
    "nonviral_proportion",
]

CONTIG_CLASSIFICATION_COLUMNS = [
    "ID",
    "call",
    "tier",
    "weighted_delta",
    "n_chunks",
    "n_confident_viral",
    "n_confident_nonviral",
    "n_ambiguous",
    "viral_proportion",
    "nonviral_proportion",
]


def _classify_contig_group(
    group,
    min_viral_proportion=0.1,
    min_nonviral_proportion=0.1,
    min_chunk_count=5,
    min_confident_score=0.8,
    max_opposing_score=0.3,
    min_ambiguous_score=0.7,
    min_weighted_delta=0.3,
    high_confidence_delta=0.6,
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

    confident_viral = (
        (group["max_score_1"] > min_confident_score)
        & (group["max_score_0"] < max_opposing_score)
    )
    confident_nonviral = (
        (group["max_score_0"] > min_confident_score)
        & (group["max_score_1"] < max_opposing_score)
    )
    ambiguous = (
        (group["max_score_1"] > min_ambiguous_score)
        & (group["max_score_0"] > min_ambiguous_score)
    )

    n_chunks = len(group)
    n_confident_viral = int(confident_viral.sum())
    n_confident_nonviral = int(confident_nonviral.sum())
    n_ambiguous = int(ambiguous.sum())
    n_effective = n_chunks - n_ambiguous
    viral_proportion = n_confident_viral / n_effective if n_effective > 0 else 0
    nonviral_proportion = n_confident_nonviral / n_effective if n_effective > 0 else 0

    if weighted_delta > min_weighted_delta:
        call = "Viral"
        if n_confident_viral >= 1 and viral_proportion >= min_viral_proportion:
            tier = (
                "high_confidence"
                if weighted_delta > high_confidence_delta
                else "moderate_confidence"
            )
        else:
            tier = "low_confidence"
    elif weighted_delta < -min_weighted_delta:
        call = "Non-viral"
        if (
            n_confident_nonviral >= 1
            and nonviral_proportion >= min_nonviral_proportion
        ):
            tier = (
                "high_confidence"
                if weighted_delta < -high_confidence_delta
                else "moderate_confidence"
            )
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
    min_confident_score=0.8,
    max_opposing_score=0.3,
    min_ambiguous_score=0.7,
    min_weighted_delta=0.3,
    high_confidence_delta=0.6,
    id_col="Modified_ID",
    id_pattern=r"(NODE\_\d+)",
):
    """Classify contigs from a VirNucPro highest-score TSV."""
    try:
        df = pd.read_csv(highestscore_tsv, sep="\t")
    except pd.errors.EmptyDataError:
        _write_empty_contig_classifications(output_tsv)
        return

    required_cols = [id_col, "max_score_0", "max_score_1"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError("Missing required columns: {}".format(", ".join(missing)))

    if df.empty:
        _write_empty_contig_classifications(output_tsv)
        return

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
            min_confident_score=min_confident_score,
            max_opposing_score=max_opposing_score,
            min_ambiguous_score=min_ambiguous_score,
            min_weighted_delta=min_weighted_delta,
            high_confidence_delta=high_confidence_delta,
        )
        row["ID"] = group["ID"].iloc[0]
        rows.append(row)

    results = pd.DataFrame(rows)
    columns = ["ID"] + [col for col in results.columns if col != "ID"]
    results = results[columns]
    results = results.sort_values("ID", key=lambda col: col.map(_natural_sort_key))

    _ensure_parent_dir(output_tsv)
    results.to_csv(output_tsv, sep="\t", index=False)


def classify_reads_by_contig(
    aligned_bam,
    contig_classifications,
    output_tsv,
    min_mapq=5,
    min_identity=90.0,
    min_query_cov=80.0,
    duckdb_memory_limit=None,
    work_dir=None,
):
    """Classify reads using minimap2 BAM alignments and VirNucPro contig calls."""
    if not os.path.isfile(aligned_bam):
        raise FileNotFoundError(aligned_bam)
    if not os.path.isfile(contig_classifications):
        raise FileNotFoundError(contig_classifications)

    min_mapq = int(min_mapq)
    min_identity = float(min_identity)
    min_query_cov = float(min_query_cov)
    contig_alignments.validate_percent_threshold("min_identity", min_identity)
    contig_alignments.validate_percent_threshold("min_query_cov", min_query_cov)

    duckdb = duckdb_utils.import_duckdb("VirNucPro read-by-contig classification")

    _ensure_parent_dir(output_tsv)
    if work_dir:
        file.mkdir_p(work_dir)

    with tempfile.TemporaryDirectory(
        prefix="virnucpro_reads_", dir=work_dir
    ) as tmp_dir:
        (
            normalized_alignments,
            n_alignments,
            n_secondary,
        ) = contig_alignments.normalize_primary_alignments(aligned_bam, tmp_dir)
        normalized_classifications = _prepare_contig_classifications_file(
            contig_classifications,
            tmp_dir,
        )
        if n_alignments == 0:
            _write_empty_read_classifications(output_tsv)
            stats = {
                "n_primary": 0,
                "n_secondary": n_secondary,
                "n_reads": 0,
                "n_well": 0,
                "n_multi": 0,
                "pct_well": 0,
            }
            if stats["n_secondary"] > 0:
                log.info(
                    "Removed %s secondary/supplementary alignments",
                    stats["n_secondary"],
                )
            log.info("%s primary alignments retained", stats["n_primary"])
            log.info("%s reads in output", stats["n_reads"])
            log.info(
                "%s reads mapped well (%.1f%%)",
                stats["n_well"],
                stats["pct_well"],
            )
            return

        duckdb_temp_dir = os.path.join(tmp_dir, "duckdb_spill")
        file.mkdir_p(duckdb_temp_dir)

        con = duckdb_utils.connect(
            duckdb_temp_dir,
            memory_limit=duckdb_memory_limit,
            duckdb_module=duckdb,
        )
        try:
            result_tsv = os.path.join(tmp_dir, "reads_classified.tsv")
            stats = _classify_reads_by_contig_duckdb(
                con,
                normalized_alignments,
                normalized_classifications,
                result_tsv,
                min_mapq=min_mapq,
                min_identity=min_identity,
                min_query_cov=min_query_cov,
                n_secondary=n_secondary,
            )
        finally:
            con.close()

        with open(result_tsv, "rt") as inf, file.open_or_gzopen(
            output_tsv, "wt"
        ) as outf:
            shutil.copyfileobj(inf, outf)

    if stats["n_secondary"] > 0:
        log.info(
            "Removed %s secondary/supplementary alignments",
            stats["n_secondary"],
        )
    log.info("%s primary alignments retained", stats["n_primary"])
    log.info("%s reads in output", stats["n_reads"])
    log.info("%s reads mapped well (%.1f%%)", stats["n_well"], stats["pct_well"])
    if stats["n_multi"] > 0:
        log.info("%s reads flagged as Multi-mapped", stats["n_multi"])


def _prepare_contig_classifications_file(contig_classifications, tmp_dir):
    normalized_classifications = os.path.join(
        tmp_dir,
        "contig_classifications.tsv",
    )
    with file.open_or_gzopen(contig_classifications, "rt") as inf, open(
        normalized_classifications, "wt"
    ) as outf:
        reader = csv.DictReader(inf, delimiter="\t")
        missing = [
            col
            for col in CONTIG_CLASSIFICATION_COLUMNS
            if col not in (reader.fieldnames or [])
        ]
        if missing:
            raise ValueError(
                "Missing required columns: {}".format(", ".join(missing))
            )

        writer = csv.DictWriter(
            outf,
            fieldnames=CONTIG_CLASSIFICATION_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in reader:
            writer.writerow(row)

    return normalized_classifications


def _classify_reads_by_contig_duckdb(
    con,
    normalized_alignments,
    normalized_classifications,
    result_tsv,
    min_mapq,
    min_identity,
    min_query_cov,
    n_secondary,
):
    n_primary = contig_alignments.load_primary_alignments(
        con, normalized_alignments, min_mapq, min_identity, min_query_cov
    )

    con.execute(
        """
        CREATE TABLE clf AS
        SELECT
            ID,
            call,
            tier,
            CAST(weighted_delta AS DOUBLE) AS weighted_delta,
            CAST(n_chunks AS BIGINT) AS n_chunks,
            CAST(n_confident_viral AS BIGINT) AS n_confident_viral,
            CAST(n_confident_nonviral AS BIGINT) AS n_confident_nonviral,
            CAST(n_ambiguous AS BIGINT) AS n_ambiguous,
            CAST(viral_proportion AS DOUBLE) AS viral_proportion,
            CAST(nonviral_proportion AS DOUBLE) AS nonviral_proportion
        FROM read_csv(
            ?,
            delim='\t',
            header=true,
            all_varchar=true
        )
        """,
        [normalized_classifications],
    )

    con.execute(
        """
        CREATE TABLE merged AS
        SELECT
            p.query_name,
            p.query_length,
            p.target_name,
            p.target_length,
            p.strand,
            p.mapping_quality,
            p.pct_identity,
            p.pct_query_cov,
            p.mapped_well,
            p.source_order,
            COALESCE(c.call, 'Unclassified') AS call,
            COALESCE(c.tier, '') AS tier,
            COALESCE(c.weighted_delta, 0) AS weighted_delta,
            COALESCE(c.n_chunks, 0) AS n_chunks,
            COALESCE(c.n_confident_viral, 0) AS n_confident_viral,
            COALESCE(c.n_confident_nonviral, 0) AS n_confident_nonviral,
            COALESCE(c.n_ambiguous, 0) AS n_ambiguous,
            COALESCE(c.viral_proportion, 0) AS viral_proportion,
            COALESCE(c.nonviral_proportion, 0) AS nonviral_proportion
        FROM alignment_filtered p
        LEFT JOIN clf c ON p.target_name = c.ID
        """
    )
    con.execute("DROP TABLE alignment_filtered")
    con.execute("DROP TABLE clf")

    con.execute(
        """
        CREATE TABLE call_counts AS
        SELECT query_name, COUNT(DISTINCT call) AS n_distinct_calls
        FROM merged
        GROUP BY query_name
        """
    )

    con.execute(
        """
        CREATE TABLE best AS
        SELECT * EXCLUDE (rn) FROM (
            SELECT *,
                ROW_NUMBER() OVER (
                    PARTITION BY query_name
                    ORDER BY mapping_quality DESC, pct_identity DESC, source_order ASC
                ) AS rn
            FROM merged
        ) WHERE rn = 1
        """
    )
    con.execute("DROP TABLE merged")

    con.execute(
        """
        CREATE TABLE result AS
        SELECT
            b.query_name AS read_id,
            b.query_length AS read_length,
            b.target_name AS contig_id,
            b.target_length AS contig_length,
            b.strand,
            b.mapping_quality,
            b.pct_identity,
            b.pct_query_cov,
            b.mapped_well,
            CASE
                WHEN cc.n_distinct_calls > 1 THEN 'Multi-mapped'
                ELSE b.call
            END AS call,
            CASE
                WHEN cc.n_distinct_calls > 1 THEN 'review'
                ELSE b.tier
            END AS tier,
            b.weighted_delta,
            b.n_chunks,
            b.n_confident_viral,
            b.n_confident_nonviral,
            b.n_ambiguous,
            b.viral_proportion,
            b.nonviral_proportion
        FROM best b
        JOIN call_counts cc ON b.query_name = cc.query_name
        """
    )
    con.execute("DROP TABLE best")
    con.execute("DROP TABLE call_counts")

    n_reads = con.execute("SELECT count(*) FROM result").fetchone()[0]
    n_well = con.execute("SELECT count(*) FROM result WHERE mapped_well").fetchone()[0]
    n_multi = con.execute(
        "SELECT count(*) FROM result WHERE call = 'Multi-mapped'"
    ).fetchone()[0]
    pct_well = (n_well / n_reads * 100) if n_reads > 0 else 0

    con.execute(
        """
        COPY (
            SELECT
                read_id,
                read_length,
                contig_id,
                contig_length,
                strand,
                mapping_quality,
                pct_identity,
                pct_query_cov,
                mapped_well,
                call,
                tier,
                weighted_delta,
                n_chunks,
                n_confident_viral,
                n_confident_nonviral,
                n_ambiguous,
                viral_proportion,
                nonviral_proportion
            FROM result
        )
        TO ? (FORMAT CSV, DELIMITER '\t', HEADER TRUE)
        """,
        [result_tsv],
    )

    return {
        "n_primary": n_primary,
        "n_secondary": n_secondary,
        "n_reads": n_reads,
        "n_well": n_well,
        "n_multi": n_multi,
        "pct_well": pct_well,
    }


def _write_empty_read_classifications(output_tsv):
    with file.open_or_gzopen(output_tsv, "wt") as outf:
        writer = csv.writer(outf, delimiter="\t", lineterminator="\n")
        writer.writerow(READ_CLASSIFICATION_COLUMNS)


def _write_empty_contig_classifications(output_tsv):
    _ensure_parent_dir(output_tsv)
    pd.DataFrame(columns=CONTIG_CLASSIFICATION_COLUMNS).to_csv(
        output_tsv,
        sep="\t",
        index=False,
    )


def _ensure_parent_dir(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        file.mkdir_p(parent)


def _natural_sort_key(value):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]
