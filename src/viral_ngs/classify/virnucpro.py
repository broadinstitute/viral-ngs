"""VirNucPro post-processing helpers."""

import csv
import logging
import os
import re
import shutil
import tempfile

import pandas as pd
import pysam

from viral_ngs.core import file

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
    _validate_percent_threshold("min_identity", min_identity)
    _validate_percent_threshold("min_query_cov", min_query_cov)

    duckdb = _import_duckdb()

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
        ) = _prepare_augmented_bam_file(aligned_bam, tmp_dir)
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

        con = _connect_duckdb(
            duckdb,
            duckdb_temp_dir=duckdb_temp_dir,
            duckdb_memory_limit=duckdb_memory_limit,
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


def _validate_percent_threshold(name, value):
    if 0 < value < 1.0:
        raise ValueError(
            "{} must be specified as a percent, not a fraction. "
            "Use {} instead of {}.".format(name, value * 100.0, value)
        )
    if value < 0 or value > 100:
        raise ValueError("{} must be between 0 and 100.".format(name))


def _import_duckdb():
    try:
        import duckdb
    except ImportError as exc:
        raise ImportError(
            "DuckDB is required for VirNucPro read-by-contig classification. "
            "Install the classify image dependencies before running this helper."
        ) from exc
    return duckdb


def _connect_duckdb(duckdb, duckdb_temp_dir, duckdb_memory_limit=None):
    """Open an in-memory DuckDB connection with a sensible memory cap.

    If `duckdb_memory_limit` is None, auto-detect from the process cgroup
    and cap DuckDB to ~75% of the available limit. Caller-supplied values
    always win. Pass an empty string ("") to explicitly opt out of any
    limit (let DuckDB use whatever it wants -- not recommended in
    container environments).
    """
    config = {"temp_directory": duckdb_temp_dir}
    if duckdb_memory_limit is None:
        duckdb_memory_limit = _default_memory_limit()
    if duckdb_memory_limit:
        config["memory_limit"] = str(duckdb_memory_limit)
        log.debug("DuckDB memory_limit set to %s", config["memory_limit"])
    return duckdb.connect(database=":memory:", config=config)


# Default cgroup file paths. Promoted to module-level constants so tests
# can monkeypatch them without having to mock os.path.isfile / open
# globally.
_CGROUP_V2_PATH = "/sys/fs/cgroup/memory.max"
_CGROUP_V1_PATH = "/sys/fs/cgroup/memory/memory.limit_in_bytes"

# cgroup v1 "no limit" sentinel: PAGE_COUNTER_MAX rounded to page size.
# Any value at or above this is treated as unlimited.
_CGROUP_V1_UNLIMITED = 0x7FFFFFFFFFFFF000

# Fraction of the cgroup limit handed to DuckDB. The remaining 25% is
# reserved for the python process, pandas dataframes, file buffers, and
# the kernel page cache for the alignment/classification spills.
_DUCKDB_MEM_FRACTION = 0.75


def _default_memory_limit():
    """Return a DuckDB-format memory cap string derived from the cgroup, or None.

    Reads cgroup v2 first, falling back to cgroup v1. Returns None when
    no cgroup limit is in effect (development hosts, dsub jobs with no
    memory cap, etc.) so the caller falls back to DuckDB's own default.
    Never raises -- on any I/O or parse error returns None so the
    wrapper is safe to call from any environment.
    """
    try:
        if os.path.isfile(_CGROUP_V2_PATH):
            with open(_CGROUP_V2_PATH, "rt") as fh:
                raw = fh.read().strip()
            if raw == "max":
                return None
            total_bytes = int(raw)
        elif os.path.isfile(_CGROUP_V1_PATH):
            with open(_CGROUP_V1_PATH, "rt") as fh:
                total_bytes = int(fh.read().strip())
            if total_bytes >= _CGROUP_V1_UNLIMITED:
                return None
        else:
            return None
    except (OSError, ValueError):
        return None

    if total_bytes <= 0:
        return None

    capped_bytes = int(total_bytes * _DUCKDB_MEM_FRACTION)
    # Express in MiB rounded down -- DuckDB accepts integer-MB strings
    # cleanly and avoids float-precision artifacts in the limit string.
    capped_mib = capped_bytes // (1024 * 1024)
    if capped_mib <= 0:
        return None
    return "{}MB".format(capped_mib)


_CIGAR_ALIGNMENT_BLOCK_OPS = {
    pysam.CMATCH,
    pysam.CINS,
    pysam.CDEL,
    pysam.CEQUAL,
    pysam.CDIFF,
}


def _prepare_augmented_bam_file(aligned_bam, tmp_dir):
    normalized_alignments = os.path.join(tmp_dir, "normalized.alignments.tsv")
    n_written = 0
    n_secondary = 0

    with pysam.AlignmentFile(aligned_bam, "rb") as bam, open(
        normalized_alignments, "wt"
    ) as outf:
        for record in bam:
            if record.is_unmapped:
                continue
            if record.is_secondary or record.is_supplementary:
                n_secondary += 1
                continue

            n_written += 1
            row = _bam_record_to_alignment_row(bam, record, n_written)
            outf.write("\t".join(row) + "\n")

    return normalized_alignments, n_written, n_secondary


def _bam_record_to_alignment_row(bam, record, source_order):
    if not record.has_tag("NM"):
        raise ValueError(
            "Mapped BAM record {} is missing required NM tag".format(
                record.query_name
            )
        )

    query_length = record.infer_read_length()
    if query_length is None:
        query_length = record.infer_query_length(always=True)
    if query_length is None:
        query_length = record.query_length
    if query_length is None or query_length == 0:
        raise ValueError(
            "Mapped BAM record {} has no inferable query length".format(
                record.query_name
            )
        )

    query_start = record.query_alignment_start
    query_end = record.query_alignment_end
    target_name = bam.get_reference_name(record.reference_id)
    target_length = bam.get_reference_length(target_name)
    target_start = record.reference_start
    target_end = record.reference_end
    if target_end is None:
        raise ValueError(
            "Mapped BAM record {} has no inferable reference end".format(
                record.query_name
            )
        )

    alignment_block_length = _alignment_block_length(record)
    if alignment_block_length == 0:
        raise ValueError(
            "Mapped BAM record {} has zero alignment block length".format(
                record.query_name
            )
        )

    nm = int(record.get_tag("NM"))
    num_matches = alignment_block_length - nm
    if num_matches < 0:
        raise ValueError(
            "Mapped BAM record {} has NM larger than alignment block length".format(
                record.query_name
            )
        )

    pct_identity = 100.0 * num_matches / alignment_block_length
    pct_query_cov = 100.0 * (query_end - query_start) / query_length

    return [
        str(source_order),
        record.query_name,
        str(query_length),
        str(query_start),
        str(query_end),
        "-" if record.is_reverse else "+",
        target_name,
        str(target_length),
        str(target_start),
        str(target_end),
        str(num_matches),
        str(alignment_block_length),
        str(record.mapping_quality),
        "{:.6f}".format(pct_identity),
        "{:.6f}".format(pct_query_cov),
    ]


def _alignment_block_length(record):
    return sum(
        length
        for op, length in (record.cigartuples or [])
        if op in _CIGAR_ALIGNMENT_BLOCK_OPS
    )


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
    con.execute(
        """
        CREATE TABLE alignment_raw AS
        SELECT * FROM read_csv(
            ?,
            delim='\t',
            header=false,
            all_varchar=true,
            columns={
                'source_order': 'VARCHAR',
                'query_name': 'VARCHAR',
                'query_length': 'VARCHAR',
                'query_start': 'VARCHAR',
                'query_end': 'VARCHAR',
                'strand': 'VARCHAR',
                'target_name': 'VARCHAR',
                'target_length': 'VARCHAR',
                'target_start': 'VARCHAR',
                'target_end': 'VARCHAR',
                'num_matches': 'VARCHAR',
                'alignment_block_length': 'VARCHAR',
                'mapping_quality': 'VARCHAR',
                'pct_identity': 'VARCHAR',
                'pct_query_cov': 'VARCHAR'
            }
        )
        """,
        [normalized_alignments],
    )

    con.execute(
        """
        CREATE TABLE alignment_typed AS
        SELECT
            CAST(source_order AS BIGINT) AS source_order,
            query_name,
            CAST(query_length AS BIGINT) AS query_length,
            CAST(query_start AS BIGINT) AS query_start,
            CAST(query_end AS BIGINT) AS query_end,
            strand,
            target_name,
            CAST(target_length AS BIGINT) AS target_length,
            CAST(target_start AS BIGINT) AS target_start,
            CAST(target_end AS BIGINT) AS target_end,
            CAST(num_matches AS BIGINT) AS num_matches,
            CAST(alignment_block_length AS BIGINT) AS alignment_block_length,
            CAST(mapping_quality AS BIGINT) AS mapping_quality,
            CAST(pct_identity AS DOUBLE) AS pct_identity,
            CAST(pct_query_cov AS DOUBLE) AS pct_query_cov
        FROM alignment_raw
        """
    )
    con.execute("DROP TABLE alignment_raw")

    con.execute(
        """
        CREATE TABLE alignment_filtered AS
        SELECT
            query_name,
            query_length,
            query_start,
            query_end,
            strand,
            target_name,
            target_length,
            target_start,
            target_end,
            num_matches,
            alignment_block_length,
            mapping_quality,
            pct_identity,
            pct_query_cov,
            source_order,
            (
                mapping_quality >= ?
                AND pct_identity >= ?
                AND pct_query_cov >= ?
            ) AS mapped_well
        FROM alignment_typed
        """,
        [min_mapq, min_identity, min_query_cov],
    )
    con.execute("DROP TABLE alignment_typed")
    n_primary = con.execute("SELECT count(*) FROM alignment_filtered").fetchone()[0]

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
