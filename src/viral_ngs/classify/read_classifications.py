"""Read-level classification annotation and joining helpers."""

import csv
import logging
import os
import tempfile

import viral_ngs.core.duckdb_utils as duckdb_utils
from viral_ngs.core import file

log = logging.getLogger(__name__)

KRAKEN2_ANNOTATED_COLUMNS = [
    "SAMPLE_ID",
    "READ_ID",
    "TAXONOMY_ID",
    "TAX_NAME",
    "KINGDOM",
    "TAX_RANK",
]

CENTRIFUGER_ANNOTATED_COLUMNS = [
    "SAMPLE_ID",
    "READ_ID",
    "SEQ_ID",
    "TAXONOMY_ID",
    "SCORE",
    "SECOND_BEST_SCORE",
    "HIT_LENGTH",
    "QUERY_LENGTH",
    "NUM_MATCHES",
    "TAX_NAME",
    "KINGDOM",
    "TAX_RANK",
]

KALLISTO_REQUIRED_COLUMNS = [
    "SAMPLE_ID",
    "READ_ID",
    "DB_ID",
    "TAXONOMY_LINEAGE",
    "TAXONOMY_NAME",
    "SEQUENCE_LENGTH",
]

VNP_REQUIRED_COLUMNS = [
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

GENOMAD_COLUMNS = {
    "seq_name": None,
    "length": "GENOMAD_LENGTH",
    "topology": "GENOMAD_TOPOLOGY",
    "coordinates": "GENOMAD_COORDINATES",
    "n_genes": "GENOMAD_N_GENES",
    "genetic_code": "GENOMAD_GENETIC_CODE",
    "virus_score": "GENOMAD_VIRUS_SCORE",
    "fdr": "GENOMAD_FDR",
    "n_hallmarks": "GENOMAD_N_HALLMARKS",
    "marker_enrichment": "GENOMAD_MARKER_ENRICHMENT",
    "taxonomy": "GENOMAD_TAXONOMY",
}

FINAL_SCHEMA = [
    ("SAMPLE_ID", "VARCHAR"),
    ("READ_ID", "VARCHAR"),
    ("KALLISTO_DB_ID", "VARCHAR"),
    ("KALLISTO_TAXONOMY_LINEAGE", "VARCHAR"),
    ("KALLISTO_TAXONOMY_NAME", "VARCHAR"),
    ("KALLISTO_SEQUENCE_LENGTH", "BIGINT"),
    ("KALLISTO_CONCORDANCE", "VARCHAR"),
    ("K2_TAXONOMY_ID", "BIGINT"),
    ("K2_TAX_NAME", "VARCHAR"),
    ("K2_KINGDOM", "VARCHAR"),
    ("VNP_CONTIG_ID", "VARCHAR"),
    ("VNP_CONTIG_LENGTH", "BIGINT"),
    ("VNP_STRAND", "VARCHAR"),
    ("VNP_READ_LENGTH", "BIGINT"),
    ("VNP_MAPPING_QUALITY", "BIGINT"),
    ("VNP_PCT_IDENTITY", "DOUBLE"),
    ("VNP_PCT_QUERY_COV", "DOUBLE"),
    ("VNP_MAPPED_WELL", "BOOLEAN"),
    ("VNP_CALL", "VARCHAR"),
    ("VNP_TIER", "VARCHAR"),
    ("VNP_WEIGHTED_DELTA", "DOUBLE"),
    ("VNP_N_CHUNKS", "BIGINT"),
    ("VNP_N_CONFIDENT_VIRAL", "BIGINT"),
    ("VNP_N_CONFIDENT_NONVIRAL", "BIGINT"),
    ("VNP_N_AMBIGUOUS", "BIGINT"),
    ("VNP_VIRAL_PROPORTION", "DOUBLE"),
    ("VNP_NONVIRAL_PROPORTION", "DOUBLE"),
    ("VNP_CONCORDANCE", "VARCHAR"),
    ("CENTRIFUGER_SEQ_ID", "VARCHAR"),
    ("CENTRIFUGER_TAXONOMY_ID", "BIGINT"),
    ("CENTRIFUGER_SCORE", "BIGINT"),
    ("CENTRIFUGER_SECOND_BEST_SCORE", "BIGINT"),
    ("CENTRIFUGER_HIT_LENGTH", "BIGINT"),
    ("CENTRIFUGER_QUERY_LENGTH", "BIGINT"),
    ("CENTRIFUGER_NUM_MATCHES", "BIGINT"),
    ("CENTRIFUGER_TAX_NAME", "VARCHAR"),
    ("CENTRIFUGER_KINGDOM", "VARCHAR"),
    ("CENTRIFUGER_TAX_RANK", "VARCHAR"),
    ("GENOMAD_LENGTH", "BIGINT"),
    ("GENOMAD_TOPOLOGY", "VARCHAR"),
    ("GENOMAD_COORDINATES", "VARCHAR"),
    ("GENOMAD_N_GENES", "BIGINT"),
    ("GENOMAD_GENETIC_CODE", "BIGINT"),
    ("GENOMAD_VIRUS_SCORE", "DOUBLE"),
    ("GENOMAD_FDR", "DOUBLE"),
    ("GENOMAD_N_HALLMARKS", "BIGINT"),
    ("GENOMAD_MARKER_ENRICHMENT", "DOUBLE"),
    ("GENOMAD_TAXONOMY", "VARCHAR"),
]

FINAL_COLUMNS = [name for name, _ in FINAL_SCHEMA]

CENTRIFUGER_NATIVE_COLUMNS = [
    "readID",
    "seqID",
    "taxID",
    "score",
    "2ndBestScore",
    "hitLength",
    "queryLength",
    "numMatches",
]


class DuckDBTaxonomyDatabase:
    """Read taxonomy annotations from a distributed DuckDB reference DB."""

    def __init__(self, taxonomy_duckdb):
        if not os.path.isfile(taxonomy_duckdb):
            raise FileNotFoundError(taxonomy_duckdb)
        duckdb = duckdb_utils.import_duckdb("taxonomy annotation")
        self._con = duckdb.connect(database=taxonomy_duckdb, read_only=True)
        self._columns = self._table_columns("taxonomy")
        required = {"taxid", "name", "kingdom", "rank"}
        missing = sorted(required - self._columns)
        if missing:
            raise ValueError(
                "taxonomy table missing required columns: {}".format(
                    ", ".join(missing)
                )
            )
        self._cache = {}

    def close(self):
        self._con.close()

    def lookup(self, taxid, resolve_strains=False):
        taxid = _parse_int(taxid, "taxid")
        cache_key = (taxid, bool(resolve_strains))
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        if taxid == 0:
            result = {
                "taxid": 0,
                "name": "Unclassified",
                "kingdom": "Unclassified",
                "rank": "unclassified",
            }
        else:
            row = self._fetch_taxon(taxid)
            if row is None:
                result = {
                    "taxid": taxid,
                    "name": "Unknown",
                    "kingdom": "Unknown",
                    "rank": "unknown",
                }
            else:
                if resolve_strains and "parent_id" in self._columns:
                    row = self._resolve_strain(row)
                result = {
                    "taxid": row["taxid"],
                    "name": row["name"],
                    "kingdom": row["kingdom"],
                    "rank": row["rank"] or "unknown",
                }

        self._cache[cache_key] = result
        return result

    def _table_columns(self, table_name):
        try:
            rows = self._con.execute("PRAGMA table_info('{}')".format(table_name)).fetchall()
        except Exception as exc:
            raise ValueError("taxonomy table is missing") from exc
        columns = {row[1] for row in rows}
        if not columns:
            raise ValueError("taxonomy table is missing")
        return columns

    def _fetch_taxon(self, taxid):
        select_cols = ["taxid", "name", "kingdom", "rank"]
        if "parent_id" in self._columns:
            select_cols.append("parent_id")
        row = self._con.execute(
            "SELECT {} FROM taxonomy WHERE taxid = ?".format(
                ", ".join(select_cols)
            ),
            [taxid],
        ).fetchone()
        if row is None:
            return None
        return dict(zip(select_cols, row))

    def _resolve_strain(self, row):
        seen = set()
        current = row
        while (
            current
            and current.get("rank") == "no rank"
            and current.get("parent_id") is not None
            and current["taxid"] not in seen
        ):
            seen.add(current["taxid"])
            parent = self._fetch_taxon(current["parent_id"])
            if parent is None:
                break
            current = parent
            if current.get("rank") and current.get("rank") != "no rank":
                return current
        return row


def kraken2_annotate_reads(
    kraken2_reads,
    taxonomy_duckdb,
    output_tsv,
    sample_id,
    resolve_strains=False,
):
    """Annotate native Kraken2 per-read output with taxonomy metadata."""
    if not os.path.isfile(kraken2_reads):
        raise FileNotFoundError(kraken2_reads)

    taxonomy = DuckDBTaxonomyDatabase(taxonomy_duckdb)
    try:
        _ensure_parent_dir(output_tsv)
        with file.open_or_gzopen(output_tsv, "wt") as outf:
            writer = csv.DictWriter(
                outf,
                KRAKEN2_ANNOTATED_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            with file.open_or_gzopen(kraken2_reads, "rt") as inf:
                for line_no, line in enumerate(inf, start=1):
                    if not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) != 5:
                        raise ValueError(
                            "Malformed Kraken2 row at line {}: expected 5 columns".format(
                                line_no
                            )
                        )
                    status, read_id, taxid, _, _ = cols
                    if status not in ("C", "U"):
                        raise ValueError(
                            "Malformed Kraken2 row at line {}: first column must be C or U".format(
                                line_no
                            )
                        )
                    taxon = taxonomy.lookup(taxid, resolve_strains=resolve_strains)
                    writer.writerow(
                        {
                            "SAMPLE_ID": sample_id,
                            "READ_ID": read_id,
                            "TAXONOMY_ID": taxon["taxid"],
                            "TAX_NAME": taxon["name"],
                            "KINGDOM": taxon["kingdom"],
                            "TAX_RANK": taxon["rank"],
                        }
                    )
    finally:
        taxonomy.close()


def centrifuger_annotate_reads(
    centrifuger_reads,
    taxonomy_duckdb,
    output_tsv,
    sample_id,
    resolve_strains=False,
):
    """Annotate native Centrifuger per-read output with taxonomy metadata."""
    if not os.path.isfile(centrifuger_reads):
        raise FileNotFoundError(centrifuger_reads)

    taxonomy = DuckDBTaxonomyDatabase(taxonomy_duckdb)
    try:
        _ensure_parent_dir(output_tsv)
        with file.open_or_gzopen(output_tsv, "wt") as outf:
            writer = csv.DictWriter(
                outf,
                CENTRIFUGER_ANNOTATED_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            with file.open_or_gzopen(centrifuger_reads, "rt") as inf:
                reader = csv.DictReader(inf, delimiter="\t")
                if reader.fieldnames is None:
                    return
                _require_columns(
                    reader.fieldnames,
                    CENTRIFUGER_NATIVE_COLUMNS,
                    "Centrifuger reads",
                )
                for line_no, row in enumerate(reader, start=2):
                    if row is None or not any(row.values()):
                        continue
                    if None in row or any(
                        row.get(col) is None for col in CENTRIFUGER_NATIVE_COLUMNS
                    ):
                        raise ValueError(
                            "Malformed Centrifuger row at line {}: wrong column count".format(
                                line_no
                            )
                        )
                    taxon = taxonomy.lookup(
                        row["taxID"],
                        resolve_strains=resolve_strains,
                    )
                    writer.writerow(
                        {
                            "SAMPLE_ID": sample_id,
                            "READ_ID": row["readID"],
                            "SEQ_ID": row["seqID"],
                            "TAXONOMY_ID": taxon["taxid"],
                            "SCORE": row["score"],
                            "SECOND_BEST_SCORE": row["2ndBestScore"],
                            "HIT_LENGTH": row["hitLength"],
                            "QUERY_LENGTH": row["queryLength"],
                            "NUM_MATCHES": row["numMatches"],
                            "TAX_NAME": taxon["name"],
                            "KINGDOM": taxon["kingdom"],
                            "TAX_RANK": taxon["rank"],
                        }
                    )
    finally:
        taxonomy.close()


def _subset_schema(prefix):
    schema = [("SAMPLE_ID", "VARCHAR"), ("READ_ID", "VARCHAR")]
    schema += [(name, typ) for name, typ in FINAL_SCHEMA if name.startswith(prefix)]
    return schema


_GENOMAD_TABLE_SCHEMA = [("CONTIG_ID", "VARCHAR")] + [
    (name, typ) for name, typ in FINAL_SCHEMA if name.startswith("GENOMAD_")
]


def join_read_classifications(
    out_parquet,
    sample_id,
    kallisto_summary=None,
    kraken2_reads=None,
    vnp_reads=None,
    genomad_virus_summary=None,
    centrifuger_reads=None,
    filter_human_only_k2=True,
    duckdb_memory_limit=None,
    work_dir=None,
):
    """Join read-level classifier outputs into a stable Parquet table.

    The join runs entirely in DuckDB so large read sets stream from disk and
    spill to a per-run temp directory rather than being materialized in Python
    memory. Each source becomes a typed table (empty when the source is absent),
    Kallisto/VirNucPro/Kraken2 are FULL OUTER JOINed as primary read sources, and
    Centrifuger/geNomad are LEFT JOINed as enrichments.
    """
    if not any(
        _source_exists(path)
        for path in (
            kallisto_summary,
            vnp_reads,
            kraken2_reads,
            centrifuger_reads,
            genomad_virus_summary,
        )
    ):
        raise ValueError("At least one classification input must exist")

    _ensure_parent_dir(out_parquet)
    if work_dir:
        file.mkdir_p(work_dir)
    with tempfile.TemporaryDirectory(
        prefix="read_classifications_", dir=work_dir
    ) as tmp_dir:
        duckdb_temp_dir = os.path.join(tmp_dir, "duckdb_spill")
        file.mkdir_p(duckdb_temp_dir)
        con = duckdb_utils.connect(duckdb_temp_dir, memory_limit=duckdb_memory_limit)
        try:
            _build_kallisto_table(con, kallisto_summary, sample_id)
            _build_vnp_table(con, vnp_reads, sample_id)
            _build_kraken2_table(con, kraken2_reads, sample_id)
            _build_centrifuger_table(con, centrifuger_reads, sample_id)
            _build_genomad_table(con, genomad_virus_summary)

            _assert_unique(con, "k2", ("SAMPLE_ID", "READ_ID"), "Kraken2 reads")
            _assert_unique(
                con, "centrifuger", ("SAMPLE_ID", "READ_ID"), "Centrifuger reads"
            )
            _assert_unique(con, "genomad", ("CONTIG_ID",), "geNomad virus summary")

            con.execute(
                "CREATE TABLE joined_kv AS "
                "SELECT * FROM kallisto FULL OUTER JOIN vnp USING (SAMPLE_ID, READ_ID)"
            )

            # Human-only Kraken2 pre-filter via a NULL-safe anti-join: drop reads
            # called Homo sapiens by Kraken2 that have no Kallisto/VNP signal.
            if filter_human_only_k2:
                con.execute(
                    "CREATE TABLE k2_kept AS "
                    "SELECT k.* FROM k2 k "
                    "LEFT JOIN (SELECT DISTINCT SAMPLE_ID, READ_ID FROM joined_kv) j "
                    "ON k.SAMPLE_ID = j.SAMPLE_ID AND k.READ_ID = j.READ_ID "
                    "WHERE NOT (k.K2_TAX_NAME = 'Homo sapiens' AND j.READ_ID IS NULL)"
                )
            else:
                con.execute("CREATE TABLE k2_kept AS SELECT * FROM k2")

            con.execute(
                "CREATE TABLE joined_all AS "
                "SELECT * FROM joined_kv "
                "FULL OUTER JOIN k2_kept USING (SAMPLE_ID, READ_ID)"
            )
            con.execute(
                "CREATE TABLE joined_enriched AS "
                "SELECT * FROM joined_all "
                "LEFT JOIN centrifuger USING (SAMPLE_ID, READ_ID)"
            )

            final_cols = ", ".join(name for name, _ in FINAL_SCHEMA)
            con.execute(
                "CREATE TABLE result AS SELECT {cols} FROM joined_enriched a "
                "LEFT JOIN genomad g ON a.VNP_CONTIG_ID = g.CONTIG_ID".format(
                    cols=final_cols
                )
            )

            total = con.execute("SELECT count(*) FROM result").fetchone()[0]
            distinct = con.execute(
                "SELECT count(*) FROM (SELECT 1 FROM result GROUP BY SAMPLE_ID, READ_ID)"
            ).fetchone()[0]
            if total != distinct:
                raise ValueError(
                    "Join produced {} rows for {} distinct (SAMPLE_ID, READ_ID) keys; "
                    "an enrichment source multiplied rows".format(total, distinct)
                )

            con.execute(
                "COPY result TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD)".format(
                    _quote_duckdb_path(out_parquet)
                )
            )
        finally:
            con.close()


def _build_kallisto_table(con, path, sample_id):
    if not _source_exists(path):
        con.execute("CREATE TABLE kallisto ({})".format(_ddl(_subset_schema("KALLISTO_"))))
        return
    con.execute(
        "CREATE TABLE kallisto_raw AS SELECT * FROM {}".format(_reader_sql(path)),
        [path],
    )
    _validate_columns(con, "kallisto_raw", KALLISTO_REQUIRED_COLUMNS, "Kallisto summary")
    con.execute(
        """
        CREATE TABLE kallisto AS
        WITH norm AS (
            SELECT
                CAST(SAMPLE_ID AS VARCHAR) AS SAMPLE_ID,
                regexp_replace(CAST(READ_ID AS VARCHAR), '/[12]$', '') AS READ_ID,
                CASE
                    WHEN CAST(READ_ID AS VARCHAR) LIKE '%/1' THEN '1'
                    WHEN CAST(READ_ID AS VARCHAR) LIKE '%/2' THEN '2'
                    ELSE ''
                END AS MATE,
                NULLIF(CAST(DB_ID AS VARCHAR), '') AS DB_ID,
                NULLIF(CAST(TAXONOMY_LINEAGE AS VARCHAR), '') AS TAXONOMY_LINEAGE,
                NULLIF(CAST(TAXONOMY_NAME AS VARCHAR), '') AS TAXONOMY_NAME,
                TRY_CAST(SEQUENCE_LENGTH AS BIGINT) AS SEQUENCE_LENGTH,
                rowid AS ord
            FROM kallisto_raw
            WHERE CAST(SAMPLE_ID AS VARCHAR) = ?
        ),
        ranked AS (
            SELECT *, ROW_NUMBER() OVER (
                PARTITION BY SAMPLE_ID, READ_ID
                ORDER BY SEQUENCE_LENGTH DESC NULLS LAST, MATE ASC, ord ASC
            ) AS rn
            FROM norm
        ),
        per_read AS (
            SELECT SAMPLE_ID, READ_ID,
                COUNT(DISTINCT NULLIF(MATE, '')) AS n_mates,
                COUNT(DISTINCT DB_ID) AS n_db
            FROM norm
            GROUP BY SAMPLE_ID, READ_ID
        )
        SELECT
            r.SAMPLE_ID,
            r.READ_ID,
            r.DB_ID AS KALLISTO_DB_ID,
            r.TAXONOMY_LINEAGE AS KALLISTO_TAXONOMY_LINEAGE,
            r.TAXONOMY_NAME AS KALLISTO_TAXONOMY_NAME,
            r.SEQUENCE_LENGTH AS KALLISTO_SEQUENCE_LENGTH,
            CASE
                WHEN p.n_mates <= 1 THEN 'Singleton'
                WHEN p.n_db <= 1 THEN 'Concordant'
                ELSE 'Discordant'
            END AS KALLISTO_CONCORDANCE
        FROM ranked r
        JOIN per_read p USING (SAMPLE_ID, READ_ID)
        WHERE r.rn = 1
        """,
        [sample_id],
    )
    con.execute("DROP TABLE kallisto_raw")


def _build_vnp_table(con, path, sample_id):
    if not _source_exists(path):
        con.execute("CREATE TABLE vnp ({})".format(_ddl(_subset_schema("VNP_"))))
        return
    con.execute(
        "CREATE TABLE vnp_raw AS SELECT * FROM {}".format(_reader_sql(path)),
        [path],
    )
    _validate_columns(con, "vnp_raw", VNP_REQUIRED_COLUMNS, "VirNucPro reads")
    con.execute(
        """
        CREATE TABLE vnp AS
        WITH norm AS (
            SELECT
                CAST(? AS VARCHAR) AS SAMPLE_ID,
                regexp_replace(CAST(read_id AS VARCHAR), '/[12]$', '') AS READ_ID,
                CASE
                    WHEN CAST(read_id AS VARCHAR) LIKE '%/1' THEN '1'
                    WHEN CAST(read_id AS VARCHAR) LIKE '%/2' THEN '2'
                    ELSE ''
                END AS MATE,
                NULLIF(CAST(contig_id AS VARCHAR), '') AS VNP_CONTIG_ID,
                TRY_CAST(contig_length AS BIGINT) AS VNP_CONTIG_LENGTH,
                NULLIF(CAST(strand AS VARCHAR), '') AS VNP_STRAND,
                TRY_CAST(read_length AS BIGINT) AS VNP_READ_LENGTH,
                TRY_CAST(mapping_quality AS BIGINT) AS VNP_MAPPING_QUALITY,
                TRY_CAST(pct_identity AS DOUBLE) AS VNP_PCT_IDENTITY,
                TRY_CAST(pct_query_cov AS DOUBLE) AS VNP_PCT_QUERY_COV,
                TRY_CAST(mapped_well AS BOOLEAN) AS VNP_MAPPED_WELL,
                NULLIF(CAST(call AS VARCHAR), '') AS VNP_CALL,
                NULLIF(CAST(tier AS VARCHAR), '') AS VNP_TIER,
                TRY_CAST(weighted_delta AS DOUBLE) AS VNP_WEIGHTED_DELTA,
                TRY_CAST(n_chunks AS BIGINT) AS VNP_N_CHUNKS,
                TRY_CAST(n_confident_viral AS BIGINT) AS VNP_N_CONFIDENT_VIRAL,
                TRY_CAST(n_confident_nonviral AS BIGINT) AS VNP_N_CONFIDENT_NONVIRAL,
                TRY_CAST(n_ambiguous AS BIGINT) AS VNP_N_AMBIGUOUS,
                TRY_CAST(viral_proportion AS DOUBLE) AS VNP_VIRAL_PROPORTION,
                TRY_CAST(nonviral_proportion AS DOUBLE) AS VNP_NONVIRAL_PROPORTION,
                rowid AS ord
            FROM vnp_raw
        ),
        ranked AS (
            SELECT *, ROW_NUMBER() OVER (
                PARTITION BY SAMPLE_ID, READ_ID
                ORDER BY VNP_MAPPING_QUALITY DESC NULLS LAST,
                         VNP_PCT_IDENTITY DESC NULLS LAST,
                         MATE ASC, ord ASC
            ) AS rn
            FROM norm
        ),
        per_read AS (
            SELECT SAMPLE_ID, READ_ID,
                COUNT(DISTINCT NULLIF(MATE, '')) AS n_mates,
                COUNT(DISTINCT VNP_CALL) AS n_calls
            FROM norm
            GROUP BY SAMPLE_ID, READ_ID
        )
        SELECT
            r.SAMPLE_ID, r.READ_ID, r.VNP_CONTIG_ID, r.VNP_CONTIG_LENGTH,
            r.VNP_STRAND, r.VNP_READ_LENGTH, r.VNP_MAPPING_QUALITY,
            r.VNP_PCT_IDENTITY, r.VNP_PCT_QUERY_COV, r.VNP_MAPPED_WELL,
            r.VNP_CALL, r.VNP_TIER, r.VNP_WEIGHTED_DELTA, r.VNP_N_CHUNKS,
            r.VNP_N_CONFIDENT_VIRAL, r.VNP_N_CONFIDENT_NONVIRAL, r.VNP_N_AMBIGUOUS,
            r.VNP_VIRAL_PROPORTION, r.VNP_NONVIRAL_PROPORTION,
            CASE
                WHEN p.n_mates <= 1 THEN 'Singleton'
                WHEN p.n_calls <= 1 THEN 'Concordant'
                ELSE 'Discordant'
            END AS VNP_CONCORDANCE
        FROM ranked r
        JOIN per_read p USING (SAMPLE_ID, READ_ID)
        WHERE r.rn = 1
        """,
        [sample_id],
    )
    con.execute("DROP TABLE vnp_raw")


def _build_kraken2_table(con, path, sample_id):
    if not _source_exists(path):
        con.execute("CREATE TABLE k2 ({})".format(_ddl(_subset_schema("K2_"))))
        return
    con.execute(
        "CREATE TABLE k2_raw AS SELECT * FROM {}".format(_reader_sql(path)),
        [path],
    )
    _validate_columns(con, "k2_raw", KRAKEN2_ANNOTATED_COLUMNS, "Kraken2 reads")
    con.execute(
        """
        CREATE TABLE k2 AS
        SELECT
            CAST(SAMPLE_ID AS VARCHAR) AS SAMPLE_ID,
            regexp_replace(CAST(READ_ID AS VARCHAR), '/[12]$', '') AS READ_ID,
            TRY_CAST(TAXONOMY_ID AS BIGINT) AS K2_TAXONOMY_ID,
            NULLIF(CAST(TAX_NAME AS VARCHAR), '') AS K2_TAX_NAME,
            NULLIF(CAST(KINGDOM AS VARCHAR), '') AS K2_KINGDOM
        FROM k2_raw
        WHERE CAST(SAMPLE_ID AS VARCHAR) = ?
        """,
        [sample_id],
    )
    con.execute("DROP TABLE k2_raw")


def _build_centrifuger_table(con, path, sample_id):
    if not _source_exists(path):
        con.execute(
            "CREATE TABLE centrifuger ({})".format(_ddl(_subset_schema("CENTRIFUGER_")))
        )
        return
    con.execute(
        "CREATE TABLE centrifuger_raw AS SELECT * FROM {}".format(_reader_sql(path)),
        [path],
    )
    _validate_columns(
        con, "centrifuger_raw", CENTRIFUGER_ANNOTATED_COLUMNS, "Centrifuger reads"
    )
    con.execute(
        """
        CREATE TABLE centrifuger AS
        SELECT
            CAST(SAMPLE_ID AS VARCHAR) AS SAMPLE_ID,
            regexp_replace(CAST(READ_ID AS VARCHAR), '/[12]$', '') AS READ_ID,
            NULLIF(CAST(SEQ_ID AS VARCHAR), '') AS CENTRIFUGER_SEQ_ID,
            TRY_CAST(TAXONOMY_ID AS BIGINT) AS CENTRIFUGER_TAXONOMY_ID,
            TRY_CAST(SCORE AS BIGINT) AS CENTRIFUGER_SCORE,
            TRY_CAST(SECOND_BEST_SCORE AS BIGINT) AS CENTRIFUGER_SECOND_BEST_SCORE,
            TRY_CAST(HIT_LENGTH AS BIGINT) AS CENTRIFUGER_HIT_LENGTH,
            TRY_CAST(QUERY_LENGTH AS BIGINT) AS CENTRIFUGER_QUERY_LENGTH,
            TRY_CAST(NUM_MATCHES AS BIGINT) AS CENTRIFUGER_NUM_MATCHES,
            NULLIF(CAST(TAX_NAME AS VARCHAR), '') AS CENTRIFUGER_TAX_NAME,
            NULLIF(CAST(KINGDOM AS VARCHAR), '') AS CENTRIFUGER_KINGDOM,
            NULLIF(CAST(TAX_RANK AS VARCHAR), '') AS CENTRIFUGER_TAX_RANK
        FROM centrifuger_raw
        WHERE CAST(SAMPLE_ID AS VARCHAR) = ?
        """,
        [sample_id],
    )
    con.execute("DROP TABLE centrifuger_raw")


def _build_genomad_table(con, path):
    if not _source_exists(path):
        con.execute("CREATE TABLE genomad ({})".format(_ddl(_GENOMAD_TABLE_SCHEMA)))
        return
    con.execute(
        "CREATE TABLE genomad_raw AS SELECT * FROM {}".format(_reader_sql(path)),
        [path],
    )
    _validate_columns(
        con, "genomad_raw", list(GENOMAD_COLUMNS.keys()), "geNomad virus summary"
    )
    con.execute(
        """
        CREATE TABLE genomad AS
        SELECT
            NULLIF(CAST(seq_name AS VARCHAR), '') AS CONTIG_ID,
            TRY_CAST(length AS BIGINT) AS GENOMAD_LENGTH,
            NULLIF(CAST(topology AS VARCHAR), '') AS GENOMAD_TOPOLOGY,
            NULLIF(CAST(coordinates AS VARCHAR), '') AS GENOMAD_COORDINATES,
            TRY_CAST(n_genes AS BIGINT) AS GENOMAD_N_GENES,
            TRY_CAST(genetic_code AS BIGINT) AS GENOMAD_GENETIC_CODE,
            TRY_CAST(virus_score AS DOUBLE) AS GENOMAD_VIRUS_SCORE,
            TRY_CAST(fdr AS DOUBLE) AS GENOMAD_FDR,
            TRY_CAST(n_hallmarks AS BIGINT) AS GENOMAD_N_HALLMARKS,
            TRY_CAST(marker_enrichment AS DOUBLE) AS GENOMAD_MARKER_ENRICHMENT,
            NULLIF(CAST(taxonomy AS VARCHAR), '') AS GENOMAD_TAXONOMY
        FROM genomad_raw
        """
    )
    con.execute("DROP TABLE genomad_raw")


def _reader_sql(path):
    if str(path).endswith(".parquet"):
        return "read_parquet(?)"
    return "read_csv(?, delim='\t', header=true, all_varchar=true, compression='auto')"


def _ddl(schema):
    return ", ".join("{} {}".format(name, typ) for name, typ in schema)


def _validate_columns(con, table, required_columns, label):
    rows = con.execute("PRAGMA table_info('{}')".format(table)).fetchall()
    columns = {row[1] for row in rows}
    missing = [col for col in required_columns if col not in columns]
    if missing:
        raise ValueError(
            "{} missing required columns: {}".format(label, ", ".join(missing))
        )


def _assert_unique(con, table, keys, label):
    key_sql = ", ".join(keys)
    duplicate = con.execute(
        "SELECT {keys} FROM {table} GROUP BY {keys} HAVING count(*) > 1 LIMIT 1".format(
            keys=key_sql, table=table
        )
    ).fetchone()
    if duplicate is not None:
        raise ValueError(
            "{} has duplicate rows for ({}): {}".format(label, key_sql, duplicate)
        )


def _require_columns(actual_columns, required_columns, label):
    missing = [col for col in required_columns if col not in actual_columns]
    if missing:
        raise ValueError(
            "{} missing required columns: {}".format(label, ", ".join(missing))
        )


def _source_exists(path):
    return bool(path) and os.path.isfile(path)


def _parse_int(value, label):
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError("Invalid integer for {}: {}".format(label, value))


def _quote_duckdb_path(path):
    return os.path.abspath(path).replace("'", "''")


def _ensure_parent_dir(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        file.mkdir_p(parent)
