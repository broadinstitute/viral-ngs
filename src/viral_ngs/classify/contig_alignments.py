"""Normalize primary BAM alignments for read-by-contig classification."""

import os

import pysam


def validate_percent_threshold(name, value):
    if 0 < value < 1.0:
        raise ValueError(
            "{} must be specified as a percent, not a fraction. "
            "Use {} instead of {}.".format(name, value * 100.0, value)
        )
    if value < 0 or value > 100:
        raise ValueError("{} must be between 0 and 100.".format(name))


_BLOCK_CIGAR_OPS = {
    pysam.CMATCH,
    pysam.CINS,
    pysam.CDEL,
    pysam.CEQUAL,
    pysam.CDIFF,
}


def normalize_primary_alignments(aligned_bam, tmp_dir):
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
            row = _record_to_row(bam, record, n_written)
            outf.write("\t".join(row) + "\n")

    return normalized_alignments, n_written, n_secondary


def _record_to_row(bam, record, source_order):
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

    alignment_block_length = _block_length(record)
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


def _block_length(record):
    return sum(
        length
        for op, length in (record.cigartuples or [])
        if op in _BLOCK_CIGAR_OPS
    )


def load_primary_alignments(con, normalized_alignments, min_mapq, min_identity, min_query_cov):
    """Create alignment_filtered with quality flags, retaining every primary alignment."""
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
    return con.execute("SELECT count(*) FROM alignment_filtered").fetchone()[0]
