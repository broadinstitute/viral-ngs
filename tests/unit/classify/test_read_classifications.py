import argparse
import csv
from unittest.mock import patch

import pytest

duckdb = pytest.importorskip("duckdb")

from viral_ngs import metagenomics
from viral_ngs.classify import read_classifications
from viral_ngs.core import file as util_file


@pytest.fixture
def taxonomy_duckdb(tmp_path):
    path = tmp_path / "taxonomy.duckdb"
    con = duckdb.connect(str(path))
    con.execute(
        """
        CREATE TABLE taxonomy (
            taxid BIGINT,
            name VARCHAR,
            kingdom VARCHAR,
            rank VARCHAR,
            parent_id BIGINT
        )
        """
    )
    con.executemany(
        "INSERT INTO taxonomy VALUES (?, ?, ?, ?, ?)",
        [
            (0, "Unclassified", "Unclassified", "unclassified", None),
            (2, "Bacteria", "Bacteria", "superkingdom", 1),
            (9606, "Homo sapiens", "Eukaryota", "species", 9605),
            (2697049, "Severe acute respiratory syndrome coronavirus 2", "Viruses", "species", 1),
            (999999, "SARS-CoV-2 strain X", "Viruses", "no rank", 2697049),
        ],
    )
    con.close()
    return str(path)


def read_tsv(path):
    with util_file.open_or_gzopen(str(path), "rt") as inf:
        return list(csv.DictReader(inf, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "wt") as outf:
        writer = csv.DictWriter(outf, header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_kraken2_annotate_reads_writes_normalized_tsv(tmp_path, taxonomy_duckdb):
    native = tmp_path / "kraken2.reads"
    native.write_text(
        "C\treadA/1\t9606\t100\tA:1\n"
        "U\treadB\t0\t100\t0:1\n"
        "C\treadC\t999999\t100\tA:1\n"
    )
    output = tmp_path / "kraken2.read_taxonomy.tsv"

    read_classifications.kraken2_annotate_reads(
        str(native),
        taxonomy_duckdb,
        str(output),
        "S1",
        resolve_strains=True,
    )

    rows = read_tsv(output)
    assert rows[0]["SAMPLE_ID"] == "S1"
    assert rows[0]["READ_ID"] == "readA/1"
    assert rows[0]["TAX_NAME"] == "Homo sapiens"
    assert rows[1]["TAX_NAME"] == "Unclassified"
    assert rows[2]["TAX_NAME"] == "Severe acute respiratory syndrome coronavirus 2"


def test_centrifuger_annotate_reads_writes_normalized_tsv(tmp_path, taxonomy_duckdb):
    native = tmp_path / "centrifuger.tsv"
    write_tsv(
        native,
        read_classifications.CENTRIFUGER_NATIVE_COLUMNS,
        [
            {
                "readID": "readA/1",
                "seqID": "seq1",
                "taxID": "2",
                "score": "42",
                "2ndBestScore": "5",
                "hitLength": "77",
                "queryLength": "100",
                "numMatches": "12",
            }
        ],
    )
    output = tmp_path / "centrifuger.read_taxonomy.tsv.zst"

    read_classifications.centrifuger_annotate_reads(
        str(native),
        taxonomy_duckdb,
        str(output),
        "S1",
    )

    rows = read_tsv(output)
    assert rows == [
        {
            "SAMPLE_ID": "S1",
            "READ_ID": "readA/1",
            "SEQ_ID": "seq1",
            "TAXONOMY_ID": "2",
            "SCORE": "42",
            "SECOND_BEST_SCORE": "5",
            "HIT_LENGTH": "77",
            "QUERY_LENGTH": "100",
            "NUM_MATCHES": "12",
            "TAX_NAME": "Bacteria",
            "KINGDOM": "Bacteria",
            "TAX_RANK": "superkingdom",
        }
    ]


def test_annotate_empty_native_outputs_header_only(tmp_path, taxonomy_duckdb):
    native = tmp_path / "empty.kraken"
    native.write_text("")
    output = tmp_path / "kraken2.read_taxonomy.tsv"

    read_classifications.kraken2_annotate_reads(
        str(native),
        taxonomy_duckdb,
        str(output),
        "S1",
    )

    assert output.read_text().strip() == "\t".join(
        read_classifications.KRAKEN2_ANNOTATED_COLUMNS
    )


def test_join_read_classifications_component_path(tmp_path, taxonomy_duckdb):
    kallisto = tmp_path / "kallisto.summary.tsv"
    write_tsv(
        kallisto,
        read_classifications.KALLISTO_REQUIRED_COLUMNS,
        [
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readA/1",
                "DB_ID": "u1",
                "TAXONOMY_LINEAGE": "Viruses; Coronaviridae",
                "TAXONOMY_NAME": "Coronaviridae",
                "SEQUENCE_LENGTH": "100",
            },
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readA/2",
                "DB_ID": "u1",
                "TAXONOMY_LINEAGE": "Viruses; Coronaviridae",
                "TAXONOMY_NAME": "Coronaviridae",
                "SEQUENCE_LENGTH": "90",
            },
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readD/1",
                "DB_ID": "u2",
                "TAXONOMY_LINEAGE": "Viruses; A",
                "TAXONOMY_NAME": "A",
                "SEQUENCE_LENGTH": "80",
            },
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readD/2",
                "DB_ID": "u3",
                "TAXONOMY_LINEAGE": "Viruses; B",
                "TAXONOMY_NAME": "B",
                "SEQUENCE_LENGTH": "70",
            },
        ],
    )

    vnp = tmp_path / "vnp.tsv"
    write_tsv(
        vnp,
        read_classifications.VNP_REQUIRED_COLUMNS,
        [
            {
                "read_id": "readA/1",
                "read_length": "100",
                "contig_id": "NODE_1",
                "contig_length": "1000",
                "strand": "+",
                "mapping_quality": "60",
                "pct_identity": "99.0",
                "pct_query_cov": "95.0",
                "mapped_well": "True",
                "call": "Viral",
                "tier": "high_confidence",
                "weighted_delta": "0.9",
                "n_chunks": "5",
                "n_confident_viral": "5",
                "n_confident_nonviral": "0",
                "n_ambiguous": "0",
                "viral_proportion": "1.0",
                "nonviral_proportion": "0.0",
            },
            {
                "read_id": "readB/1",
                "read_length": "100",
                "contig_id": "NODE_2",
                "contig_length": "900",
                "strand": "+",
                "mapping_quality": "50",
                "pct_identity": "98.0",
                "pct_query_cov": "90.0",
                "mapped_well": "True",
                "call": "Non-viral",
                "tier": "high_confidence",
                "weighted_delta": "-0.9",
                "n_chunks": "5",
                "n_confident_viral": "0",
                "n_confident_nonviral": "5",
                "n_ambiguous": "0",
                "viral_proportion": "0.0",
                "nonviral_proportion": "1.0",
            },
        ],
    )

    kraken_native = tmp_path / "kraken2.reads"
    kraken_native.write_text(
        "C\treadA/1\t9606\t100\tA:1\n"
        "C\treadC/1\t2\t100\tA:1\n"
        "C\treadHuman/1\t9606\t100\tA:1\n"
    )
    kraken = tmp_path / "S1.kraken2.read_taxonomy.tsv"
    read_classifications.kraken2_annotate_reads(
        str(kraken_native),
        taxonomy_duckdb,
        str(kraken),
        "S1",
    )

    centrifuger_native = tmp_path / "centrifuger.tsv"
    write_tsv(
        centrifuger_native,
        read_classifications.CENTRIFUGER_NATIVE_COLUMNS,
        [
            {
                "readID": "readA/1",
                "seqID": "seq1",
                "taxID": "2",
                "score": "42",
                "2ndBestScore": "5",
                "hitLength": "77",
                "queryLength": "100",
                "numMatches": "12",
            },
            {
                "readID": "readZ/1",
                "seqID": "seqZ",
                "taxID": "2",
                "score": "40",
                "2ndBestScore": "5",
                "hitLength": "77",
                "queryLength": "100",
                "numMatches": "12",
            },
        ],
    )
    centrifuger = tmp_path / "S1.centrifuger.read_taxonomy.tsv"
    read_classifications.centrifuger_annotate_reads(
        str(centrifuger_native),
        taxonomy_duckdb,
        str(centrifuger),
        "S1",
    )

    genomad = tmp_path / "virus_summary.tsv"
    write_tsv(
        genomad,
        [
            "seq_name",
            "length",
            "topology",
            "coordinates",
            "n_genes",
            "genetic_code",
            "virus_score",
            "fdr",
            "n_hallmarks",
            "marker_enrichment",
            "taxonomy",
        ],
        [
            {
                "seq_name": "NODE_1",
                "length": "1000",
                "topology": "linear",
                "coordinates": "1-1000",
                "n_genes": "10",
                "genetic_code": "11",
                "virus_score": "0.99",
                "fdr": "0.01",
                "n_hallmarks": "2",
                "marker_enrichment": "1.5",
                "taxonomy": "Viruses",
            }
        ],
    )

    out_parquet = tmp_path / "S1.read_classifications.parquet"
    read_classifications.join_read_classifications(
        str(out_parquet),
        "S1",
        kallisto_summary=str(kallisto),
        kraken2_reads=str(kraken),
        vnp_reads=str(vnp),
        genomad_virus_summary=str(genomad),
        centrifuger_reads=str(centrifuger),
    )

    con = duckdb.connect(database=":memory:")
    rows = con.execute(
        "SELECT * FROM read_parquet(?) ORDER BY READ_ID",
        [str(out_parquet)],
    ).fetchdf()
    con.close()

    assert rows["READ_ID"].tolist() == ["readA", "readB", "readC", "readD"]
    read_a = rows[rows["READ_ID"] == "readA"].iloc[0]
    assert read_a["KALLISTO_CONCORDANCE"] == "Concordant"
    assert read_a["K2_TAX_NAME"] == "Homo sapiens"
    assert read_a["CENTRIFUGER_SEQ_ID"] == "seq1"
    assert read_a["GENOMAD_TAXONOMY"] == "Viruses"
    read_d = rows[rows["READ_ID"] == "readD"].iloc[0]
    assert read_d["KALLISTO_CONCORDANCE"] == "Discordant"
    assert "readHuman" not in set(rows["READ_ID"])


def test_join_all_inputs_missing_is_fatal(tmp_path):
    with pytest.raises(ValueError, match="At least one classification input"):
        read_classifications.join_read_classifications(
            str(tmp_path / "out.parquet"),
            "S1",
        )


def test_join_rejects_duplicate_centrifuger_keys(tmp_path):
    existing = tmp_path / "kraken.tsv"
    write_tsv(
        existing,
        read_classifications.KRAKEN2_ANNOTATED_COLUMNS,
        [
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readA",
                "TAXONOMY_ID": "2",
                "TAX_NAME": "Bacteria",
                "KINGDOM": "Bacteria",
                "TAX_RANK": "superkingdom",
            }
        ],
    )
    centrifuger = tmp_path / "centrifuger.tsv"
    write_tsv(
        centrifuger,
        read_classifications.CENTRIFUGER_ANNOTATED_COLUMNS,
        [
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readA/1",
                "SEQ_ID": "seq1",
                "TAXONOMY_ID": "2",
                "SCORE": "1",
                "SECOND_BEST_SCORE": "0",
                "HIT_LENGTH": "10",
                "QUERY_LENGTH": "10",
                "NUM_MATCHES": "1",
                "TAX_NAME": "Bacteria",
                "KINGDOM": "Bacteria",
                "TAX_RANK": "superkingdom",
            },
            {
                "SAMPLE_ID": "S1",
                "READ_ID": "readA/2",
                "SEQ_ID": "seq2",
                "TAXONOMY_ID": "2",
                "SCORE": "1",
                "SECOND_BEST_SCORE": "0",
                "HIT_LENGTH": "10",
                "QUERY_LENGTH": "10",
                "NUM_MATCHES": "1",
                "TAX_NAME": "Bacteria",
                "KINGDOM": "Bacteria",
                "TAX_RANK": "superkingdom",
            },
        ],
    )

    with pytest.raises(ValueError, match="Centrifuger reads has duplicate"):
        read_classifications.join_read_classifications(
            str(tmp_path / "out.parquet"),
            "S1",
            kraken2_reads=str(existing),
            centrifuger_reads=str(centrifuger),
        )


def test_metagenomics_parsers_dispatch_new_commands():
    with patch(
        "viral_ngs.metagenomics.read_classifications.kraken2_annotate_reads",
        autospec=True,
    ) as mock_kraken:
        args = metagenomics.parser_kraken2_annotate_reads(
            argparse.ArgumentParser()
        ).parse_args([
            "k2.tsv",
            "tax.duckdb",
            "out.tsv",
            "--sample-id",
            "S1",
            "--resolve-strains",
        ])
        args.func_main(args)
        mock_kraken.assert_called_once_with(
            "k2.tsv",
            "tax.duckdb",
            "out.tsv",
            "S1",
            resolve_strains=True,
        )

    with patch(
        "viral_ngs.metagenomics.read_classifications.centrifuger_annotate_reads",
        autospec=True,
    ) as mock_centrifuger:
        args = metagenomics.parser_centrifuger_annotate_reads(
            argparse.ArgumentParser()
        ).parse_args([
            "cfr.tsv",
            "tax.duckdb",
            "out.tsv",
            "--sample-id",
            "S1",
        ])
        args.func_main(args)
        mock_centrifuger.assert_called_once_with(
            "cfr.tsv",
            "tax.duckdb",
            "out.tsv",
            "S1",
            resolve_strains=False,
        )

    with patch(
        "viral_ngs.metagenomics.read_classifications.join_read_classifications",
        autospec=True,
    ) as mock_join:
        args = metagenomics.parser_join_read_classifications(
            argparse.ArgumentParser()
        ).parse_args([
            "--sample-id",
            "S1",
            "--out-parquet",
            "out.parquet",
            "--kraken2-reads",
            "k2.tsv",
            "--no-filter-human-only-k2",
            "--duckdb-memory-limit",
            "4GB",
            "--work-dir",
            "/tmp/work",
        ])
        args.func_main(args)
        mock_join.assert_called_once_with(
            "out.parquet",
            "S1",
            kallisto_summary=None,
            kraken2_reads="k2.tsv",
            vnp_reads=None,
            genomad_virus_summary=None,
            centrifuger_reads=None,
            filter_human_only_k2=False,
            duckdb_memory_limit="4GB",
            work_dir="/tmp/work",
        )

def _vnp_row(read_id, contig_id="NODE_1", call="Viral"):
    return {
        "read_id": read_id,
        "read_length": "100",
        "contig_id": contig_id,
        "contig_length": "1000",
        "strand": "+",
        "mapping_quality": "60",
        "pct_identity": "99.0",
        "pct_query_cov": "95.0",
        "mapped_well": "True",
        "call": call,
        "tier": "high_confidence",
        "weighted_delta": "0.9",
        "n_chunks": "5",
        "n_confident_viral": "5",
        "n_confident_nonviral": "0",
        "n_ambiguous": "0",
        "viral_proportion": "1.0",
        "nonviral_proportion": "0.0",
    }


_GENOMAD_HEADER = [
    "seq_name", "length", "topology", "coordinates", "n_genes", "genetic_code",
    "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy",
]


def _genomad_row(seq_name):
    return {
        "seq_name": seq_name, "length": "1000", "topology": "linear",
        "coordinates": "1-1000", "n_genes": "10", "genetic_code": "11",
        "virus_score": "0.99", "fdr": "0.01", "n_hallmarks": "2",
        "marker_enrichment": "1.5", "taxonomy": "Viruses",
    }


def _parquet_rows(path, columns="*", order="READ_ID"):
    con = duckdb.connect(database=":memory:")
    try:
        return con.execute(
            "SELECT {} FROM read_parquet(?) ORDER BY {}".format(columns, order),
            [str(path)],
        ).fetchall()
    finally:
        con.close()


def test_kallisto_single_mate_multiple_db_is_singleton(tmp_path):
    kallisto = tmp_path / "kallisto.tsv"
    write_tsv(
        kallisto,
        read_classifications.KALLISTO_REQUIRED_COLUMNS,
        [
            {"SAMPLE_ID": "S1", "READ_ID": "readX/1", "DB_ID": "u1",
             "TAXONOMY_LINEAGE": "L", "TAXONOMY_NAME": "N", "SEQUENCE_LENGTH": "100"},
            {"SAMPLE_ID": "S1", "READ_ID": "readX/1", "DB_ID": "u2",
             "TAXONOMY_LINEAGE": "L", "TAXONOMY_NAME": "N", "SEQUENCE_LENGTH": "90"},
        ],
    )
    out = tmp_path / "out.parquet"
    read_classifications.join_read_classifications(
        str(out), "S1", kallisto_summary=str(kallisto)
    )
    # One mate but two DB hits must collapse to Singleton, not Discordant.
    assert _parquet_rows(out, "READ_ID, KALLISTO_CONCORDANCE") == [("readX", "Singleton")]


def test_join_rejects_duplicate_kraken2_keys(tmp_path):
    kraken = tmp_path / "k2.tsv"
    write_tsv(
        kraken,
        read_classifications.KRAKEN2_ANNOTATED_COLUMNS,
        [
            {"SAMPLE_ID": "S1", "READ_ID": "readA/1", "TAXONOMY_ID": "2",
             "TAX_NAME": "Bacteria", "KINGDOM": "Bacteria", "TAX_RANK": "superkingdom"},
            {"SAMPLE_ID": "S1", "READ_ID": "readA/2", "TAXONOMY_ID": "2",
             "TAX_NAME": "Bacteria", "KINGDOM": "Bacteria", "TAX_RANK": "superkingdom"},
        ],
    )
    with pytest.raises(ValueError, match="Kraken2 reads has duplicate"):
        read_classifications.join_read_classifications(
            str(tmp_path / "out.parquet"), "S1", kraken2_reads=str(kraken)
        )


def test_join_rejects_duplicate_genomad_contigs(tmp_path):
    vnp = tmp_path / "vnp.tsv"
    write_tsv(vnp, read_classifications.VNP_REQUIRED_COLUMNS, [_vnp_row("readA/1")])
    genomad = tmp_path / "g.tsv"
    write_tsv(genomad, _GENOMAD_HEADER, [_genomad_row("NODE_1"), _genomad_row("NODE_1")])
    with pytest.raises(ValueError, match="geNomad virus summary has duplicate"):
        read_classifications.join_read_classifications(
            str(tmp_path / "out.parquet"), "S1",
            vnp_reads=str(vnp), genomad_virus_summary=str(genomad),
        )


def test_join_genomad_enriches_without_multiplying_rows(tmp_path):
    vnp = tmp_path / "vnp.tsv"
    write_tsv(
        vnp,
        read_classifications.VNP_REQUIRED_COLUMNS,
        [_vnp_row("readA/1", "NODE_1"), _vnp_row("readB/1", "NODE_1")],
    )
    genomad = tmp_path / "g.tsv"
    write_tsv(genomad, _GENOMAD_HEADER, [_genomad_row("NODE_1")])
    out = tmp_path / "out.parquet"
    read_classifications.join_read_classifications(
        str(out), "S1", vnp_reads=str(vnp), genomad_virus_summary=str(genomad)
    )
    rows = _parquet_rows(out, "READ_ID, GENOMAD_TAXONOMY")
    # Two reads on one contig stay two rows, both annotated.
    assert rows == [("readA", "Viruses"), ("readB", "Viruses")]


def test_join_header_only_input_emits_empty_typed_parquet(tmp_path):
    kallisto = tmp_path / "kallisto.tsv"
    write_tsv(kallisto, read_classifications.KALLISTO_REQUIRED_COLUMNS, [])
    out = tmp_path / "out.parquet"
    read_classifications.join_read_classifications(
        str(out), "S1", kallisto_summary=str(kallisto)
    )
    con = duckdb.connect(database=":memory:")
    try:
        result = con.execute("SELECT * FROM read_parquet(?)", [str(out)])
        columns = [desc[0] for desc in result.description]
        n_rows = len(result.fetchall())
    finally:
        con.close()
    assert n_rows == 0
    assert columns == read_classifications.FINAL_COLUMNS


def test_annotate_missing_taxonomy_db_is_fatal(tmp_path):
    native = tmp_path / "k2.reads"
    native.write_text("C\treadA\t2\t100\tA:1\n")
    with pytest.raises(FileNotFoundError):
        read_classifications.kraken2_annotate_reads(
            str(native), str(tmp_path / "missing.duckdb"),
            str(tmp_path / "o.tsv"), "S1",
        )


def test_centrifuger_annotate_missing_taxonomy_db_is_fatal(tmp_path):
    native = tmp_path / "cfr.tsv"
    write_tsv(native, read_classifications.CENTRIFUGER_NATIVE_COLUMNS, [])
    with pytest.raises(FileNotFoundError):
        read_classifications.centrifuger_annotate_reads(
            str(native), str(tmp_path / "missing.duckdb"),
            str(tmp_path / "o.tsv"), "S1",
        )


def test_annotate_missing_native_output_is_fatal(tmp_path, taxonomy_duckdb):
    with pytest.raises(FileNotFoundError):
        read_classifications.kraken2_annotate_reads(
            str(tmp_path / "nope.reads"), taxonomy_duckdb,
            str(tmp_path / "o.tsv"), "S1",
        )


def test_centrifuger_annotate_missing_native_output_is_fatal(tmp_path, taxonomy_duckdb):
    with pytest.raises(FileNotFoundError):
        read_classifications.centrifuger_annotate_reads(
            str(tmp_path / "nope.tsv"), taxonomy_duckdb,
            str(tmp_path / "o.tsv"), "S1",
        )


def test_kraken2_annotate_malformed_layout_is_fatal(tmp_path, taxonomy_duckdb):
    native = tmp_path / "k2.reads"
    native.write_text("C\treadA\t9606\n")  # 3 columns, not 5
    with pytest.raises(ValueError, match="Malformed Kraken2 row"):
        read_classifications.kraken2_annotate_reads(
            str(native), taxonomy_duckdb, str(tmp_path / "o.tsv"), "S1"
        )


def test_centrifuger_annotate_malformed_layout_is_fatal(tmp_path, taxonomy_duckdb):
    native = tmp_path / "cfr.tsv"
    native.write_text(
        "\t".join(read_classifications.CENTRIFUGER_NATIVE_COLUMNS) + "\n"
        + "readA\tseq1\t2\n"  # too few columns
    )
    with pytest.raises(ValueError, match="Malformed Centrifuger row"):
        read_classifications.centrifuger_annotate_reads(
            str(native), taxonomy_duckdb, str(tmp_path / "o.tsv"), "S1"
        )


def test_annotate_unknown_taxid_falls_back(tmp_path, taxonomy_duckdb):
    native = tmp_path / "k2.reads"
    native.write_text("C\treadA\t424242\t100\tA:1\n")
    out = tmp_path / "o.tsv"
    read_classifications.kraken2_annotate_reads(
        str(native), taxonomy_duckdb, str(out), "S1"
    )
    rows = read_tsv(out)
    assert rows[0]["TAX_NAME"] == "Unknown"
    assert rows[0]["KINGDOM"] == "Unknown"
    assert rows[0]["TAX_RANK"] == "unknown"


def test_taxonomy_lookup_is_cached(taxonomy_duckdb):
    db = read_classifications.DuckDBTaxonomyDatabase(taxonomy_duckdb)
    try:
        first = db.lookup(2)
        second = db.lookup(2)
        assert first is second  # memoized result reused
        assert first["name"] == "Bacteria"
        assert db.lookup(0)["name"] == "Unclassified"
    finally:
        db.close()
