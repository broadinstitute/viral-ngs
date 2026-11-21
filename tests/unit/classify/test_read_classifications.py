import argparse
import csv
import gzip
from unittest.mock import patch

import pytest

duckdb = pytest.importorskip("duckdb")

from viral_ngs import metagenomics
from viral_ngs.classify import genomad, read_classifications
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

    genomad_reads = tmp_path / "genomad.reads.tsv"
    _write_genomad_reads(genomad_reads, [_genomad_row("readA", "NODE_1")])

    out_parquet = tmp_path / "S1.read_classifications.parquet"
    read_classifications.join_read_classifications(
        str(out_parquet),
        "S1",
        kallisto_summary=str(kallisto),
        kraken2_reads=str(kraken),
        vnp_reads=str(vnp),
        genomad_reads=str(genomad_reads),
        centrifuger_reads=str(centrifuger),
    )

    con = duckdb.connect(database=":memory:")
    rows = con.execute(
        "SELECT * FROM read_parquet(?) ORDER BY READ_ID",
        [str(out_parquet)],
    ).fetchdf()
    con.close()

    assert rows["READ_ID"].tolist() == ["readA", "readB", "readC", "readD", "readZ"]
    read_a = rows[rows["READ_ID"] == "readA"].iloc[0]
    assert read_a["KALLISTO_CONCORDANCE"] == "Concordant"
    assert read_a["K2_TAX_NAME"] == "Homo sapiens"
    assert read_a["CENTRIFUGER_SEQ_ID"] == "seq1"
    assert read_a["GENOMAD_TAXONOMY"] == "Viruses"
    read_d = rows[rows["READ_ID"] == "readD"].iloc[0]
    assert read_d["KALLISTO_CONCORDANCE"] == "Discordant"
    assert "readHuman" not in set(rows["READ_ID"])


def test_join_all_inputs_missing_is_fatal(tmp_path):
    with pytest.raises(ValueError, match="At least one classification input must be provided"):
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
            "--genomad-reads",
            "g.tsv",
            "--no-filter-human-only",
            "--drop-unclassified",
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
            genomad_reads="g.tsv",
            centrifuger_reads=None,
            filter_human_only=False,
            drop_unclassified=True,
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


def _write_genomad_reads(path, rows):
    write_tsv(path, genomad.READ_LABEL_COLUMNS, rows)


def _genomad_row(read_id, seq_name="NODE_1", mapq=60):
    return {
        "read_id": read_id, "read_length": "100", "contig_id": seq_name.partition("|")[0],
        "mapping_quality": mapq, "pct_identity": "99.0", "pct_query_cov": "95.0", "mapped_well": "True",
        "seq_name": seq_name, "length": "1000", "topology": "linear",
        "coordinates": "1-1000", "n_genes": "10", "genetic_code": "11",
        "virus_score": "0.99", "fdr": "0.01", "n_hallmarks": "2",
        "marker_enrichment": "1.5", "taxonomy": "Viruses",
    }


def _taxon_row(read_id, taxid):
    name, kingdom, rank = {
        0: ("Unclassified", "Unclassified", "unclassified"),
        2: ("Bacteria", "Bacteria", "superkingdom"),
        9606: ("Homo sapiens", "Eukaryota", "species"),
    }[taxid]
    return {"SAMPLE_ID": "S1", "READ_ID": read_id, "TAXONOMY_ID": taxid,
            "TAX_NAME": name, "KINGDOM": kingdom, "TAX_RANK": rank}


def _centrifuger_row(read_id, taxid):
    return dict(_taxon_row(read_id, taxid), SEQ_ID="seq1", SCORE=42, SECOND_BEST_SCORE=5,
                HIT_LENGTH=77, QUERY_LENGTH=100, NUM_MATCHES=12)


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


def test_join_genomad_reads_collapse_mates(tmp_path):
    labels = tmp_path / "g.tsv"
    _write_genomad_reads(labels, [_genomad_row("readA/1", "NODE_1", 10),
                                 _genomad_row("readA/2", "NODE_2|provirus_201_400", 60)])
    out = tmp_path / "out.parquet"
    read_classifications.join_read_classifications(str(out), "S1", genomad_reads=str(labels))
    assert _parquet_rows(out, "READ_ID, GENOMAD_CONTIG_ID, GENOMAD_SEQ_NAME") == [
        ("readA", "NODE_2", "NODE_2|provirus_201_400")
    ]


@pytest.mark.parametrize(("source", "prefix", "header", "row"), [
    ("kallisto_summary", "KALLISTO_", read_classifications.KALLISTO_REQUIRED_COLUMNS,
     {"SAMPLE_ID": "S1", "READ_ID": "readA/1", "DB_ID": "u1", "TAXONOMY_LINEAGE": "L",
      "TAXONOMY_NAME": "N", "SEQUENCE_LENGTH": 100}),
    ("kraken2_reads", "K2_", read_classifications.KRAKEN2_ANNOTATED_COLUMNS, _taxon_row("readA/1", 2)),
    ("centrifuger_reads", "CENTRIFUGER_", read_classifications.CENTRIFUGER_ANNOTATED_COLUMNS,
     _centrifuger_row("readA/1", 2)),
    ("vnp_reads", "VNP_", read_classifications.VNP_REQUIRED_COLUMNS, _vnp_row("readA/1")),
    ("genomad_reads", "GENOMAD_", genomad.READ_LABEL_COLUMNS, _genomad_row("readA/1")),
], ids=["kallisto", "kraken2", "centrifuger", "virnucpro", "genomad"])
def test_join_each_source_alone_produces_rows(tmp_path, source, prefix, header, row):
    path, out = tmp_path / "source.tsv", tmp_path / "out.parquet"
    write_tsv(path, header, [row])
    read_classifications.join_read_classifications(str(out), "S1", **{source: str(path)})
    result, = _parquet_rows(out)
    assert result[:2] == ("S1", "readA")
    for name, value in zip(read_classifications.FINAL_COLUMNS[2:], result[2:]):
        assert (value is not None) == name.startswith(prefix), name


def _joined_ids(tmp_path, **kwargs):
    out = tmp_path / "out.parquet"
    read_classifications.join_read_classifications(str(out), "S1", **kwargs)
    return {read_id for read_id, in _parquet_rows(out, "READ_ID")}


def test_join_human_only_filter(tmp_path):
    k2, cfr, gen = (tmp_path / name for name in ("k2.tsv", "cfr.tsv", "gen.tsv"))
    write_tsv(k2, read_classifications.KRAKEN2_ANNOTATED_COLUMNS, [
        _taxon_row("k2-human", 9606), _taxon_row("human-zero", 9606),
        _taxon_row("human-other", 9606), _taxon_row("human-gen", 9606),
        dict(_taxon_row("no-name", 2), TAX_NAME=""),
    ])
    write_tsv(cfr, read_classifications.CENTRIFUGER_ANNOTATED_COLUMNS, [
        _centrifuger_row("cfr-human", 9606), _centrifuger_row("human-zero", 0),
        _centrifuger_row("human-other", 2),
    ])
    _write_genomad_reads(gen, [_genomad_row("human-gen")])
    sources = dict(kraken2_reads=str(k2), centrifuger_reads=str(cfr), genomad_reads=str(gen))
    retained = {"human-other", "human-gen", "no-name"}
    assert _joined_ids(tmp_path, **sources) == retained
    assert _joined_ids(tmp_path, **sources, filter_human_only=False) == retained | {"k2-human", "cfr-human", "human-zero"}


def test_join_drop_unclassified(tmp_path):
    k2, cfr, vnp = (tmp_path / name for name in ("k2.tsv", "cfr.tsv", "vnp.tsv"))
    write_tsv(k2, read_classifications.KRAKEN2_ANNOTATED_COLUMNS,
              [_taxon_row("k2-zero", 0), _taxon_row("with-vnp", 0)])
    write_tsv(cfr, read_classifications.CENTRIFUGER_ANNOTATED_COLUMNS, [_centrifuger_row("cfr-zero", 0)])
    write_tsv(vnp, read_classifications.VNP_REQUIRED_COLUMNS, [_vnp_row("with-vnp")])
    sources = dict(kraken2_reads=str(k2), centrifuger_reads=str(cfr), vnp_reads=str(vnp))
    assert _joined_ids(tmp_path, **sources) == {"k2-zero", "cfr-zero", "with-vnp"}
    assert _joined_ids(tmp_path, **sources, drop_unclassified=True) == {"with-vnp"}


def test_join_sample_id_mismatch_is_fatal(tmp_path):
    k2 = tmp_path / "k2.tsv"
    write_tsv(k2, read_classifications.KRAKEN2_ANNOTATED_COLUMNS, [_taxon_row("r1", 2)])
    with pytest.raises(ValueError, match="none for SAMPLE_ID 'S2'"):
        read_classifications.join_read_classifications(str(tmp_path / "out.parquet"), "S2", kraken2_reads=str(k2))


def test_join_mixed_samples_only_uses_requested_sample(tmp_path, caplog):
    cfr = tmp_path / "cfr.tsv"
    write_tsv(cfr, read_classifications.CENTRIFUGER_ANNOTATED_COLUMNS, [
        _centrifuger_row("r1", 2), dict(_centrifuger_row("r2", 2), SAMPLE_ID="S2"),
    ])
    assert _joined_ids(tmp_path, centrifuger_reads=str(cfr)) == {"r1"}
    assert "ignoring 1 rows" in caplog.text


def test_join_missing_path_is_fatal(tmp_path):
    k2 = tmp_path / "k2.tsv"
    write_tsv(k2, read_classifications.KRAKEN2_ANNOTATED_COLUMNS, [_taxon_row("readA", 2)])
    with pytest.raises(FileNotFoundError, match="missing.tsv"):
        _joined_ids(tmp_path, kraken2_reads=str(k2), genomad_reads=str(tmp_path / "missing.tsv"))


def test_join_zero_byte_input_is_empty_source(tmp_path):
    k2, gen, out = (tmp_path / name for name in ("k2.tsv", "empty.tsv", "out.parquet"))
    write_tsv(k2, read_classifications.KRAKEN2_ANNOTATED_COLUMNS, [_taxon_row("readA", 2)])
    gen.touch()
    read_classifications.join_read_classifications(str(out), "S1", kraken2_reads=str(k2), genomad_reads=str(gen))
    gen_cols = [c for c in read_classifications.FINAL_COLUMNS if c.startswith("GENOMAD_")]
    assert _parquet_rows(out, ", ".join(["READ_ID", "K2_TAXONOMY_ID", *gen_cols])) == [
        ("readA", 2, *([None] * len(gen_cols)))
    ]


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
        types = [str(desc[1]) for desc in result.description]
        n_rows = len(result.fetchall())
    finally:
        con.close()
    assert n_rows == 0
    assert columns == read_classifications.FINAL_COLUMNS
    assert list(zip(columns, types)) == read_classifications.FINAL_SCHEMA


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


def _write_taxdump(tmp_path, compressed=False):
    suffix = ".gz" if compressed else ""
    nodes, names = (tmp_path / (name + suffix) for name in ("nodes.dmp", "names.dmp"))
    opener = gzip.open if compressed else open
    with opener(nodes, "wt") as nf, opener(names, "wt") as sf:
        for taxid, parent, rank, name in [
            (1, 1, "no rank", "root"),
            (10239, 1, "acellular root", "Viruses"),
            (2559587, 10239, "realm", "Riboviria"),
            (2697049, 2559587, "species", "SARS-CoV-2"),
            (999999, 2697049, "no rank", "strain X"),
            (131567, 1, "cellular root", "cellular organisms"),
            (2759, 131567, "domain", "Eukaryota"),
            (9605, 2759, "genus", "Homo"),
            (9606, 9605, "species", "Homo sapiens"),
            (2157, 131567, "superkingdom", "Archaea"),
        ]:
            nf.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\n")
            sf.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    return str(nodes), str(names)


def test_build_taxonomy_duckdb_kingdoms_and_round_trip(tmp_path):
    nodes, names = _write_taxdump(tmp_path)
    db = str(tmp_path / "taxonomy.duckdb")
    read_classifications.build_taxonomy_duckdb(nodes, names, db)
    with duckdb.connect(db, read_only=True) as con:
        assert dict(con.execute("SELECT taxid, kingdom FROM taxonomy").fetchall()) == {
            1: "Other", 10239: "Viruses", 2559587: "Viruses", 2697049: "Viruses",
            999999: "Viruses", 131567: "Other", 2759: "Eukaryota", 9605: "Eukaryota",
            9606: "Eukaryota", 2157: "Archaea",
        }
        assert con.execute("SELECT parent_id FROM taxonomy WHERE taxid=1").fetchone() == (1,)
    native, out = tmp_path / "k2.reads", tmp_path / "k2.tsv"
    native.write_text("C\tr1\t999999\t100\t0:1\n")
    read_classifications.kraken2_annotate_reads(str(native), db, str(out), "S1", resolve_strains=True)
    row, = read_tsv(out)
    assert (row["TAXONOMY_ID"], row["TAX_RANK"], row["KINGDOM"]) == ("2697049", "species", "Viruses")


def test_taxonomy_duckdb_build_parser_runs(tmp_path):
    nodes, names = _write_taxdump(tmp_path, compressed=True)
    out = tmp_path / "reference" / "tax.duckdb"
    parser = metagenomics.parser_taxonomy_duckdb_build(argparse.ArgumentParser())
    args = parser.parse_args([nodes, names, str(out)])
    args.func_main(args)
    with duckdb.connect(str(out), read_only=True) as con:
        assert con.execute("SELECT name, kingdom FROM taxonomy WHERE taxid=9606").fetchone() == ("Homo sapiens", "Eukaryota")
