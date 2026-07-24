import json
import os
import subprocess
import sys

import pysam
import pytest

import viral_ngs.core.file as util_file
from viral_ngs import metagenomics
from viral_ngs.classify import lyra


_BAM_HEADER = {
    "HD": {"VN": "1.6", "SO": "unsorted"},
    "RG": [{"ID": "fixture", "SM": "sample matrix"}],
    "PG": [{"ID": "fixture", "PN": "lyra-integration"}],
}
FORBIDDEN_TOP_LEVEL_PREFIXES = (
    "lyra",
    "torch",
    "pytorch_lightning",
    "cuda",
    "cupy",
    "nvidia",
    "tensorflow",
    "transformers",
)
CHILD_CODE = r"""
import json
import os
import sys

from viral_ngs import metagenomics

score_path, bam_path, normalized, summary, viral_bam, modules_path = sys.argv[1:]
args = metagenomics.full_parser().parse_args([
    "lyra_postprocess",
    score_path,
    bam_path,
    "runtime sample",
    normalized,
    summary,
    viral_bam,
])
result = args.func_main(args)
if result is not None:
    raise RuntimeError("lyra_postprocess dispatcher returned a value")
if not all(os.path.isfile(path) for path in (normalized, summary, viral_bam)):
    raise RuntimeError("lyra_postprocess did not create every artifact")
forbidden_prefixes = (
    "lyra",
    "torch",
    "pytorch_lightning",
    "cuda",
    "cupy",
    "nvidia",
    "tensorflow",
    "transformers",
)
loaded = sorted(
    name
    for name in sys.modules
    if any(
        name == prefix or name.startswith(prefix + ".")
        for prefix in forbidden_prefixes
    )
)
with open(modules_path, "w", encoding="utf-8") as output:
    json.dump(loaded, output)
"""


def _segment(query_name, flag=0x4):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = flag
    record.query_sequence = "A" * 50
    record.query_qualities = pysam.qualitystring_to_array("I" * 50)
    record.reference_id = -1
    record.reference_start = -1
    record.next_reference_id = -1
    record.next_reference_start = -1
    return record


def _write_bam(path, records, header=None):
    with pysam.AlignmentFile(
        str(path),
        "wb",
        header=header or _BAM_HEADER,
    ) as bam:
        for record in records:
            bam.write(record)
    return path


def _write_score_table(path, rows):
    contents = "read_id\tscore\tcall\n" + "".join(
        "{}\t{}\t{}\n".format(*row) for row in rows
    )
    with util_file.open_or_gzopen(str(path), "wb") as output:
        output.write(contents.encode("utf-8"))
    return path


def _logical_tsv_bytes(path):
    with util_file.open_or_gzopen(str(path), "rb") as stream:
        return stream.read()


def _read_summary(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    header = tuple(lines[0].split("\t"))
    values = tuple(lines[1].split("\t"))
    assert header == lyra.SUMMARY_HEADER
    assert len(values) == len(header)
    return dict(zip(header, values))


def _record_snapshot(record):
    return (
        record.query_name,
        record.flag,
        record.reference_id,
        record.reference_start,
        record.mapping_quality,
        record.cigarstring,
        record.next_reference_id,
        record.next_reference_start,
        record.template_length,
        record.query_sequence,
        tuple(record.query_qualities),
        record.get_tags(with_value_type=True),
    )


def _read_bam_records(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
        header = bam.header.to_dict()
        records = [
            _record_snapshot(record)
            for record in bam.fetch(until_eof=True)
        ]
    return header, records


def _dispatch_lyra_postprocess(argv):
    args = metagenomics.full_parser().parse_args(["lyra_postprocess"] + argv)
    return args.func_main(args)


def _assert_summary_equations(row, normalized_rows, output_bam_records):
    counts = {
        name: int(row[name])
        for name in lyra.SUMMARY_HEADER[2:]
    }
    assert counts["LYRA_ELIGIBLE_BAM_RECORDS"] == counts["LYRA_SCORE_RECORDS"]
    assert counts["LYRA_SCORE_RECORDS"] == (
        counts["LYRA_SINGLE_END_FRAGMENTS"]
        + counts["LYRA_INCOMPLETE_PAIR_FRAGMENTS"]
        + 2 * counts["LYRA_COMPLETE_PAIR_FRAGMENTS"]
    )
    assert counts["LYRA_FRAGMENTS"] == (
        counts["LYRA_SINGLE_END_FRAGMENTS"]
        + counts["LYRA_COMPLETE_PAIR_FRAGMENTS"]
        + counts["LYRA_INCOMPLETE_PAIR_FRAGMENTS"]
    )
    assert counts["LYRA_FRAGMENTS"] == (
        counts["LYRA_VIRAL_FRAGMENT_CALLS"]
        + counts["LYRA_NONVIRAL_FRAGMENT_CALLS"]
    )
    assert counts["LYRA_FRAGMENTS"] == normalized_rows
    assert 0 <= counts["LYRA_OUTPUT_BAM_RECORDS"] <= counts[
        "LYRA_INPUT_BAM_RECORDS"
    ]
    assert (counts["LYRA_OUTPUT_BAM_RECORDS"] == 0) == (
        counts["LYRA_VIRAL_FRAGMENT_CALLS"] == 0
    )
    assert counts["LYRA_OUTPUT_BAM_RECORDS"] == output_bam_records


def _output_paths(tmp_path, prefix, normalized_suffix=".tsv"):
    return (
        tmp_path / (prefix + "-normalized" + normalized_suffix),
        tmp_path / (prefix + "-summary.tsv"),
        tmp_path / (prefix + "-viral.bam"),
    )


def _assert_paths_absent(paths):
    assert not any(os.path.lexists(path) for path in paths)


def test_full_parser_codec_matrix_has_identical_logical_content(tmp_path, capsys):
    input_suffixes = (".tsv", ".tsv.gz", ".tsv.zst")
    output_suffixes = (".tsv", ".tsv.gz", ".tsv.zst")
    score_rows = [
        ("pair-incomplete", "1", "1"),
        ("pair-viral", "0.9", "0"),
        ("single-threshold", "0.8", "0"),
        ("pair-nonviral", "0.79", "1"),
        ("pair-viral", "0.8", "0"),
        ("pair-nonviral", "0.9", "1"),
    ]
    source_bam = _write_bam(
        tmp_path / "matrix-source.bam",
        [
            _segment("pair-nonviral", flag=0x4 | 0x1 | 0x40),
            _segment("pair-viral", flag=0x4 | 0x1 | 0x80),
            _segment("single-threshold"),
            _segment("pair-incomplete", flag=0x4 | 0x1 | 0x40),
            _segment("pair-viral", flag=0x4 | 0x1 | 0x40),
            _segment("pair-nonviral", flag=0x4 | 0x1 | 0x80),
        ],
    )
    source_header, source_records = _read_bam_records(source_bam)
    expected_bam_records = [
        record
        for record in source_records
        if record[0] in {"single-threshold", "pair-viral"}
    ]
    expected_normalized = (
        "SAMPLE_ID\tREAD_ID\tLYRA_N_SCORES\tLYRA_PAIRING\t"
        "LYRA_MIN_SCORE\tLYRA_MAX_SCORE\tLYRA_THRESHOLD\tLYRA_CALL\n"
        "sample matrix\tpair-incomplete\t1\tPaired-incomplete\t1\t1\t0.8\tNon-viral\n"
        "sample matrix\tpair-nonviral\t2\tPaired-complete\t0.79\t0.9\t0.8\tNon-viral\n"
        "sample matrix\tpair-viral\t2\tPaired-complete\t0.8\t0.9\t0.8\tViral\n"
        "sample matrix\tsingle-threshold\t1\tSingle-end\t0.8\t0.8\t0.8\tViral\n"
    ).encode("utf-8")
    expected_summary = (
        "\t".join(lyra.SUMMARY_HEADER)
        + "\n"
        + "sample matrix\t0.8\t6\t6\t6\t4\t1\t2\t1\t2\t2\t3\n"
    ).encode("utf-8")

    score_paths = {
        suffix: _write_score_table(
            tmp_path / ("matrix-input" + suffix),
            score_rows,
        )
        for suffix in input_suffixes
    }
    normalized_contents = []
    for input_index, input_suffix in enumerate(input_suffixes):
        for output_index, output_suffix in enumerate(output_suffixes):
            cell = "{}-{}".format(input_index, output_index)
            normalized = tmp_path / ("matrix-normalized-" + cell + output_suffix)
            summary = tmp_path / ("matrix-summary-" + cell + ".tsv")
            viral_bam = tmp_path / ("matrix-viral-" + cell + ".bam")

            result = _dispatch_lyra_postprocess([
                str(score_paths[input_suffix]),
                str(source_bam),
                "sample matrix",
                str(normalized),
                str(summary),
                str(viral_bam),
            ])

            assert result is None
            assert capsys.readouterr().out == ""
            logical_content = _logical_tsv_bytes(normalized)
            normalized_contents.append(logical_content)
            assert logical_content == expected_normalized
            assert not logical_content.startswith(b"\xef\xbb\xbf")
            assert summary.read_bytes() == expected_summary
            summary_row = _read_summary(summary)
            assert summary_row == dict(zip(
                lyra.SUMMARY_HEADER,
                expected_summary.decode("utf-8").splitlines()[1].split("\t"),
            ))

            output_header, output_records = _read_bam_records(viral_bam)
            assert output_header == source_header
            assert output_records == expected_bam_records
            assert len(output_records) == 3
            assert not {
                "pair-nonviral",
                "pair-incomplete",
            }.intersection(record[0] for record in output_records)
            _assert_summary_equations(
                summary_row,
                normalized_rows=4,
                output_bam_records=len(output_records),
            )

    assert len(normalized_contents) == 9
    assert all(content == normalized_contents[0] for content in normalized_contents)


def test_full_parser_explicit_threshold_exact_score_is_viral(tmp_path, capsys):
    score_path = _write_score_table(
        tmp_path / "threshold-scores.tsv",
        [("threshold", "0.8", "0")],
    )
    source_bam = _write_bam(
        tmp_path / "threshold-source.bam",
        [_segment("threshold")],
    )
    outputs = _output_paths(tmp_path, "threshold")

    result = _dispatch_lyra_postprocess([
        str(score_path),
        str(source_bam),
        "threshold sample",
        *(str(path) for path in outputs),
        "--threshold",
        "0.8",
    ])

    assert result is None
    assert capsys.readouterr().out == ""
    assert _logical_tsv_bytes(outputs[0]) == (
        "\t".join(lyra.NORMALIZED_HEADER)
        + "\nthreshold sample\tthreshold\t1\tSingle-end\t0.8\t0.8\t0.8\tViral\n"
    ).encode("utf-8")
    summary = _read_summary(outputs[1])
    assert summary == dict(zip(
        lyra.SUMMARY_HEADER,
        ("threshold sample", "0.8", "1", "1", "1", "1", "1", "0", "0", "1", "0", "1"),
    ))
    source_header, source_records = _read_bam_records(source_bam)
    output_header, output_records = _read_bam_records(outputs[2])
    assert output_header == source_header
    assert output_records == source_records
    _assert_summary_equations(summary, normalized_rows=1, output_bam_records=1)


def test_full_parser_empty_input_writes_complete_artifacts(tmp_path, capsys):
    score_path = _write_score_table(tmp_path / "empty-scores.tsv", [])
    source_bam = _write_bam(tmp_path / "empty-source.bam", [])
    outputs = _output_paths(tmp_path, "empty")

    result = _dispatch_lyra_postprocess([
        str(score_path),
        str(source_bam),
        "empty sample",
        *(str(path) for path in outputs),
    ])

    assert result is None
    assert capsys.readouterr().out == ""
    assert _logical_tsv_bytes(outputs[0]) == (
        "\t".join(lyra.NORMALIZED_HEADER) + "\n"
    ).encode("utf-8")
    expected_summary = (
        "\t".join(lyra.SUMMARY_HEADER)
        + "\nempty sample\t0.8\t"
        + "\t".join(["0"] * 10)
        + "\n"
    ).encode("utf-8")
    assert outputs[1].read_bytes() == expected_summary
    summary = _read_summary(outputs[1])
    source_header, _ = _read_bam_records(source_bam)
    output_header, output_records = _read_bam_records(outputs[2])
    assert output_header == source_header
    assert output_records == []
    _assert_summary_equations(summary, normalized_rows=0, output_bam_records=0)


def test_full_parser_no_hit_retains_nonviral_row(tmp_path, capsys):
    score_path = _write_score_table(
        tmp_path / "no-hit-scores.tsv",
        [("nonviral", "0.1", "1")],
    )
    source_bam = _write_bam(
        tmp_path / "no-hit-source.bam",
        [_segment("nonviral")],
    )
    outputs = _output_paths(tmp_path, "no-hit", normalized_suffix=".tsv.zst")

    result = _dispatch_lyra_postprocess([
        str(score_path),
        str(source_bam),
        "no hit sample",
        *(str(path) for path in outputs),
    ])

    assert result is None
    assert capsys.readouterr().out == ""
    assert _logical_tsv_bytes(outputs[0]) == (
        "\t".join(lyra.NORMALIZED_HEADER)
        + "\nno hit sample\tnonviral\t1\tSingle-end\t0.1\t0.1\t0.8\tNon-viral\n"
    ).encode("utf-8")
    summary = _read_summary(outputs[1])
    assert summary == dict(zip(
        lyra.SUMMARY_HEADER,
        ("no hit sample", "0.8", "1", "1", "1", "1", "1", "0", "0", "0", "1", "0"),
    ))
    source_header, _ = _read_bam_records(source_bam)
    output_header, output_records = _read_bam_records(outputs[2])
    assert output_header == source_header
    assert output_records == []
    _assert_summary_equations(summary, normalized_rows=1, output_bam_records=0)


def test_full_parser_malformed_native_input_has_stable_context(tmp_path, capsys):
    score_path = tmp_path / "malformed-scores.tsv"
    score_path.write_bytes(b"read_id\tscore\tcall\nmalformed\t0.9\tinvalid\n")
    source_bam = _write_bam(
        tmp_path / "malformed-source.bam",
        [_segment("malformed")],
    )
    outputs = _output_paths(tmp_path, "malformed")

    with pytest.raises(lyra.LyraInputError) as exc_info:
        _dispatch_lyra_postprocess([
            str(score_path),
            str(source_bam),
            "malformed sample",
            *(str(path) for path in outputs),
        ])

    error = exc_info.value
    assert error.category == "call"
    assert error.path == str(score_path)
    assert error.line_number == 2
    assert error.field == "call"
    assert capsys.readouterr().out == ""
    _assert_paths_absent(outputs)


def test_full_parser_score_bam_mismatch_has_stable_context(tmp_path, capsys):
    score_path = _write_score_table(
        tmp_path / "mismatch-scores.tsv",
        [("mismatch", "0.9", "1")],
    )
    source_bam = _write_bam(
        tmp_path / "mismatch-source.bam",
        [
            _segment("mismatch", flag=0x4 | 0x1 | 0x40),
            _segment("mismatch", flag=0x4 | 0x1 | 0x80),
        ],
    )
    outputs = _output_paths(tmp_path, "mismatch")

    with pytest.raises(lyra.LyraReconciliationError) as exc_info:
        _dispatch_lyra_postprocess([
            str(score_path),
            str(source_bam),
            "mismatch sample",
            *(str(path) for path in outputs),
        ])

    error = exc_info.value
    assert error.category == "record_count_mismatch"
    assert error.read_id == "mismatch"
    assert error.score_count == 1
    assert error.eligible_bam_count == 2
    assert error.score_line_numbers == (2,)
    assert capsys.readouterr().out == ""
    _assert_paths_absent(outputs)


def test_full_parser_existing_output_preserves_caller_bytes(tmp_path, capsys):
    score_path = _write_score_table(
        tmp_path / "existing-scores.tsv",
        [("read", "0.9", "1")],
    )
    source_bam = _write_bam(
        tmp_path / "existing-source.bam",
        [_segment("read")],
    )
    outputs = _output_paths(tmp_path, "existing")
    caller_bytes = b"caller-owned normalized artifact\n"
    outputs[0].write_bytes(caller_bytes)

    with pytest.raises(lyra.LyraPathError) as exc_info:
        _dispatch_lyra_postprocess([
            str(score_path),
            str(source_bam),
            "existing sample",
            *(str(path) for path in outputs),
        ])

    error = exc_info.value
    assert error.category == "output_exists"
    assert error.role == "normalized"
    assert error.path == str(outputs[0])
    assert outputs[0].read_bytes() == caller_bytes
    assert capsys.readouterr().out == ""
    _assert_paths_absent(outputs[1:])


def test_full_parser_real_command_does_not_load_lyra_runtime(tmp_path):
    assert FORBIDDEN_TOP_LEVEL_PREFIXES == (
        "lyra",
        "torch",
        "pytorch_lightning",
        "cuda",
        "cupy",
        "nvidia",
        "tensorflow",
        "transformers",
    )
    score_path = _write_score_table(
        tmp_path / "runtime-scores.tsv",
        [("runtime", "0.9", "0")],
    )
    source_bam = _write_bam(
        tmp_path / "runtime-source.bam",
        [_segment("runtime")],
    )
    outputs = _output_paths(tmp_path, "runtime")
    modules_path = tmp_path / "runtime-modules.json"

    child = subprocess.run(
        [
            sys.executable,
            "-c",
            CHILD_CODE,
            str(score_path),
            str(source_bam),
            *(str(path) for path in outputs),
            str(modules_path),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    assert child.stdout == ""
    assert json.loads(modules_path.read_text(encoding="utf-8")) == []
    assert _logical_tsv_bytes(outputs[0]) == (
        "\t".join(lyra.NORMALIZED_HEADER)
        + "\nruntime sample\truntime\t1\tSingle-end\t0.9\t0.9\t0.8\tViral\n"
    ).encode("utf-8")
    summary = _read_summary(outputs[1])
    assert summary == dict(zip(
        lyra.SUMMARY_HEADER,
        ("runtime sample", "0.8", "1", "1", "1", "1", "1", "0", "0", "1", "0", "1"),
    ))
    source_header, source_records = _read_bam_records(source_bam)
    output_header, output_records = _read_bam_records(outputs[2])
    assert output_header == source_header
    assert output_records == source_records
    _assert_summary_equations(summary, normalized_rows=1, output_bam_records=1)
