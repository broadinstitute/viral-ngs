import os
from unittest.mock import patch

import pandas as pd
import pytest

from viral_ngs.classify import virnucpro
from viral_ngs.core import misc as util_misc


@pytest.fixture
def virnucpro_tool():
    return virnucpro.VirNucPro(python_executable="python")


def _write(path, contents):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wt") as outf:
        outf.write(contents)


def test_predict_invokes_virnucpro_and_writes_both_tsv_outputs(tmp_path, virnucpro_tool):
    in_fasta = tmp_path / "sample.fa"
    in_fasta.write_text(">seq1\nACGT\n")
    output_dir = tmp_path / "virnucpro_out"
    out_predictions = tmp_path / "sample.virnucpro.predictions.tsv"
    out_highestscore = tmp_path / "sample.virnucpro.highestscore.tsv"

    predictions_contents = "Sequence_ID\tPrediction\tscore1\tscore2\nseq1\tvirus\t0.1\t0.9\n"
    highestscore_contents = "Modified_ID\tIs_Virus\tmax_score_0\tmax_score_1\nseq1\tTrue\t0.1\t0.9\n"

    def fake_execute(args):
        assert args[0] == "predict"
        merged_dir = output_dir / "sample_merged"
        _write(str(merged_dir / "prediction_results.txt"), predictions_contents)
        _write(str(merged_dir / "prediction_results_highestscore.csv"), highestscore_contents)

    with patch.object(virnucpro_tool, "execute", side_effect=fake_execute) as mock_execute:
        virnucpro_tool.predict(
            str(in_fasta),
            str(out_predictions),
            str(out_highestscore),
            model_type="300",
            model_path="/models/300_model.pth",
            expected_length=300,
            output_dir=str(output_dir),
            device="cuda",
            parallel=True,
            persistent_models=True,
            resume=True,
            keep_intermediate=True,
            batch_size=64,
            num_workers=2,
            esm_batch_size=1024,
            dnabert_batch_size=16,
            gpus="0,1",
            num_threads=3,
        )

    args = mock_execute.call_args.kwargs["args"]
    assert args[:2] == ["predict", str(in_fasta)]
    assert _contains_pair(args, "--output-dir", str(output_dir))
    assert _contains_pair(args, "--model-type", "300")
    assert _contains_pair(args, "--model-path", "/models/300_model.pth")
    assert _contains_pair(args, "--expected-length", "300")
    assert _contains_pair(args, "--device", "cuda")
    assert _contains_pair(args, "--batch-size", "64")
    assert _contains_pair(args, "--num-workers", "2")
    assert _contains_pair(args, "--esm-batch-size", "1024")
    assert _contains_pair(args, "--dnabert-batch-size", "16")
    assert _contains_pair(args, "--gpus", "0,1")
    assert _contains_pair(args, "--threads", str(util_misc.sanitize_thread_count(3)))
    assert "--force" in args
    assert "--no-progress" in args
    assert "--parallel" in args
    assert "--persistent-models" in args
    assert "--resume" in args
    assert "--keep-intermediate" in args
    assert out_predictions.read_text() == predictions_contents
    assert out_highestscore.read_text() == highestscore_contents


def test_predict_defaults_expected_length_from_model_type(tmp_path, virnucpro_tool):
    in_fasta = tmp_path / "sample.fa"
    in_fasta.write_text(">seq1\nACGT\n")
    output_dir = tmp_path / "virnucpro_out"

    def fake_execute(args):
        merged_dir = output_dir / "sample_merged"
        _write(str(merged_dir / "prediction_results.txt"), "Sequence_ID\tPrediction\tscore1\tscore2\n")
        _write(
            str(merged_dir / "prediction_results_highestscore.csv"),
            "Modified_ID\tIs_Virus\tmax_score_0\tmax_score_1\n",
        )

    with patch.object(virnucpro_tool, "execute", side_effect=fake_execute) as mock_execute:
        virnucpro_tool.predict(
            str(in_fasta),
            str(tmp_path / "predictions.tsv"),
            str(tmp_path / "highestscore.tsv"),
            model_type="500",
            output_dir=str(output_dir),
        )

    args = mock_execute.call_args.kwargs["args"]
    assert _contains_pair(args, "--expected-length", "500")


def test_predict_rejects_missing_input(tmp_path, virnucpro_tool):
    with pytest.raises(FileNotFoundError):
        virnucpro_tool.predict(
            str(tmp_path / "missing.fa"),
            str(tmp_path / "predictions.tsv"),
            str(tmp_path / "highestscore.tsv"),
        )


def test_predict_rejects_invalid_expected_length(tmp_path, virnucpro_tool):
    in_fasta = tmp_path / "sample.fa"
    in_fasta.write_text(">seq1\nACGT\n")

    with pytest.raises(ValueError, match="expected_length must be 300 or 500"):
        virnucpro_tool.predict(
            str(in_fasta),
            str(tmp_path / "predictions.tsv"),
            str(tmp_path / "highestscore.tsv"),
            expected_length=301,
        )


def test_predict_rejects_custom_model_without_model_path(tmp_path, virnucpro_tool):
    in_fasta = tmp_path / "sample.fa"
    in_fasta.write_text(">seq1\nACGT\n")

    with pytest.raises(ValueError, match="model_path is required"):
        virnucpro_tool.predict(
            str(in_fasta),
            str(tmp_path / "predictions.tsv"),
            str(tmp_path / "highestscore.tsv"),
            model_type="custom",
            expected_length=500,
        )


def test_predict_raises_when_expected_outputs_are_missing(tmp_path, virnucpro_tool):
    in_fasta = tmp_path / "sample.fa"
    in_fasta.write_text(">seq1\nACGT\n")

    with patch.object(virnucpro_tool, "execute", return_value=None):
        with pytest.raises(FileNotFoundError, match="Missing expected VirNucPro output"):
            virnucpro_tool.predict(
                str(in_fasta),
                str(tmp_path / "predictions.tsv"),
                str(tmp_path / "highestscore.tsv"),
                output_dir=str(tmp_path / "virnucpro_out"),
            )


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


def _contains_pair(args, key, value):
    return key in args and args[args.index(key) + 1] == value
