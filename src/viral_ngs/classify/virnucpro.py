"""VirNucPro viral nucleotide sequence classifier."""

import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

from viral_ngs import core
from viral_ngs.core import file
from viral_ngs.core import misc

log = logging.getLogger(__name__)


class VirNucPro(core.Tool):
    """Wrapper for the VirNucPro prediction CLI."""

    VALID_MODEL_TYPES = ("300", "500", "custom")
    VALID_EXPECTED_LENGTHS = (300, 500)
    PREDICTIONS_FILENAME = "prediction_results.txt"
    HIGHESTSCORE_FILENAME = "prediction_results_highestscore.csv"

    def __init__(self, install_methods=None, python_executable=None):
        self.python_executable = python_executable or sys.executable
        if not install_methods:
            install_methods = [
                core.PrexistingUnixCommand(
                    shutil.which("python") or self.python_executable,
                    require_executability=False,
                )
            ]
        super(VirNucPro, self).__init__(install_methods=install_methods)

    def version(self):
        """Return VirNucPro CLI version if available."""
        try:
            result = subprocess.run(
                [self.python_executable, "-m", "virnucpro", "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "unknown"
        return (result.stdout or result.stderr).strip() or "unknown"

    def execute(self, args=None):
        """Run a VirNucPro CLI command."""
        args = args or []
        cmd = [self.python_executable, "-m", "virnucpro"] + args
        log.debug("Calling VirNucPro: %s", " ".join(cmd))
        subprocess.check_call(cmd)

    def predict(
        self,
        in_fasta,
        out_predictions_tsv,
        out_highestscore_tsv,
        model_type="500",
        model_path=None,
        expected_length=None,
        output_dir=None,
        device="auto",
        parallel=False,
        persistent_models=False,
        resume=False,
        keep_intermediate=False,
        batch_size=None,
        num_workers=None,
        esm_batch_size=None,
        dnabert_batch_size=None,
        gpus=None,
        num_threads=None,
    ):
        """Classify nucleotide sequences and write stable TSV outputs."""
        if not os.path.isfile(in_fasta):
            raise FileNotFoundError(in_fasta)

        model_type = str(model_type)
        if model_type not in self.VALID_MODEL_TYPES:
            raise ValueError(
                "model_type must be one of {}".format(", ".join(self.VALID_MODEL_TYPES))
            )
        if model_type == "custom" and not model_path:
            raise ValueError("model_path is required when model_type='custom'")

        if expected_length is None:
            if model_type == "custom":
                raise ValueError("expected_length is required when model_type='custom'")
            expected_length = int(model_type)
        expected_length = int(expected_length)
        if expected_length not in self.VALID_EXPECTED_LENGTHS:
            raise ValueError("expected_length must be 300 or 500")

        cleanup_output_dir = output_dir is None
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="virnucpro-")

        try:
            args = [
                "predict",
                in_fasta,
                "--output-dir",
                output_dir,
                "--model-type",
                model_type,
                "--expected-length",
                expected_length,
                "--device",
                device,
                "--force",
                "--no-progress",
            ]
            if model_path:
                args.extend(["--model-path", model_path])
            if parallel:
                args.append("--parallel")
            if persistent_models:
                args.append("--persistent-models")
            if resume:
                args.append("--resume")
            if keep_intermediate:
                args.append("--keep-intermediate")
            if batch_size is not None:
                args.extend(["--batch-size", batch_size])
            if num_workers is not None:
                args.extend(["--num-workers", num_workers])
            if esm_batch_size is not None:
                args.extend(["--esm-batch-size", esm_batch_size])
            if dnabert_batch_size is not None:
                args.extend(["--dnabert-batch-size", dnabert_batch_size])
            if gpus is not None:
                args.extend(["--gpus", gpus])
            if num_threads is not None:
                args.extend(["--threads", misc.sanitize_thread_count(num_threads)])

            self.execute(args=[str(arg) for arg in args])

            predictions_path, highestscore_path = self._result_paths(output_dir, in_fasta)
            _copy_tsv(predictions_path, out_predictions_tsv)
            _copy_tsv(highestscore_path, out_highestscore_tsv)
        finally:
            if cleanup_output_dir:
                shutil.rmtree(output_dir, ignore_errors=True)

    @classmethod
    def _result_paths(cls, output_dir, in_fasta):
        merged_dir = os.path.join(
            output_dir,
            "{}_merged".format(Path(in_fasta).stem),
        )
        predictions_path = os.path.join(merged_dir, cls.PREDICTIONS_FILENAME)
        highestscore_path = os.path.join(merged_dir, cls.HIGHESTSCORE_FILENAME)
        missing = [
            path
            for path in (predictions_path, highestscore_path)
            if not os.path.isfile(path)
        ]
        if missing:
            raise FileNotFoundError(
                "Missing expected VirNucPro output file(s): {}".format(
                    ", ".join(missing)
                )
            )
        return predictions_path, highestscore_path


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


def _copy_tsv(src, dest):
    _ensure_parent_dir(dest)
    if os.path.abspath(src) == os.path.abspath(dest):
        return
    shutil.copyfile(src, dest)


def _ensure_parent_dir(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        file.mkdir_p(parent)


def _natural_sort_key(value):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]
