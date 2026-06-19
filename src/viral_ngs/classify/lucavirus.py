"""LucaVirus preflight and post-processing helpers."""

import csv
import os

from Bio import SeqIO

from viral_ngs.core import file


OUTPUT_HEADER = ["seq_id", "seq", "prob", "label_index", "label"]
LUCAVIRUS_INPUT_HEADER = ["seq_id", "seq_type", "seq"]
PREPARE_STATS_HEADER = ["n_sequences", "has_lucavirus_input"]
TASK_PROFILES = ("rdrp", "viral_capsid", "virus_ec4")


def validate_task_profile(task_profile):
    """Validate the configured LucaVirus task profile."""
    if task_profile not in TASK_PROFILES:
        raise ValueError(
            "Unsupported LucaVirus task profile '{}'. Expected one of: {}".format(
                task_profile, ", ".join(TASK_PROFILES)
            )
        )


def prepare_contigs(contigs_fasta, lucavirus_input_csv, stats_tsv, seq_type="prot"):
    """Convert input FASTA records to LucaVirus CSV and write preflight stats."""
    if not os.path.isfile(contigs_fasta):
        raise FileNotFoundError(contigs_fasta)
    if seq_type != "prot":
        raise ValueError(
            "Unsupported LucaVirus seq_type '{}'. Expected 'prot'.".format(seq_type)
        )

    _ensure_parent_dir(lucavirus_input_csv)
    _ensure_parent_dir(stats_tsv)

    n_sequences = 0
    with file.open_or_gzopen(lucavirus_input_csv, "wt", newline="") as out_fh:
        writer = csv.writer(out_fh, lineterminator="\n")
        writer.writerow(LUCAVIRUS_INPUT_HEADER)
        if not _is_blank_text_file(contigs_fasta):
            with file.open_or_gzopen(contigs_fasta, "rt") as in_fh:
                for record in SeqIO.parse(in_fh, "fasta"):
                    sequence = str(record.seq).strip()
                    if not sequence:
                        continue
                    writer.writerow([record.id, seq_type, sequence])
                    n_sequences += 1

    stats = {
        "n_sequences": n_sequences,
        "has_lucavirus_input": n_sequences > 0,
    }
    _write_prepare_stats(stats_tsv, stats)
    return stats


def write_empty_predictions(output_tsv):
    """Write a header-only LucaVirus prediction table."""
    _ensure_parent_dir(output_tsv)
    with file.open_or_gzopen(output_tsv, "wt", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_HEADER)


def normalize_output(input_tsv, output_tsv, task_profile="rdrp"):
    """Validate lucavirus-cuda output and copy it to the durable TSV artifact."""
    validate_task_profile(task_profile)
    if not os.path.isfile(input_tsv):
        raise FileNotFoundError(input_tsv)
    if os.path.getsize(input_tsv) == 0:
        raise ValueError("LucaVirus output is zero bytes: {}".format(input_tsv))

    rows = []
    with file.open_or_gzopen(input_tsv, "rt", newline="") as in_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("LucaVirus output is empty: {}".format(input_tsv))

        if header != OUTPUT_HEADER:
            raise ValueError(
                "Unexpected LucaVirus output header in {}. Expected {}, observed {}".format(
                    input_tsv, OUTPUT_HEADER, header
                )
            )

        for line_number, row in enumerate(reader, start=2):
            if not row or all(value == "" for value in row):
                raise ValueError(
                    "Blank LucaVirus output row at line {} in {}".format(
                        line_number, input_tsv
                    )
                )
            if len(row) != len(OUTPUT_HEADER):
                raise ValueError(
                    "Malformed LucaVirus output row at line {} in {}. "
                    "Expected {} columns, observed {}".format(
                        line_number, input_tsv, len(OUTPUT_HEADER), len(row)
                    )
                )
            rows.append(row)

    _ensure_parent_dir(output_tsv)
    with file.open_or_gzopen(output_tsv, "wt", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_HEADER)
        writer.writerows(rows)


def _write_prepare_stats(stats_tsv, stats):
    with file.open_or_gzopen(stats_tsv, "wt", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(PREPARE_STATS_HEADER)
        writer.writerow(
            [
                stats["n_sequences"],
                "true" if stats["has_lucavirus_input"] else "false",
            ]
        )


def _is_blank_text_file(path):
    with file.open_or_gzopen(path, "rt") as in_fh:
        for chunk in iter(lambda: in_fh.read(1024 * 1024), ""):
            if chunk.strip():
                return False
    return True


def _ensure_parent_dir(path):
    parent_dir = os.path.dirname(os.path.abspath(path))
    if parent_dir:
        file.mkdir_p(parent_dir)
