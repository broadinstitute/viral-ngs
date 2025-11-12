import csv
import glob
import json
import logging
import os
import shutil
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import tools
import tools.samtools
import util.file

TOOL_NAME = 'splitcode'

log = logging.getLogger(__name__)

class SplitCodeTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(SplitCodeTool, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '--version']).decode('UTF-8').strip()

    def execute(    self,
                    n_fastqs,
                    threads,
                    config_file,
                    keep_file,
                    unassigned_r1,
                    unassigned_r2,
                    summary_stats,
                    r1,
                    r2,
                    splitcode_opts=None,
                    gzip_output=False,
                    keep_r1_r2_suffixes=True,
                    predemux_r1_trim_5prime_num_bp=None,
                    predemux_r1_trim_3prime_num_bp=None,
                    predemux_r2_trim_5prime_num_bp=None,
                    predemux_r2_trim_3prime_num_bp=None
                ):
        """
        Execute the splitcode command with the provided parameters.

        :param n_fastqs: Number of FASTQ files (e.g., 2 for paired-end reads).
        :param threads: Number of threads to use.
        :param config_file: Path to the configuration file.
        :param keep_file: Path to the keep file.
        :param unassigned_r1: Path for unassigned R1 reads.
        :param unassigned_r2: Path for unassigned R2 reads.
        :param summary_stats: Path to the summary statistics file.
        :param r1: Input FASTQ file for R1 reads.
        :param r2: Input FASTQ file for R2 reads.
        :param splitcode_opts: additional parameters to pass to splitcode
        """
        splitcode_opts = splitcode_opts or []

        # documentation: 
        #   https://splitcode.readthedocs.io/en/latest/reference_guide.html#command-line-config-optional

        tool_cmd = [
            self.install_and_get_path(),
            f"--nFastqs={n_fastqs}",
            f"-t", str(threads),
            f"-c", config_file,
            f"--keep={keep_file}",
            f"--unassigned={unassigned_r1},{unassigned_r2}",
            f"--summary", summary_stats
        ]
        if gzip_output and "--gzip" not in splitcode_opts:
            tool_cmd.append("--gzip")
        if keep_r1_r2_suffixes and "--keep-r1-r2" not in splitcode_opts:
            tool_cmd.append("--keep-r1-r2")
        if predemux_r1_trim_3prime_num_bp or predemux_r2_trim_3prime_num_bp:
            r1_3prime_bp_to_remove = predemux_r1_trim_3prime_num_bp or 0
            r2_3prime_bp_to_remove = predemux_r2_trim_3prime_num_bp or 0
            tool_cmd.extend(["--trim-3",
                                f"{r1_3prime_bp_to_remove},{r2_3prime_bp_to_remove}"])
        if predemux_r1_trim_5prime_num_bp or predemux_r2_trim_5prime_num_bp:
            r1_5prime_bp_to_remove = predemux_r1_trim_5prime_num_bp or 0
            r2_5prime_bp_to_remove = predemux_r2_trim_5prime_num_bp or 0
            tool_cmd.extend(["--trim-5", 
                                f"{r1_5prime_bp_to_remove},{r2_5prime_bp_to_remove}"])

        tool_cmd.extend(splitcode_opts)
        tool_cmd.extend([r1,r2])

        log.debug('Running splitcode with command: %s', ' '.join(tool_cmd))
        try:
            return subprocess.check_call(tool_cmd)
        except subprocess.CalledProcessError as e:
            log.error(f'Splitcode failed with return code {e.returncode}')
            log.error(f'Command was: {" ".join(tool_cmd)}')
            raise

    def check_installation(self):
        """Ensure the tool is properly installed."""
        if not shutil.which(TOOL_NAME):
            raise FileNotFoundError(f"{TOOL_NAME} is not installed or not in PATH.")

    def run_splitcode(self, **kwargs):
        """
        Wrapper method to execute the splitcode command with named parameters.

        Expected kwargs:
        - n_fastqs, threads, config_file, keep_file, unassigned_r1, unassigned_r2, summary_stats, r1, r2
        """
        self.execute(
            kwargs.get('n_fastqs', 2),
            kwargs.get('threads', 1),
            kwargs['config_file'],
            kwargs['keep_file'],
            kwargs['unassigned_r1'],
            kwargs['unassigned_r2'],
            kwargs['summary_stats'],
            kwargs['r1'],
            kwargs['r2']
        )


# ==================
# ***  Splitcode helper functions   ***
# ==================


def create_splitcode_lookup_table(sample_sheet_or_dataframe, csv_out, unmatched_name, pool_ids=None, append_run_id=None, check_sample_sheet_consistency=False):
    """
    Create a lookup table (LUT) consolidating splitcode demux results with sample metadata.

    After splitcode demultiplexes pooled samples by inner barcodes (barcode_3), this function
    reads the splitcode summary JSON files, joins them with the sample sheet, and creates a
    unified CSV mapping samples to barcodes to read counts. This LUT is used by downstream
    plotting and metrics functions.

    The function performs 3-barcode demux integration by:
    1. Loading sample sheet with all barcode mappings (barcode_1, barcode_2, barcode_3)
    2. Building pool identifiers from outer barcode pairs (barcode_1 + barcode_2)
    3. Loading splitcode summary JSONs for each pool
    4. Extracting read counts at Hamming distance 0 (perfect match) and 1 (1-mismatch)
    5. Joining splitcode counts with sample metadata based on run identifiers
    6. Creating unmatched read entries for reads that didn't match any barcode
    7. Outputting a CSV with unified sample-to-barcode-to-count mappings

    Parameters
    ----------
    sample_sheet_or_dataframe : str or pd.DataFrame
        Either:
        - Path to TSV file with barcode mappings (str), OR
        - DataFrame from SampleSheet.inner_demux_mapper() (pd.DataFrame)

        If DataFrame, must already contain columns:
        - barcode_1, barcode_2, barcode_3, sample, library_id_per_sample
        - run: Sample run identifier (sample.lLibrary[.FlowcellID])
        - muxed_run: Pool identifier (barcode_1-barcode_2.lLibrary[.FlowcellID])

        If str (file path), must be TSV with columns:
        - barcode_1, barcode_2, barcode_3, sample, library_id_per_sample
    csv_out : str
        Output path for the lookup table CSV file
    unmatched_name : str
        Name prefix for unmatched/unassigned reads (e.g., "Unmatched")
    pool_ids : list, optional
        List of pool IDs to process. If None/empty, processes all pools found.
    append_run_id : str, optional
        Suffix to append to run identifiers (typically flowcell ID).
        Only used when sample_sheet_or_dataframe is a file path (str).
        Ignored when a DataFrame is provided (which already has run/muxed_run).
    check_sample_sheet_consistency : bool, optional
        If True, validates sample sheet for duplicate barcode combinations

    Returns
    -------
    str
        Path to the output CSV file containing the lookup table

    Output CSV Schema
    -----------------
    The output CSV contains columns:
    - sample: Sample identifier
    - library_id: Library/pool identifier
    - barcode_1: Outer barcode (i7 index)
    - barcode_2: Outer barcode (i5 index)
    - inline_barcode: Inner barcode (barcode_3)
    - run: Sample run identifier (sample.lLibrary[.FlowcellID])
    - muxed_pool: Pool identifier (barcode_1-barcode_2.lLibrary)
    - num_reads_hdistance0: Read count with perfect barcode match
    - num_reads_hdistance1: Read count with 1-mismatch to barcode
    - num_reads_total: Total reads (hdistance0 + hdistance1)

    Notes
    -----
    - Expects splitcode summary JSON files named "{pool}_summary.json" in csv_out directory
    - Validates that (barcode_1, barcode_2) pairs don't have duplicate barcode_3 values
    - Handles pools with 0 reads gracefully by creating empty metrics
    - Unmatched reads get barcode_3 set to all "N"s matching expected barcode length
    """
    pool_ids = pool_ids or []

    outDir=os.path.dirname(csv_out)

    # Load or reuse barcode data
    # If we receive a DataFrame (from SampleSheet.inner_demux_mapper()), use it directly
    # This ensures pool IDs match exactly what was used to create summary JSON files
    if isinstance(sample_sheet_or_dataframe, pd.DataFrame):
        log.debug("Using provided DataFrame for barcode lookup table")
        barcodes_df = sample_sheet_or_dataframe.copy()
        # Reset index to make 'sample' a regular column (inner_demux_mapper sets it as index)
        barcodes_df.reset_index(inplace=True)

        # Rename 'muxed_run' to 'muxed_pool' for consistency with legacy code
        # The inner_demux_mapper() creates 'muxed_run' column, but this function expects 'muxed_pool'
        if 'muxed_run' in barcodes_df.columns and 'muxed_pool' not in barcodes_df.columns:
            barcodes_df['muxed_pool'] = barcodes_df['muxed_run']
    else:
        # Legacy path: read from file and build columns
        log.debug("Reading sample sheet from file: %s", sample_sheet_or_dataframe)
        barcodes_df = pd.read_csv(sample_sheet_or_dataframe, sep="\t", dtype=str)

        # Build "run" column = "<sample>.l<library_id_per_sample>[.<append_run_id>]"
        def build_sample_library_id_string(row):
            sample_val = row.get("sample", "")
            lib_val    = row.get("library_id_per_sample", "")
            run_str = f"{sample_val}.l{lib_val}"
            if append_run_id:
                run_str += f".{append_run_id}"
            return run_str
        barcodes_df["run"] = barcodes_df.apply(build_sample_library_id_string, axis=1)

        def build_muxed_run_string(row):
            b1 = row.get("barcode_1", "")
            b2 = row.get("barcode_2", None)
            if b2 and b2.strip():
                muxed_sample = f"{b1}-{b2}"
            else:
                muxed_sample = b1

            lib_val = row.get("library_id_per_sample", "")
            muxed_run_str = f"{muxed_sample}.l{lib_val}"
            if append_run_id:
                muxed_run_str += f".{append_run_id}"
            return muxed_run_str
        barcodes_df["muxed_pool"] = barcodes_df.apply(build_muxed_run_string, axis=1)

    df_csv_out = f"{csv_out}"

    i7_barcodes = list(set(barcodes_df["barcode_1"].values))
    barcodes = sorted(list(set(barcodes_df["barcode_3"].values)))

    def duplication_check(df, primary_cols, secondary_col, error_message_header=None, error_message=None):
        default_error_message_header = "Error: More than one '{column_to_check_for_duplicates}' value present for distinct combinations of the columns {affected_column_names}:"
        default_error_message        = "'{duplicated_values}' appears {duplicate_count} times, for {affected_values}."

        error_message_header = error_message_header or default_error_message_header
        error_message        = error_message        or default_error_message

        # Check which placeholders are actually present in the template
        err_header_has_dup_check_col         = '{column_to_check_for_duplicates}' in error_message_header
        err_header_has_affected_column_names = '{affected_column_names}' in error_message_header

        err_msg_has_affected_values          = '{affected_values}' in error_message
        err_msg_has_dup_val                  = '{duplicated_values}' in error_message
        err_msg_has_count                    = '{duplicate_count}' in error_message


        # To store all generated error messages
        errors = []

        # Group the dataframe by the primary columns
        grouped = df.groupby(primary_cols, dropna=False)

        error_header_out=None
        for group_key, subdf in grouped:
            # If there's only one primary column, group_key is a single value
            # Otherwise, it's a tuple of values
            if isinstance(group_key, tuple):
                # E.g.: ("North", "A")
                group_str = ','.join(
                    f"{col}={val}" for col, val in zip(primary_cols, group_key)
                )
            else:
                # E.g.: "North" when there's only one primary column
                group_str = f"{primary_cols[0]}={group_key}"

            # Count occurrences of each value in secondary_col
            value_counts = subdf[secondary_col].value_counts(dropna=True)
            # Identify duplicates: values with count > 1
            duplicates = value_counts[value_counts > 1]


            # Build a separate message for each duplicated value
            for val, count_ in duplicates.items():
                msg_dict        = {}
                header_msg_dict = {}

                if err_msg_has_affected_values:
                    msg_dict['affected_values'] = group_str
                if err_msg_has_dup_val:
                    msg_dict['duplicated_values'] = val
                if err_msg_has_count:
                    msg_dict['duplicate_count'] = count_
                if err_header_has_dup_check_col:
                    header_msg_dict['column_to_check_for_duplicates'] = secondary_col
                if err_header_has_affected_column_names:
                    header_msg_dict['affected_column_names'] = f"{'+'.join([f'{chr(39)+c+chr(39)}' for c in primary_cols])}"

                if error_header_out is None:
                    error_header_out = error_message_header.format(**header_msg_dict)
                message = error_message.format(**msg_dict)
                errors.append(message)
        return (error_header_out, errors)


    duplication_check_conditions = [
        {
            "columns": ["barcode_1","barcode_2"],
            "column_to_check_for_duplicates": "barcode_3",
        }
    ]

    for dup_condition in duplication_check_conditions:

        problem_header, problems = duplication_check( barcodes_df,
                                                        dup_condition["columns"],
                                                        dup_condition["column_to_check_for_duplicates"],
                                                        dup_condition.get("error_message_header", None),
                                                        dup_condition.get("error_message", None) )
        problem_found = False
        if len(problems):
            log.warning(problem_header)
            for problem in problems:
                log.warning(f"\t{problem}")
            problem_found = True
        if problem_found:
            raise ValueError("Problem(s) found in sample sheet; see above for details.")

    pool_dfs      = []
    unmatched_dfs = []

    for pool in barcodes_df["muxed_pool"].unique():
        # Get and load splitcode stats report json
        # Use the full pool name (including run suffix) to match the JSON filename created by splitcode
        pool_for_file_lookup = pool

        # Try to find and load the splitcode summary JSON file
        # Add robust error handling since missing/misplaced JSON files are a common issue
        try:
            summary_pattern = f"{outDir}/{pool_for_file_lookup}_summary.json"
            matching_files = glob.glob(summary_pattern)

            if not matching_files:
                # JSON file not found - list directory contents for debugging
                log.error(f"Splitcode summary JSON not found for pool '{pool_for_file_lookup}'")
                log.error(f"  Expected pattern: {summary_pattern}")
                log.error(f"  Searching in directory: {outDir}")

                # List all files in the output directory to help debug
                try:
                    dir_contents = os.listdir(outDir)
                    log.error(f"  Directory contents ({len(dir_contents)} files):")
                    # List JSON files first (most relevant)
                    json_files = [f for f in dir_contents if f.endswith('.json')]
                    if json_files:
                        log.error(f"    JSON files found ({len(json_files)}):")
                        for f in sorted(json_files):
                            log.error(f"      - {f}")
                    else:
                        log.error(f"    No JSON files found in directory")

                    # List first 20 other files for context
                    other_files = [f for f in dir_contents if not f.endswith('.json')]
                    if other_files:
                        log.error(f"    Other files (showing first 20 of {len(other_files)}):")
                        for f in sorted(other_files)[:20]:
                            log.error(f"      - {f}")
                except OSError as list_err:
                    log.error(f"  Could not list directory contents: {list_err}")

                raise FileNotFoundError(
                    f"Splitcode summary JSON not found for pool '{pool_for_file_lookup}'. "
                    f"Expected file: {summary_pattern}. "
                    f"Check logs above for directory contents."
                )

            splitcode_summary_file = matching_files[0]

            # Warn if multiple matches found (shouldn't happen but good to catch)
            if len(matching_files) > 1:
                log.warning(f"Multiple summary JSON files match pattern '{summary_pattern}':")
                for f in matching_files:
                    log.warning(f"  - {f}")
                log.warning(f"Using first match: {splitcode_summary_file}")

            log.debug(f"Loading splitcode summary from: {splitcode_summary_file}")

            with open(splitcode_summary_file, "r") as f:
                splitcode_summary = json.load(f)

        except (FileNotFoundError, IndexError) as e:
            # Re-raise with more context (directory listing already logged above)
            raise
        except json.JSONDecodeError as e:
            log.error(f"Failed to parse JSON from {splitcode_summary_file}")
            log.error(f"  JSON decode error: {e}")
            # Try to show first few lines of the file for debugging
            try:
                with open(splitcode_summary_file, "r") as f:
                    lines = f.readlines()
                    log.error(f"  File contents (first 10 lines):")
                    for i, line in enumerate(lines[:10], 1):
                        log.error(f"    {i}: {line.rstrip()}")
                    if len(lines) > 10:
                        log.error(f"    ... ({len(lines) - 10} more lines)")
            except Exception as read_err:
                log.error(f"  Could not read file for debugging: {read_err}")
            raise
        except Exception as e:
            log.error(f"Unexpected error loading splitcode summary for pool '{pool_for_file_lookup}'")
            log.error(f"  File: {splitcode_summary_file if 'splitcode_summary_file' in locals() else 'not determined'}")
            log.error(f"  Error type: {type(e).__name__}")
            log.error(f"  Error message: {e}")
            raise

        samplesheet_rows_for_pool_df = barcodes_df[barcodes_df["muxed_pool"] == pool]

        # Parse splitcode summary JSON
        # IMPORTANT: The tag_qc array has MULTIPLE entries per barcode tag!
        # Each barcode appears once for each hamming distance level (0, 1, 2, 3).
        # Example tag_qc structure:
        #   [
        #     {"tag": "Sample1_R1", "distance": 0, "count": 5},   # Perfect matches
        #     {"tag": "Sample1_R1", "distance": 1, "count": 2},   # 1 mismatch
        #     {"tag": "Sample1_R1", "distance": 2, "count": 0},   # 2 mismatches
        #     {"tag": "Sample1_R1", "distance": 3, "count": 0},   # 3 mismatches
        #     {"tag": "Sample2_R1", "distance": 0, "count": 3},
        #     ...
        #   ]
        #
        # Here we filter to distance=0 (perfect matches) and distance=1 (1-mismatch) separately
        # for downstream metrics and QC analysis.
        if len(splitcode_summary.get("tag_qc", [])) > 0:
            splitcode_summary_df = pd.DataFrame.from_records(splitcode_summary["tag_qc"])
            # Convert only the tag column to string, keep count/distance as numeric
            splitcode_summary_df['tag'] = splitcode_summary_df['tag'].astype(str)

            splitcode_summary_df['run'] = splitcode_summary_df['tag'].copy()
            splitcode_summary_df['run'] = splitcode_summary_df['run'].str.removesuffix('_R1')

            # Extract perfect matches (hamming distance = 0)
            splitcode_summary_df_h0_df  = splitcode_summary_df[splitcode_summary_df["distance"] == 0]
            # Extract 1-mismatch reads (hamming distance = 1)
            splitcode_summary_df_h1_df  = splitcode_summary_df[splitcode_summary_df["distance"] == 1].copy()

            splitcode_summary_df_h1_df  = splitcode_summary_df_h1_df.rename(columns={"count": "count_h1"})

            samplesheet_rows_for_pool_hx_df = samplesheet_rows_for_pool_df.join(
                                                splitcode_summary_df_h0_df.set_index('run'),
                                                on='run')

            samplesheet_rows_for_pool_hx_df = pd.merge(samplesheet_rows_for_pool_hx_df,
                                                splitcode_summary_df_h1_df[['run','count_h1']].rename(columns={'run':'run_h1'}),
                                                 left_on  = 'run',
                                                 right_on = 'run_h1',
                                                 how      = 'left')
            samplesheet_rows_for_pool_hx_df = samplesheet_rows_for_pool_hx_df.drop(columns=['run_h1','distance','tag']) # dropping distance since we've added a col with different distance (as indicated by _h1 suffix)
            # fil NA values in 'count_h1' and cast to int
            samplesheet_rows_for_pool_hx_df["count_h1"] = samplesheet_rows_for_pool_hx_df["count_h1"].fillna(0).astype(int)
        else:
            # No reads were processed by splitcode for this pool
            # Create a dataframe with the expected schema but all counts set to 0
            log.warning(f"Pool {pool} has 0 reads processed by splitcode. Creating empty metrics.")
            samplesheet_rows_for_pool_hx_df = samplesheet_rows_for_pool_df.copy()
            samplesheet_rows_for_pool_hx_df['count'] = 0
            samplesheet_rows_for_pool_hx_df['count_h1'] = 0

        pool_dfs.append(samplesheet_rows_for_pool_hx_df)

        unmatched_dict = {
            "sample"                : f"{unmatched_name}.{pool}",
            "library_id_per_sample" : list(set(samplesheet_rows_for_pool_hx_df["library_id_per_sample"]))[0],
            "run"                   : f"{unmatched_name}.{pool}",
            "muxed_pool"            : pool,
            "count"                 : splitcode_summary["n_processed"] - splitcode_summary["n_assigned"],
            "count_h1"              : 0,
            "barcode_1"             : list(samplesheet_rows_for_pool_hx_df["barcode_1"])[0],
            "barcode_2"             : list(samplesheet_rows_for_pool_hx_df["barcode_2"])[0],
            "barcode_3"             : "N" * len(list(samplesheet_rows_for_pool_hx_df["barcode_3"])[0]),
        }
        unmatched_df = pd.DataFrame.from_dict([unmatched_dict], dtype=str)
        unmatched_dfs.append(unmatched_df)

    all_pools_and_unmatched_df = pd.concat(pool_dfs + unmatched_dfs, ignore_index=True)

    df_lut = all_pools_and_unmatched_df.copy()
    # rename columns to values expected by downstream plotting code
    df_lut = df_lut.rename(columns={
                                    "count": "num_reads_hdistance0",
                                    "count_h1": "num_reads_hdistance1",
                                    "barcode_3":"inline_barcode",
                                    "library_id_per_sample":"library_id"
                                    })

    df_lut["num_reads_total"] = (
        df_lut["num_reads_hdistance0"] + df_lut["num_reads_hdistance1"]
    )

    # TO-DO: Adjust the table outputted by this function to match Picard output

    df_lut.to_csv(df_csv_out, index=False)

    return df_csv_out

def plot_read_counts(df_csv_path, outDir):
    df_lut = pd.read_csv(df_csv_path, dtype=str)
    # Convert numeric columns explicitly
    numeric_cols = ['num_reads_total', 'num_reads_hdistance0', 'num_reads_hdistance1']
    for col in numeric_cols:
        if col in df_lut.columns:
            df_lut[col] = pd.to_numeric(df_lut[col], errors='coerce').fillna(0)

    fig, axs = plt.subplots(figsize=(10, 10), nrows=3, sharex=True)
    fontsize = 14

    df_grouped = (
        df_lut.groupby(["inline_barcode", "library_id"])["num_reads_total"]
        .sum()
        .unstack(fill_value=0)
    )
    df_grouped_fracs = df_grouped.div(df_grouped.sum(axis=0), axis=1)

    bar_width = 0.2
    bar_positions = np.arange(len(df_grouped))

    # Define colors
    unique_library_ids = df_lut["library_id"].nunique()
    tab20_colors = plt.cm.tab20.colors
    pool_colors = (tab20_colors * (unique_library_ids // 20 + 1))[:unique_library_ids]

    for i, pool in enumerate(df_grouped.columns):
        axs[0].bar(
            bar_positions + i * bar_width,
            df_grouped[pool],
            width=bar_width,
            label=f'{pool.split("_")[-1]}',
            color=pool_colors[i],
        )
        axs[1].bar(
            bar_positions + i * bar_width,
            df_grouped[pool],
            width=bar_width,
            color=pool_colors[i],
        )
        axs[2].bar(
            bar_positions + i * bar_width,
            df_grouped_fracs[pool],
            width=bar_width,
            color=pool_colors[i],
        )

    # Only set log scale if there are positive values
    if df_grouped.values.max() > 0:
        axs[1].set_yscale("log")
    else:
        log.warning("All read counts are zero; skipping log scale for plot")

    for ax in axs[:2]:
        ax.set_ylabel("# Reads", fontsize=fontsize)
    axs[2].set_ylabel("Fraction of Reads", fontsize=fontsize)

    for ax in axs:
        ax.set_xticks(bar_positions + 1.5 * bar_width)
        ax.set_xticklabels(df_grouped.index, rotation=45, ha="right")
        ax.tick_params(axis="both", labelsize=fontsize - 2)
        ax.grid(True, axis="y", linestyle="--", alpha=0.7)
        ax.margins(x=0.01)

    axs[0].legend(title="Pool", fontsize=fontsize - 2, title_fontsize=fontsize)
    axs[0].set_title("# Reads per inline barcode", fontsize=fontsize)
    axs[2].set_xlabel("Inline Barcode", fontsize=fontsize)

    plt.tight_layout()

    fig.savefig(f"{outDir}/reads_per_pool.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{outDir}/reads_per_pool.png", bbox_inches="tight", dpi=300)


def plot_sorted_curve(df_csv_path, out_dir, unmatched_name, out_basename=None):
    out_basename = out_basename or "reads_per_pool_sorted_curve"

    df_lut = pd.read_csv(df_csv_path, dtype=str)
    # Convert numeric columns explicitly
    numeric_cols = ['num_reads_total', 'num_reads_hdistance0', 'num_reads_hdistance1']
    for col in numeric_cols:
        if col in df_lut.columns:
            df_lut[col] = pd.to_numeric(df_lut[col], errors='coerce').fillna(0)

    log.debug(f"Reading in metrics file for plotting barcode read counts per pool: {df_csv_path}")
    log.debug(f"unmatched_name: {unmatched_name}")

    fig, axs = plt.subplots(figsize=(10, 10), nrows=4, sharex=True)
    fontsize = 14

    # Define colors
    unique_library_ids = df_lut["library_id"].nunique()
    log.debug(f"Number of distinct library_id values (pools) present: {unique_library_ids}")
    log.debug(f"library_id values (pools): {", ".join(sorted(list(set(df_lut['library_id'].astype(str)))))}")

    tab20_colors = plt.cm.tab20.colors
    pool_colors  = (tab20_colors * (unique_library_ids // 20 + 1))[:unique_library_ids]

    for i, pool in enumerate(sorted(df_lut["library_id"].unique())):
        log.debug(f"Processing read counts to plot for library_id (pool): {pool}")
        # pool_metrics = df_lut[
        #     (df_lut["library_id"] == str(pool))
        #     & (df_lut["inline_barcode"] != unmatched_name)
        # ]
        pool_metrics = df_lut[
            (df_lut["library_id"] == pool)
            & (~df_lut['inline_barcode'].str.match(r'^N+$', na=False))
        ]
        num_reads = pool_metrics["num_reads_total"]
        num_reads = sorted(num_reads, reverse=True)
        for ax in axs[:2]:
            ax.scatter(
                np.arange(len(num_reads)),
                num_reads,
                label=f'{pool.split("_")[-1]}',
                color=pool_colors[i],
            )
            ax.plot(np.arange(len(num_reads)), num_reads, color=pool_colors[i])

        # Calculate total and fractions (handle division by zero)
        total_reads = sum(num_reads)
        if total_reads > 0:
            fractions = np.array([read / total_reads for read in num_reads])
        else:
            fractions = np.zeros(len(num_reads))
            log.warning(f"Pool {pool} has 0 total reads; setting all fractions to 0")

        for ax in axs[2:]:
            ax.scatter(np.arange(len(fractions)), fractions * 100, color=pool_colors[i])
            ax.plot(np.arange(len(fractions)), fractions * 100, color=pool_colors[i])

    axs[1].set_yscale("symlog")
    axs[1].set_ylim(bottom=0)
    axs[3].set_ylim(bottom=0, top=0.5)
    axs[3].set_xlabel("Inline Barcode", fontsize=fontsize)
    axs[0].set_title("Reads per Inline Barcode (Sorted Curve)", fontsize=fontsize)
    axs[0].legend(title="Pool", fontsize=fontsize - 2, title_fontsize=fontsize)

    for ax in axs[:2]:
        ax.set_ylabel("# Reads", fontsize=fontsize)

    for ax in axs[2:]:
        ax.set_ylabel("% Reads", fontsize=fontsize)

    for ax in axs:
        ax.tick_params(axis="both", labelsize=fontsize - 2)
        ax.grid(True, axis="y", linestyle="--", alpha=0.7)
        ax.margins(x=0.01)

    plt.tight_layout()

    fig.savefig(
        f"{out_dir}/{out_basename}.pdf", bbox_inches="tight", dpi=300
    )
    fig.savefig(
        f"{out_dir}/{out_basename}.png", bbox_inches="tight", dpi=300
    )

# this function is called in new processes
# and must remain at the top level (global scope) of this file
# to be picklable and thus compatible
# with concurrent.futures.ProcessPoolExecutor()
# see: https://docs.python.org/3/library/concurrent.futures.html#processpoolexecutor
#      https://stackoverflow.com/a/72776044
def run_splitcode_on_pool(  pool_id,
                            pool_bam_file,
                            splitcode_config,
                            splitcode_keepfile,
                            out_demux_dir_path,
                            unmatched_name         = "unmatched",
                            threads_per_worker     = None,
                            out_dir_path           = None,
                            out_demux_dir_path_tmp = None,
                            string_to_log          = None,
                            predemux_r1_trim_5prime_num_bp=None,
                            predemux_r1_trim_3prime_num_bp=None,
                            predemux_r2_trim_5prime_num_bp=None,
                            predemux_r2_trim_3prime_num_bp=None
                            ):
    """
    Execute splitcode demultiplexing on a single BAM pool.

    This function converts the input BAM to FASTQ, runs splitcode to demultiplex based on
    inline barcodes, and generates a summary JSON with barcode matching statistics.

    Args:
        pool_id: Identifier for this pool (used in output filenames)
        pool_bam_file: Path to input BAM file containing multiplexed reads
        splitcode_config: Path to splitcode config file (TSV format with columns:
                         tag, id, locations, distance, left, right)
        splitcode_keepfile: Path to keep file (TSV: barcode_id<tab>output_prefix)
        out_demux_dir_path: Directory for demultiplexed output FASTQs
        unmatched_name: Prefix for unmatched/unassigned reads output
        out_dir_path: Directory for summary JSON (defaults to out_demux_dir_path)
        out_demux_dir_path_tmp: Temporary directory for intermediate files

    Returns:
        tuple: (return_code, pool_id) where return_code is 0 on success

    Output files:
        - Summary JSON: {out_dir_path}/{pool_id}_summary.json
          Contains 'tag_qc' array with barcode matching statistics. Each barcode has
          multiple entries (one per hamming distance level 0-3) with 'tag', 'distance',
          and 'count' fields.
        - Demuxed FASTQs: {output_prefix}_R1.fastq, {output_prefix}_R2.fastq
          (output_prefix comes from keepfile, one pair per barcode)
        - Unmatched FASTQs: {unmatched_name}.{pool_id}_R1.fastq, {unmatched_name}.{pool_id}_R2.fastq

    Notes:
        - Splitcode config locations format: FILE_NUMBER:START_BP:END_BP
          Example: "0:0:8" means file 0 (R1), positions 0-8 (8bp barcode)
        - The 'id' column in config file must match the first column in keepfile
        - By convention, we use f"{sample_library_id}_R1" as the ID
        - tag_qc in summary JSON has multiple entries per tag (one for each distance level),
          so downstream code must sum counts across all distance levels if needed

    See test/unit/test_illumina_splitcode.py for detailed examples and validation of
    splitcode behavior assumptions.
    """
    with tools.samtools.SamtoolsTool().bam2fq_tmp(pool_bam_file) as (fqin1, fqin2):
        n_fastqs = 2
        splitcode_tool = tools.splitcode.SplitCodeTool()

        # output args
        unmapped_r1   = f"{out_demux_dir_path_tmp}/{unmatched_name}.{pool_id}_R1.fastq"
        unmapped_r2   = f"{out_demux_dir_path_tmp}/{unmatched_name}.{pool_id}_R2.fastq"
        # write the stats to the output directory rather than tmp
        summary_stats = f"{out_dir_path}/{pool_id}_summary.json"

        # Dummy output files for --output parameter
        # Note: --no-output flag prevents --keep files from being written, so we must
        # provide --output even though we don't use these files
        dummy_output_r1 = f"{out_demux_dir_path_tmp}/_dummy_main_output_{pool_id}_R1.fastq"
        dummy_output_r2 = f"{out_demux_dir_path_tmp}/_dummy_main_output_{pool_id}_R2.fastq"

        # Optional: Pass 'predemux_r2_trim_3prime_num_bp':8 for '--trim-3 0,8'
        # to remove 8 bases from end of R2 read
        # (note: trimming happens before demuxing)
        # valid options:
        #   predemux_r1_trim_5prime_num_bp
        #   predemux_r1_trim_3prime_num_bp
        #   predemux_r2_trim_5prime_num_bp
        #   predemux_r2_trim_3prime_num_bp
        # see splitcode.py wrapper for more details, and also:
        #   https://splitcode.readthedocs.io/en/latest/reference_guide.html#command-line-config-optional
        splitcode_kwargs={
            "n_fastqs"                       : n_fastqs,
            "threads"                        : threads_per_worker,
            "config_file"                    : splitcode_config,
            "keep_file"                      : splitcode_keepfile,
            "unassigned_r1"                  : unmapped_r1,
            "unassigned_r2"                  : unmapped_r2,
            "summary_stats"                  : summary_stats,
            "r1"                             : fqin1,
            "r2"                             : fqin2,
            "predemux_r1_trim_5prime_num_bp" : predemux_r1_trim_5prime_num_bp,
            "predemux_r1_trim_3prime_num_bp" : predemux_r1_trim_3prime_num_bp,
            "predemux_r2_trim_5prime_num_bp" : predemux_r2_trim_5prime_num_bp,
            "predemux_r2_trim_3prime_num_bp" : predemux_r2_trim_3prime_num_bp,
            "keep_r1_r2_suffixes"            : True,  # Add --keep-r1-r2 flag for _R1/_R2 suffixes
            "splitcode_opts"                 : [f"--output={dummy_output_r1},{dummy_output_r2}"],  # Required for --keep to work
        }
        if string_to_log:
            log.info(string_to_log)
        return (splitcode_tool.execute(**splitcode_kwargs), pool_id)

def generate_splitcode_config_and_keep_files(
    inner_demux_barcode_map_df,
    pool_id,
    output_dir,
    max_hamming_dist=1,
    r1_trim_bp_right_of_barcode=None,
    r2_trim_bp_left_of_barcode=None
):
    """
    Generate splitcode config and keep files for a single pool.

    This function creates the TSV files needed by splitcode to demultiplex reads
    based on inline barcodes (barcode_3). It takes sample metadata from a DataFrame
    and generates both the config file (which specifies barcode sequences and matching
    parameters) and the keep file (which maps barcodes to output file prefixes).

    Args:
        inner_demux_barcode_map_df (pd.DataFrame): DataFrame with sample metadata.
            Required columns:
            - 'barcode_3': Inline barcode sequence
            - 'run': Sample library ID (e.g., "Sample1.lLib1.FLOWCELL.1")
            - 'muxed_run': Pool ID for grouping samples
            Index: sample names
        pool_id (str): Pool identifier to filter samples and name output files
        output_dir (str): Directory where demultiplexed FASTQs will be written
        max_hamming_dist (int): Maximum Hamming distance for fuzzy barcode matching.
            Default: 1 (allows 1 mismatch)
        r1_trim_bp_right_of_barcode (int, optional): Additional bases to trim from R1
            after removing the barcode. If None, only the barcode is trimmed.
        r2_trim_bp_left_of_barcode (int, optional): Bases to trim from R2 left side
            (currently unused - R2 barcodes not implemented)

    Returns:
        tuple: (config_file_path, keep_file_path, sample_library_ids)
            - config_file_path (str): Path to generated splitcode config TSV
            - keep_file_path (str): Path to generated splitcode keep TSV
            - sample_library_ids (list): List of sample_library_id values for this pool

    Config file format (TSV with header):
        tag        - Barcode sequence (e.g., "AAAAAAAA")
        id         - Barcode identifier (e.g., "Sample1_R1") - MUST match keep file
        locations  - FILE_NUMBER:START_BP:END_BP (e.g., "0:0:8" for R1 positions 0-8)
        distance   - Max hamming distance for fuzzy matching (0=exact, 1=1 mismatch, etc.)
        left       - Trim from left: "1" removes barcode, "1:N" removes barcode + N more bp
        right      - Trim from right (typically "0" for R1 barcodes)

    Keep file format (TSV, NO header):
        Column 1: barcode_id - MUST match "id" from config file
        Column 2: output_prefix - Path prefix for output FASTQs (without _R1/_R2.fastq)

    Example:
        >>> df = pd.DataFrame({
        ...     'barcode_3': ['AAAAAAAA', 'CCCCCCCC'],
        ...     'run': ['Sample1.lLib1', 'Sample2.lLib1'],
        ...     'muxed_run': ['Pool1', 'Pool1']
        ... }, index=['Sample1', 'Sample2'])
        >>> config, keep, ids = generate_splitcode_config_and_keep_files(
        ...     df, 'Pool1', '/output', max_hamming_dist=1
        ... )
        >>> # config file contains:
        >>> # tag  id  locations  distance  left  right
        >>> # AAAAAAAA  Sample1.lLib1_R1  0:0:8  1  1  0
        >>> # CCCCCCCC  Sample2.lLib1_R1  0:0:8  1  1  0
        >>> # keep file contains:
        >>> # Sample1.lLib1_R1  /output/Sample1.lLib1
        >>> # Sample2.lLib1_R1  /output/Sample2.lLib1

    See test/unit/test_illumina_splitcode.py for validation of splitcode's behavior
    with these file formats.
    """
    # Filter DataFrame to only samples in this pool
    pool_samples_df = inner_demux_barcode_map_df[
        inner_demux_barcode_map_df['muxed_run'] == pool_id
    ]

    if len(pool_samples_df) == 0:
        raise ValueError(f"No samples found for pool_id '{pool_id}'")

    # ======== Create splitcode config file ========
    # For more information on config format and parameters see:
    #   https://splitcode.readthedocs.io/en/latest/reference_guide.html#table-options
    #   https://splitcode.readthedocs.io/en/latest/user_guide_tags.html#locations
    #   https://splitcode.readthedocs.io/en/latest/user_guide_tags.html#left-and-right-trimming

    config_file = util.file.mkstempfname(f'splitcode_{pool_id}_config.txt')

    with open(config_file, "w") as config_fh:
        config_tsv_writer = csv.writer(config_fh, delimiter="\t", lineterminator='\n')

        # Write header columns
        config_header = [
            "tag",        # Barcode sequence
            "id",         # Barcode identifier (must match keep file)
            "locations",  # FILE:START:END format
            "distance",   # Max hamming distance
            "left",       # Trim from left
            "right"       # Trim from right
        ]
        config_tsv_writer.writerow(config_header)

        # Write config lines for each sample in this pool
        for sample_name, sample_row in pool_samples_df.iterrows():
            barcode_sequence = sample_row["barcode_3"]
            barcode_len = len(barcode_sequence)
            sample_library_id = sample_row["run"]

            # R1 barcode configuration
            # The "left" column controls trimming from the left side of R1
            # Format: "1" (trim barcode only) or "1:N" (trim barcode + N more bp)
            left_trim = "1" if r1_trim_bp_right_of_barcode is None else f"1:{r1_trim_bp_right_of_barcode}"

            # When using --keep-r1-r2 flag, splitcode automatically adds _R1/_R2 suffixes
            # to output filenames, so the ID should NOT include _R1 suffix
            config_line_r1 = [
                barcode_sequence,              # The barcode sequence to search for
                sample_library_id,             # ID (MUST match keep file column 1)
                f"0:0:{barcode_len}",          # locations: FILE:START:END (0=R1, 0-barcode_len)
                str(max_hamming_dist),         # Maximum hamming distance
                left_trim,                     # Trim from left (barcode + optional extra bp)
                "0"                            # Trim from right (not used for R1)
            ]
            config_tsv_writer.writerow(config_line_r1)

            # TODO: R2 barcode support (currently commented out in splitcode_demux)
            # If needed, would add config_line_r2 here using r2_trim_bp_left_of_barcode

    # ======== Create splitcode keep file ========
    # Keep file maps barcode IDs to output file prefixes
    # NO header row (unlike config file)

    keep_file = util.file.mkstempfname(f"splitcode_{pool_id}_keepfile.txt")
    sample_library_ids = []

    with open(keep_file, "w") as keep_fh:
        keep_tsv_writer = csv.writer(keep_fh, delimiter="\t", lineterminator='\n')

        for sample_name, sample_row in pool_samples_df.iterrows():
            sample_library_id = sample_row["run"]
            sample_library_ids.append(sample_library_id)

            # Output prefix: /output/dir/Sample1.lLib1
            # When using --keep-r1-r2, splitcode will create:
            #   /output/dir/Sample1.lLib1_R1.fastq, ...R2.fastq
            sample_output_prefix = f"{output_dir}/{sample_library_id}"

            keep_row = [
                sample_library_id,       # MUST match config file "id" column (no _R1 suffix)
                sample_output_prefix,    # Output prefix (without _R1/_R2.fastq)
            ]
            keep_tsv_writer.writerow(keep_row)

    return (config_file, keep_file, sample_library_ids)

def convert_splitcode_demux_metrics_to_picard_style(
    in_splitcode_csv_path,
    out_splitcode_tsv_metrics_path,
    demux_function                 = "viral-core.SplitcodeMetrics",
    catchall_name                  = "unmatched",
    combine_innerbarcode_unmatched = False,
    report_within_pools            = True
):
    """
    Convert a custom Splitcode demux CSV file into a Picard-style
    ExtractIlluminaBarcodes.BarcodeMetric TSV file.

    :param in_splitcode_csv_path: Path to the Splitcode-format CSV input
    :param out_splitcode_tsv_metrics_path: Path to the desired Picard-style metrics TSV output
    :param demux_function: The string to insert after 'ExtractIlluminaBarcodes$' in
                           the '## METRICS CLASS' comment line
    :param catchall_name: A placeholder name for collapsed 'N' rows (not strictly used here)
    :param combine_innerbarcode_unmatched: If True, we collapse "all-N" rows into one single row.
                                           If False, we leave them as-is (no collapsing).
    :param report_within_pools: If True, the ratio/percentage columns are computed for each unique
                                (barcode_1, barcode_2) group. If False, they are computed globally.
                                Also, if True, we check only the *last* barcode segment for 'N's.
                                If False, we check the entire combined barcode for 'N's.
    """

    # Required columns for the first set of metrics
    required_cols = [
        "sample",                    # -> BARCODE_NAME
        "barcode_1",                 # -> part of BARCODE
        "barcode_2",                 # -> part of BARCODE
        "num_reads_hdistance0",      # -> PERFECT_MATCHES, PF_PERFECT_MATCHES
        "num_reads_hdistance1",      # -> ONE_MISMATCH_MATCHES, PF_ONE_MISMATCH_MATCHES
        "num_reads_total"            # -> READS, PF_READS
    ]

    # Read the CSV
    df = pd.read_csv(in_splitcode_csv_path, dtype=str)

    # Check for presence of required columns
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' is missing in the input file.")
    # Also check that these columns are not all empty
    for col in required_cols:
        if df[col].dropna().empty:
            raise ValueError(
                f"Required column '{col}' is present but entirely empty (no usable data)."
            )

    # Convert numeric columns from string to numeric
    # These columns are read as strings but need to be numeric for calculations
    numeric_cols = ["num_reads_hdistance0", "num_reads_hdistance1", "num_reads_total"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    # If 'inline_barcode' is present, treat it as a potential "third" barcode
    barcode_3_col = "inline_barcode" if "inline_barcode" in df.columns else None

    # Combine up to 3 barcodes into 'BARCODE' (with '-') and 'BARCODE_WITHOUT_DELIMITER' (concatenated)
    def combine_barcodes(row):
        bcs = []
        for bc_col in ["barcode_1", "barcode_2"]:
            val = row.get(bc_col, "")
            if pd.notnull(val) and val != "":
                bcs.append(val)
        if barcode_3_col and pd.notnull(row.get(barcode_3_col, "")) and row[barcode_3_col] != "":
            bcs.append(row[barcode_3_col])
        barcode = "-".join(bcs)
        barcode_no_delim = "".join(bcs)
        return pd.Series([barcode, barcode_no_delim])

    df[["BARCODE", "BARCODE_WITHOUT_DELIMITER"]] = df.apply(combine_barcodes, axis=1)

    # Assign 'BARCODE_NAME' and 'LIBRARY_NAME'
    df["BARCODE_NAME"] = df["sample"]
    if "run" in df.columns and not df["run"].dropna().empty:
        df["LIBRARY_NAME"] = df["run"]
    elif "library_id" in df.columns and not df["library_id"].dropna().empty:
        df["LIBRARY_NAME"] = df["library_id"]
    else:
        df["LIBRARY_NAME"] = df["BARCODE_NAME"]

    # Fill the first set of columns (READS, PF_READS, etc.)
    df["READS"]                   = df["num_reads_total"]
    df["PF_READS"]                = df["num_reads_total"]
    df["PERFECT_MATCHES"]         = df["num_reads_hdistance0"]
    df["PF_PERFECT_MATCHES"]      = df["num_reads_hdistance0"]
    df["ONE_MISMATCH_MATCHES"]    = df["num_reads_hdistance1"]
    df["PF_ONE_MISMATCH_MATCHES"] = df["num_reads_hdistance1"]

    # Define the "all N" checkers
    def is_all_N(barcode_str: str) -> bool:
        """
        Return True if the entire (combined) barcode is all 'N' and non-empty.
        Ignores delimiter characters ('-').
        """
        if not isinstance(barcode_str, str):
            return False
        # Remove delimiters before checking
        barcode_no_delim = barcode_str.replace('-', '')
        return len(barcode_no_delim) > 0 and all(ch == 'N' for ch in barcode_no_delim)

    def last_barcode_is_all_N(barcode_str: str) -> bool:
        """
        Return True if the final piece (e.g. inline barcode) is non-empty and all 'N'.
        For example, 'ACTGATCG-NNNNNNNNN' => True for that last piece.
        """
        if not isinstance(barcode_str, str):
            return False
        parts = barcode_str.split("-")
        if not parts:
            return False
        return is_all_N(parts[-1])

    # Row Collapsing for "all-N" if combine_innerbarcode_unmatched == True
    if combine_innerbarcode_unmatched:
        # The definition of "all-N" changes depending on report_within_pools:
        #  - If report_within_pools=True => check ONLY the last barcode piece
        #  - Otherwise => check entire combined barcode
        check_all_n_func = last_barcode_is_all_N if report_within_pools else is_all_N

        # Identify any rows that are "all N" by that definition
        df["IS_ALL_N"] = df["BARCODE"].apply(check_all_n_func)

        # Separate them out
        df_all_n     = df[df["IS_ALL_N"]]
        df_not_all_n = df[~df["IS_ALL_N"]]

        if not df_all_n.empty:
            # Sum their metric columns
            sum_reads               = df_all_n["READS"].sum()
            sum_pf_reads            = df_all_n["PF_READS"].sum()
            sum_perfect_matches     = df_all_n["PERFECT_MATCHES"].sum()
            sum_pf_perfect_matches  = df_all_n["PF_PERFECT_MATCHES"].sum()
            sum_one_mismatch        = df_all_n["ONE_MISMATCH_MATCHES"].sum()
            sum_pf_one_mismatch     = df_all_n["PF_ONE_MISMATCH_MATCHES"].sum()

            collapsed = df_all_n.iloc[0].copy()
            collapsed["BARCODE"]                   = "N"
            collapsed["BARCODE_WITHOUT_DELIMITER"] = "N"
            collapsed["BARCODE_NAME"]              = ""
            collapsed["LIBRARY_NAME"]              = ""
            collapsed["READS"]                     = sum_reads
            collapsed["PF_READS"]                  = sum_pf_reads
            collapsed["PERFECT_MATCHES"]           = sum_perfect_matches
            collapsed["PF_PERFECT_MATCHES"]        = sum_pf_perfect_matches
            collapsed["ONE_MISMATCH_MATCHES"]      = sum_one_mismatch
            collapsed["PF_ONE_MISMATCH_MATCHES"]   = sum_pf_one_mismatch

            # Recombine the "not all-N" rows plus this one collapsed row
            df = pd.concat([df_not_all_n, pd.DataFrame([collapsed], dtype=str)], ignore_index=True)

            # After concat with dtype=str, re-convert numeric columns
            numeric_cols_to_reconvert = [
                "READS", "PF_READS", "PERFECT_MATCHES", "PF_PERFECT_MATCHES",
                "ONE_MISMATCH_MATCHES", "PF_ONE_MISMATCH_MATCHES"
            ]
            for col in numeric_cols_to_reconvert:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

        # Drop helper column
        df.drop(columns=["IS_ALL_N"], inplace=True, errors="ignore")

    #
    # Now define how we compute the ratio/percentage columns
    #
    # Because you specifically want "if report_within_pools=True, only the last barcode is considered
    # for 'N's in the PF_NORMALIZED_MATCHES exclusion," we will conditionally pick which function
    # to use in that step too.
    #
    exclude_for_mean_func = last_barcode_is_all_N if report_within_pools else is_all_N

    def compute_stats_per_group(group: pd.DataFrame) -> pd.DataFrame:
        # Compute sums, maxima, etc. *within* the group, fill ratio/percentage columns.
        sum_of_reads    = group["READS"].sum()
        max_of_reads    = group["READS"].max() if len(group) > 0 else 1
        sum_of_pf_reads = group["PF_READS"].sum()
        max_of_pf_reads = group["PF_READS"].max() if len(group) > 0 else 1

        # For PF_NORMALIZED_MATCHES, exclude rows that pass `exclude_for_mean_func`
        # (which is either last_barcode_is_all_N or is_all_N, depending on 'report_within_pools').
        group_for_mean = group[~group["BARCODE"].apply(exclude_for_mean_func)]
        mean_pf_reads = group_for_mean["PF_READS"].mean() if len(group_for_mean) > 0 else 1

        out = group.copy()
        if sum_of_reads == 0:
            out["PCT_MATCHES"] = 0
        else:
            out["PCT_MATCHES"] = out["PERFECT_MATCHES"] / sum_of_reads

        if max_of_reads == 0:
            out["RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT"] = 0
        else:
            out["RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT"] = out["READS"] / max_of_reads

        if sum_of_pf_reads == 0:
            out["PF_PCT_MATCHES"] = 0
        else:
            out["PF_PCT_MATCHES"] = out["PF_PERFECT_MATCHES"] / sum_of_pf_reads

        if max_of_pf_reads == 0:
            out["PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT"] = 0
        else:
            out["PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT"] = out["PF_READS"] / max_of_pf_reads

        if mean_pf_reads == 0:
            out["PF_NORMALIZED_MATCHES"] = 0
        else:
            out["PF_NORMALIZED_MATCHES"] = out["PF_PERFECT_MATCHES"] / mean_pf_reads

        return out

    # Perform the group-based or global stats
    if report_within_pools:
        # Group by (barcode_1, barcode_2)
        grouped = df.groupby(["barcode_1", "barcode_2"], group_keys=False)
        df = grouped[df.columns].apply(compute_stats_per_group)
    else:
        # Global
        df = compute_stats_per_group(df)

    # Add DEMUX_TYPE column to all rows
    df["DEMUX_TYPE"] = "inline_splitcode"

    # Prepare the output
    columns_out = [
        "BARCODE",
        "BARCODE_WITHOUT_DELIMITER",
        "BARCODE_NAME",
        "LIBRARY_NAME",
        "READS",
        "PF_READS",
        "PERFECT_MATCHES",
        "PF_PERFECT_MATCHES",
        "ONE_MISMATCH_MATCHES",
        "PF_ONE_MISMATCH_MATCHES",
        "PCT_MATCHES",
        "RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT",
        "PF_PCT_MATCHES",
        "PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT",
        "PF_NORMALIZED_MATCHES",
        "DEMUX_TYPE",
    ]

    with open(out_splitcode_tsv_metrics_path, "w") as f:
        # Write the Picard-style "## METRICS CLASS" header line
        f.write(f"## METRICS CLASS\t{demux_function}\n")
        # Write column header
        f.write("\t".join(columns_out) + "\n")

        # Write data rows
        for _, row in df.iterrows():
            row_values = [str(row.get(col, "")) for col in columns_out]
            f.write("\t".join(row_values) + "\n")
