#!/usr/bin/env python3
"""
Utilities for demultiplexing Illumina data.
"""

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import concurrent.futures
import csv
import gc
import glob
import itertools
import json
import logging
import os
import os.path
import re
import shutil
import sqlite3
import subprocess
import tempfile
import xml.etree.ElementTree
from collections import defaultdict

import arrow
import matplotlib.pyplot as plt
import numpy  as np
import pandas as pd

import tools.picard
import tools.samtools
from tools.splitcode import SplitCodeTool
import read_utils
import util.cmd
import util.file
import util.misc
from util.illumina_indices import IlluminaIndexReference, IlluminaBarcodeHelper

log = logging.getLogger(__name__)


# =========================
# ***  illumina_demux   ***
# =========================


def parser_illumina_demux(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "inDir",
        help="Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.",
    )
    parser.add_argument("lane", help="Lane number.", type=int)
    parser.add_argument("outDir", help="Output directory for BAM files.")

    parser.add_argument(
        "--outMetrics",
        help="Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file.",
        default=None,
    )
    parser.add_argument(
        "--commonBarcodes",
        help='''Write a TSV report of all barcode counts, in descending order. 
                                Only applicable for read structures containing "B"''',
        default=None,
    )
    parser.add_argument(
        "--max_barcodes",
        help="""Cap the commonBarcodes report length to this size (default: %(default)s)""",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "--sampleSheet",
        default=None,
        help="""Override SampleSheet. Input TSV or CSV file w/header and four named columns:
                                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2.
                                Default is to look for a SampleSheet.csv in the inDir.""",
    )
    parser.add_argument(
        "--runInfo",
        default=None,
        dest="runinfo",
        help="""Override RunInfo. Input xml file.
                                Default is to look for a RunInfo.xml file in the inDir.""",
    )
    parser.add_argument(
        "--flowcell",
        help="Override flowcell ID (default: read from RunInfo.xml).",
        default=None,
    )
    parser.add_argument(
        "--read_structure",
        help="Override read structure (default: read from RunInfo.xml).",
        default=None,
    )
    parser.add_argument(
        "--append_run_id",
        help="If specified, output filenames will include the flowcell ID and lane number.",
        action="store_true",
    )
    parser.add_argument(
        "--collapse_duplicated_barcodes",
        help="""If specified, reads from samples with duplicated barcodes or barcode pairs 
                will be collapsed into a single output for each distinct barcode (or distinct barcode pair). 
                Intended for protocols allowing additional demultiplexing downstream by other means 
                (ex. breaking out samples based on a third, inner barcode, added via swift-seq). 
                If not specified, an error will be raised if duplicated barcodes (or barcode pairs) 
                are present in the sample sheet. If a value is specified, it will be used as 
                the path to store output sample sheet with barcodes collapsed""",
        default=False,
        const=None,
        nargs='?',
        #action="store_true"
    )
    parser.add_argument(
        "--rev_comp_barcodes_before_demux",
        help="""Reverse complement barcodes before demultiplexing.
                If specified without setting a value, 
                    the "barcode_2" column will be reverse-complemented.
                If one or more values are specified, 
                    the columns with those names will be reverse-complemented.
                (and if not specified, barcodes will not be reverse-complemented)""",
        nargs='*',
        action=util.cmd.storeMultiArgsOrFallBackToConst,
        type=str,
        const=["barcode_2"],
    )
    parser.add_argument(
        "--out_meta_by_sample", help="Output json metadata by sample", default=None
    )
    parser.add_argument(
        "--out_meta_by_filename",
        help="Output json metadata by bam file basename",
        default=None,
    )
    parser.add_argument(
        "--out_runinfo", help="Output json metadata about the run", default=None
    )

    for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list:
        if opt not in ("read_structure", "num_processors"):
            parser.add_argument(
                "--" + opt,
                help="Picard ExtractIlluminaBarcodes "
                + opt.upper()
                + " (default: %(default)s)",
                default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt),
            )
    for opt in tools.picard.IlluminaBasecallsToSamTool.option_list:
        if opt == "adapters_to_check":
            parser.add_argument(
                "--" + opt,
                nargs="*",
                help="Picard IlluminaBasecallsToSam "
                + opt.upper()
                + " (default: %(default)s)",
                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt),
            )
        elif opt in ("read_structure", "num_processors"):
            pass
        else:
            parser.add_argument(
                "--" + opt,
                help="Picard IlluminaBasecallsToSam "
                + opt.upper()
                + " (default: %(default)s)",
                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt),
            )

    parser.add_argument(
        "--JVMmemory",
        help="JVM virtual memory size (default: %(default)s)",
        default=tools.picard.IlluminaBasecallsToSamTool.jvmMemDefault,
    )
    util.cmd.common_args(
        parser,
        (
            (
                "threads",
                tools.picard.IlluminaBasecallsToSamTool.defaults["num_processors"],
            ),
            ("loglevel", None),
            ("version", None),
            ("tmp_dir", None),
        ),
    )
    util.cmd.attach_main(parser, main_illumina_demux)
    return parser


def main_illumina_demux(args):
    """Read Illumina runs & produce BAM files, demultiplexing to one bam per sample, or
    for simplex runs, a single bam will be produced bearing the flowcell ID.
    Wraps together Picard's ExtractBarcodes (for multiplexed samples) and IlluminaBasecallsToSam
    while handling the various required input formats. Also can
    read Illumina BCL directories, tar.gz BCL directories.
    """

    # prepare
    illumina = IlluminaDirectory(args.inDir)
    illumina.load()

    if args.runinfo:
        runinfo = RunInfo(args.runinfo)
    else:
        runinfo = illumina.get_RunInfo()
    if args.flowcell:
        flowcell = args.flowcell
    else:
        flowcell = runinfo.get_flowcell()
    if args.run_start_date:
        run_date = args.run_start_date
    else:
        run_date = runinfo.get_rundate_american()
    if args.read_structure:
        read_structure = args.read_structure
    else:
        read_structure = runinfo.get_read_structure()
    if args.append_run_id:
        run_id = "{}.{}".format(flowcell, args.lane)
    else:
        run_id = None
    if args.sampleSheet:
        samples = SampleSheet(
            args.sampleSheet,
            only_lane=args.lane,
            append_run_id=run_id,
            allow_non_unique=True if args.collapse_duplicated_barcodes!=False else False,
            # barcode_columns_to_revcomp:
            #   For --rev_comp_barcodes_before_demux: 
            #    1) None if not passed
            #    2) 'barcode_2' if passed without value
            #    3) barcodes specified if values passed with --rev_comp_barcodes_before_demux
            barcode_columns_to_revcomp = args.rev_comp_barcodes_before_demux 
        )
    else:
        samples = illumina.get_SampleSheet(
            only_lane=args.lane,
            append_run_id=run_id,
            allow_non_unique=True if args.collapse_duplicated_barcodes!=False else False,
            # barcode_columns_to_revcomp:
            #   For --rev_comp_barcodes_before_demux: 
            #    1) None if not passed
            #    2) 'barcode_2' if passed without value
            #    3) barcodes specified if values passed with --rev_comp_barcodes_before_demux
            barcode_columns_to_revcomp = args.rev_comp_barcodes_before_demux
        )

    # === debugging ===
    # print sample rows BEFORE collapsing duplicates
    #samples.print_rows()
    # =================
    collapse_requested = (args.collapse_duplicated_barcodes is not False)
    if not samples.can_be_collapsed:
        if collapse_requested:
            log.warning(
                "'--collapse_duplicated_barcodes' specified, but no duplicated barcodes "
                "(or barcode pairs) were found in the sample sheet."
            )
    else:
        log.warning("Duplicated barcodes (or barcode pairs) found in the sample sheet.")
        if collapse_requested:
            log.info(
                "'--collapse_duplicated_barcodes' specified; collapsing each duplicated barcode into "
                "a distinct barcode (or barcode pair)..."
            )
            samples.collapse_sample_index_duplicates(output_tsv=args.collapse_duplicated_barcodes)

    # === debugging ===
    # print sample rows AFTER collapsing duplicates
    # samples.print_rows()
    #
    # write library params file earlier than it would normally be
    # so we can inspect it for debugging purposes
    #
    # basecalls_input = util.file.mkstempfname(
    #         ".txt", prefix=".".join(["library_params", flowcell, str(args.lane)])
    #     )
    # samples.make_params_file(args.outDir, basecalls_input)
    # log.error("library params file written")
    # exit(0)
    # =================

    link_locs = False
    # For HiSeq-4000/X runs, If Picard's CheckIlluminaDirectory is
    # called with LINK_LOCS=true, symlinks with absolute paths
    # may be created, pointing from tile-specific *.locs to the
    # single s.locs file in the Intensities directory.
    # These links may break if the run directory is moved.
    # We should begin by removing broken links, if present,
    # and call CheckIlluminaDirectory ourselves if a 's.locs'
    # file is present, but only if the directory check fails
    # since link_locs=true tries to create symlinks even if they
    # (or the files) already exist
    try:
        tools.picard.CheckIlluminaDirectoryTool().execute(
            illumina.get_BCLdir(), args.lane, read_structure, link_locs=link_locs
        )
    except subprocess.CalledProcessError as e:
        log.warning("CheckIlluminaDirectory failed for %s", illumina.get_BCLdir())
        if os.path.exists(os.path.join(illumina.get_intensities_dir(), "s.locs")):
            # recurse to remove broken links in directory
            log.info(
                "This run has an 's.locs' file; checking for and removing broken per-tile symlinks..."
            )
            broken_links = util.file.find_broken_symlinks(
                illumina.get_intensities_dir()
            )
            if len(broken_links):
                for lpath in broken_links:
                    log.info("Removing broken symlink: %s", lpath)
                    os.unlink(lpath)

            # call CheckIlluminaDirectory with LINK_LOCS=true
            link_locs = True

            log.info("Checking run directory with Picard...")
            tools.picard.CheckIlluminaDirectoryTool().execute(
                illumina.get_BCLdir(), args.lane, read_structure, link_locs=link_locs
            )
        else:
            log.error("CheckIlluminaDirectory failed for %s", illumina.get_BCLdir())

    multiplexed_samples = True if "B" in read_structure else False

    if multiplexed_samples:
        assert (
            samples is not None
        ), "This looks like a multiplexed run since 'B' is in the read_structure: a SampleSheet must be given."
    else:
        assert (
            samples == None
        ), "A SampleSheet may not be provided unless 'B' is present in the read_structure"
        if args.commonBarcodes:
            log.warning(
                "--commonBarcodes was set but 'B' is not present in the read_structure; emitting an empty file."
            )
            util.file.touch(args.commonBarcodes)

    # B in read structure indicates barcoded multiplexed samples
    if multiplexed_samples:
        # Picard ExtractIlluminaBarcodes
        extract_input = util.file.mkstempfname(
            ".txt", prefix=".".join(["barcodeData", flowcell, str(args.lane)])
        )
        barcodes_tmpdir = tempfile.mkdtemp(prefix="extracted_barcodes-")
        samples.make_barcodes_file(extract_input)
        out_metrics = (
            (args.outMetrics is None)
            and util.file.mkstempfname(".metrics.txt")
            or args.outMetrics
        )
        picardOpts = dict(
            (opt, getattr(args, opt))
            for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
            if hasattr(args, opt) and getattr(args, opt) != None
        )
        picardOpts["read_structure"] = read_structure
        tools.picard.ExtractIlluminaBarcodesTool().execute(
            illumina.get_BCLdir(),
            args.lane,
            extract_input,
            barcodes_tmpdir,
            out_metrics,
            picardOptions=picardOpts,
            JVMmemory=args.JVMmemory,
        )

        if args.commonBarcodes:
            barcode_lengths = re.findall(r"(\d+)B", read_structure)
            try:
                barcode1_len = int(barcode_lengths[0])
            except IndexError:
                barcode1_len = 0
            try:
                barcode2_len = int(barcode_lengths[1])
            except IndexError:
                barcode2_len = 0

            # this step can take > 2 hours on a large high-output flowcell
            # so kick it to the background while we demux
            # count_and_sort_barcodes(barcodes_tmpdir, args.commonBarcodes)
            executor = concurrent.futures.ProcessPoolExecutor()
            executor.submit(
                count_and_sort_barcodes,
                barcodes_tmpdir,
                args.commonBarcodes,
                barcode1_len,
                barcode2_len,
                truncateToLength=args.max_barcodes,
                threads=util.misc.sanitize_thread_count(args.threads),
            )

        # Picard IlluminaBasecallsToSam
        basecalls_input = util.file.mkstempfname(
            ".txt", prefix=".".join(["library_params", flowcell, str(args.lane)])
        )
        samples.make_params_file(args.outDir, basecalls_input)

    picardOpts = dict(
        (opt, getattr(args, opt))
        for opt in tools.picard.IlluminaBasecallsToSamTool.option_list
        if hasattr(args, opt) and getattr(args, opt) != None
    )
    picardOpts["run_start_date"] = run_date
    picardOpts["read_structure"] = read_structure
    if args.threads:
        picardOpts["num_processors"] = args.threads
    if not picardOpts.get("sequencing_center") and runinfo:
        picardOpts["sequencing_center"] = runinfo.get_machine()

    if picardOpts.get("sequencing_center"):
        picardOpts["sequencing_center"] = util.file.string_to_file_name(
            picardOpts["sequencing_center"]
        )

    if args.out_runinfo:
        with open(args.out_runinfo, "wt") as outf:
            json.dump(
                {
                    "sequencing_center": picardOpts["sequencing_center"],
                    "run_start_date": runinfo.get_rundate_iso(),
                    "read_structure": picardOpts["read_structure"],
                    "indexes": str(samples.indexes),
                    "run_id": runinfo.get_run_id(),
                    "lane": str(args.lane),
                    "flowcell": str(runinfo.get_flowcell()),
                    "lane_count": str(runinfo.get_lane_count()),
                    "surface_count": str(runinfo.get_surface_count()),
                    "swath_count": str(runinfo.get_swath_count()),
                    "tile_count": str(runinfo.get_tile_count()),
                    "total_tile_count": str(runinfo.tile_count()),
                    "sequencer_model": runinfo.get_machine_model(),
                },
                outf,
                indent=2,
            )

    # manually garbage collect to make sure we have as much RAM free as possible
    gc.collect()
    if multiplexed_samples:
        tools.picard.IlluminaBasecallsToSamTool().execute(
            illumina.get_BCLdir(),
            barcodes_tmpdir,
            flowcell,
            args.lane,
            basecalls_input,
            picardOptions=picardOpts,
            JVMmemory=args.JVMmemory,
        )

        # organize samplesheet metadata as json
        sample_meta = list(samples.get_rows())
        for row in sample_meta:
            row["lane"] = str(args.lane)
        if args.out_meta_by_sample:
            with open(args.out_meta_by_sample, "wt") as outf:
                json.dump(dict((r["sample"], r) for r in sample_meta), outf, indent=2)
        if args.out_meta_by_filename:
            with open(args.out_meta_by_filename, "wt") as outf:
                json.dump(dict((r["run"], r) for r in sample_meta), outf, indent=2)

    else:
        tools.picard.IlluminaBasecallsToSamTool().execute_single_sample(
            illumina.get_BCLdir(),
            os.path.join(args.outDir, flowcell + ".bam"),
            flowcell,
            args.lane,
            flowcell,
            picardOptions=picardOpts,
            JVMmemory=args.JVMmemory,
        )

    # clean up
    if multiplexed_samples:
        if args.commonBarcodes:
            log.info("waiting for commonBarcodes output to finish...")
            executor.shutdown(wait=True)
        os.unlink(extract_input)
        os.unlink(basecalls_input)
        shutil.rmtree(barcodes_tmpdir)
    illumina.close()
    log.info("illumina_demux complete")
    return 0


__commands__.append(("illumina_demux", parser_illumina_demux))

# ==========================
# ***  flowcell_metadata   ***
# ==========================


def parser_flowcell_metadata(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "outMetadataFile", help="path of file to which metadata will be written."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--inDir",
        dest="in_dir",
        help="Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.",
        default=None,
    )
    group.add_argument(
        "--runInfo", default=None, dest="run_info", help="""RunInfo.xml file."""
    )
    group.add_argument(
        "--flowcellID",
        dest="flowcell_id",
        help="flowcell ID (default: read from RunInfo.xml).",
        default=None,
    )
    util.cmd.common_args(
        parser,
        (
            (
                "threads",
                tools.picard.IlluminaBasecallsToSamTool.defaults["num_processors"],
            ),
            ("loglevel", None),
            ("version", None),
            ("tmp_dir", None),
        ),
    )
    util.cmd.attach_main(parser, main_flowcell_metadata)
    return parser


def main_flowcell_metadata(args):
    """Writes run metadata to a file"""

    if args.flowcell_id:
        machine_matches = RunInfo.get_machines_for_flowcell_id(args.flowcell_id)
        machine_match = None
        if len(machine_matches) > 1:
            raise LookupError(
                "Multiple sequencers found for flowcell ID: %s"
                % " ".join([m["machine"] for m in machine_matches])
            )
        if len(machine_matches) == 0:
            raise LookupError(
                "No sequencers found for flowcell ID '%s' " % args.flowcell_id
            )
        machine_match = machine_matches[0]
    if args.run_info:
        runinfo = RunInfo(args.run_info)
        machine_match = runinfo.infer_sequencer_model()
    if args.in_dir:
        illumina = IlluminaDirectory(args.in_dir)
        illumina.load()
        runinfo = illumina.get_RunInfo()
        machine_match = runinfo.infer_sequencer_model()

    with open(args.outMetadataFile, "w") as outf:
        for k, v in machine_match.items():
            if type(v) == str and len(v) > 0 or type(v) != str:
                outline = "{k}\t{v}\n".format(k=k, v=v)
                print(outline, end="")
                outf.writelines([outline])


__commands__.append(("flowcell_metadata", parser_flowcell_metadata))

# ==========================
# ***  lane_metrics   ***
# ==========================


def parser_lane_metrics(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "inDir",
        help="Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.",
    )
    parser.add_argument(
        "outPrefix",
        help="""Prefix path to the *.illumina_lane_metrics and *.illumina_phasing_metrics files.""",
    )
    parser.add_argument(
        "--read_structure",
        help="Override read structure (default: read from RunInfo.xml).",
        default=None,
    )
    parser.add_argument(
        "--JVMmemory",
        help="JVM virtual memory size (default: %(default)s)",
        default=tools.picard.ExtractIlluminaBarcodesTool.jvmMemDefault,
    )
    util.cmd.common_args(
        parser, (("loglevel", None), ("version", None), ("tmp_dir", None))
    )
    util.cmd.attach_main(parser, main_lane_metrics)
    return parser


def main_lane_metrics(args):
    """
    Write out lane metrics to a tsv file.
    """
    # prepare
    illumina = IlluminaDirectory(args.inDir)
    illumina.load()
    if args.read_structure:
        read_structure = args.read_structure
    else:
        read_structure = illumina.get_RunInfo().get_read_structure()

    # Picard CollectIlluminaLaneMetrics
    output_dir = os.path.dirname(os.path.realpath(args.outPrefix))
    output_prefix = os.path.basename(os.path.realpath(args.outPrefix))

    picardOpts = dict(
        (opt, getattr(args, opt))
        for opt in tools.picard.CollectIlluminaLaneMetricsTool.option_list
        if hasattr(args, opt) and getattr(args, opt) != None
    )
    picardOpts["read_structure"] = read_structure
    tools.picard.CollectIlluminaLaneMetricsTool().execute(
        illumina.path,
        output_dir,
        output_prefix,
        picardOptions=picardOpts,
        JVMmemory=args.JVMmemory,
    )

    illumina.close()
    return 0


__commands__.append(("lane_metrics", parser_lane_metrics))


# ==========================
# ***  common_barcodes   ***
# ==========================


def parser_common_barcodes(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "inDir",
        help="Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.",
    )
    parser.add_argument("lane", help="Lane number.", type=int)
    parser.add_argument(
        "outSummary",
        help="""Path to the summary file (.tsv format). It includes several columns: 
                                            (barcode1, likely_index_name1, barcode2, likely_index_name2, count), 
                                            where likely index names are either the exact match index name for the barcode 
                                            sequence, or those Hamming distance of 1 away.""",
    )

    parser.add_argument(
        "--truncateToLength",
        help="If specified, only this number of barcodes will be returned. Useful if you only want the top N barcodes.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--omitHeader",
        help="If specified, a header will not be added to the outSummary tsv file.",
        action="store_true",
    )
    parser.add_argument(
        "--includeNoise",
        help='If specified, barcodes with periods (".") will be included.',
        action="store_true",
    )
    parser.add_argument(
        "--outMetrics",
        help="Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file.",
        default=None,
    )
    parser.add_argument(
        "--sampleSheet",
        default=None,
        help="""Override SampleSheet. Input tab or CSV file w/header and four named columns:
                                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2.
                                Default is to look for a SampleSheet.csv in the inDir.""",
    )
    parser.add_argument(
        "--flowcell",
        help="Override flowcell ID (default: read from RunInfo.xml).",
        default=None,
    )
    parser.add_argument(
        "--read_structure",
        help="Override read structure (default: read from RunInfo.xml).",
        default=None,
    )

    for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list:
        if opt not in ("read_structure", "num_processors"):
            parser.add_argument(
                "--" + opt,
                help="Picard ExtractIlluminaBarcodes "
                + opt.upper()
                + " (default: %(default)s)",
                default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt),
            )

    parser.add_argument(
        "--JVMmemory",
        help="JVM virtual memory size (default: %(default)s)",
        default=tools.picard.ExtractIlluminaBarcodesTool.jvmMemDefault,
    )
    util.cmd.common_args(
        parser,
        (("threads", None), ("loglevel", None), ("version", None), ("tmp_dir", None)),
    )
    util.cmd.attach_main(parser, main_common_barcodes)
    return parser


def main_common_barcodes(args):
    """
    Extract Illumina barcodes for a run and write a TSV report
    of the barcode counts in descending order
    """

    # prepare
    illumina = IlluminaDirectory(args.inDir)
    illumina.load()
    if args.flowcell:
        flowcell = args.flowcell
    else:
        flowcell = illumina.get_RunInfo().get_flowcell()
    if args.read_structure:
        read_structure = args.read_structure
    else:
        read_structure = illumina.get_RunInfo().get_read_structure()
    if args.sampleSheet:
        samples = SampleSheet(args.sampleSheet, only_lane=args.lane)
    else:
        samples = illumina.get_SampleSheet(only_lane=args.lane)

    # Picard ExtractIlluminaBarcodes
    barcode_file = util.file.mkstempfname(
        ".txt", prefix=".".join(["barcodeData", flowcell, str(args.lane)])
    )
    barcodes_tmpdir = tempfile.mkdtemp(prefix="extracted_barcodes-")
    samples.make_barcodes_file(barcode_file)
    out_metrics = (
        (args.outMetrics is None)
        and util.file.mkstempfname(".metrics.txt")
        or args.outMetrics
    )
    picardOpts = dict(
        (opt, getattr(args, opt))
        for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
        if hasattr(args, opt) and getattr(args, opt) != None
    )
    picardOpts["read_structure"] = read_structure
    tools.picard.ExtractIlluminaBarcodesTool().execute(
        illumina.get_BCLdir(),
        args.lane,
        barcode_file,
        barcodes_tmpdir,
        out_metrics,
        picardOptions=picardOpts,
        JVMmemory=args.JVMmemory,
    )

    barcode_lengths = re.findall(r"(\d+)B", read_structure)
    try:
        barcode1_len = int(barcode_lengths[0])
    except IndexError:
        barcode1_len = 0
    try:
        barcode2_len = int(barcode_lengths[1])
    except IndexError:
        barcode2_len = 0

    count_and_sort_barcodes(
        barcodes_tmpdir,
        args.outSummary,
        barcode1_len,
        barcode2_len,
        args.truncateToLength,
        args.includeNoise,
        args.omitHeader,
        args.threads,
    )

    # clean up
    os.unlink(barcode_file)
    shutil.rmtree(barcodes_tmpdir)
    illumina.close()
    return 0


__commands__.append(("common_barcodes", parser_common_barcodes))


def count_and_sort_barcodes(
    barcodes_dir,
    outSummary,
    barcode1_len=8,
    barcode2_len=8,
    truncateToLength=None,
    includeNoise=False,
    omitHeader=False,
    threads=None,
):
    # collect the barcode file paths for all tiles
    tile_barcode_files = [
        os.path.join(barcodes_dir, filename) for filename in os.listdir(barcodes_dir)
    ]

    # count all of the barcodes present in the tile files
    log.info("reading barcodes in all tile files")

    with util.file.CountDB() as reduce_db:
        barcodefile_tempfile_tuples = [
            (tile_barcode_file, util.file.mkstempfname("sqlite_.db"))
            for tile_barcode_file in tile_barcode_files
        ]

        # scatter tile-specific barcode files among workers to store barcode counts in SQLite
        workers = util.misc.sanitize_thread_count(threads)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(
                    util.file.count_occurrences_in_tsv_sqlite_backed,
                    tf,
                    bf,
                    include_noise=includeNoise,
                )
                for bf, tf in barcodefile_tempfile_tuples
            ]
            for future in concurrent.futures.as_completed(futures):
                tmp_db, barcode_file = future.result()
                log.debug(
                    "done reading barcodes from %s; adding to total...", barcode_file
                )
                # gather and reduce counts from separate SQLite databases into one
                reduce_db.add_counts_from_other_db(tmp_db)
                os.unlink(tmp_db)

        illumina_reference = IlluminaIndexReference()

        log.info("Number of barcodes seen %s", reduce_db.get_num_IDS())

        # write the barcodes and their corresponding counts
        with open(outSummary, "w") as tsvfile:
            log.info("sorting counts...")
            log.info("writing output...")
            writer = csv.writer(tsvfile, delimiter="\t")
            # write the header unless the user has specified not to do so
            if not omitHeader:
                writer.writerow(
                    (
                        "Barcode1",
                        "Likely_Index_Names1",
                        "Barcode2",
                        "Likely_Index_Names2",
                        "Count",
                    )
                )

            for num_processed, row in enumerate(reduce_db.get_counts_descending()):

                if truncateToLength and num_processed > truncateToLength:
                    break

                barcode, count = row

                writer.writerow(
                    (
                        barcode[:barcode1_len],
                        ",".join(
                            [
                                x
                                for x in illumina_reference.guess_index(
                                    barcode[:barcode1_len], distance=1
                                )
                            ]
                            or ["Unknown"]
                        ),
                        barcode[barcode1_len : len(barcode)],
                        ",".join(
                            [
                                x
                                for x in illumina_reference.guess_index(
                                    barcode[barcode1_len : len(barcode)], distance=1
                                )
                            ]
                            or ["Unknown"]
                        ),
                        count,
                    )
                )

                if num_processed % 50000 == 0:
                    log.debug(
                        "written %s barcode summaries to output file", num_processed
                    )

    log.info("done")


# ======================================
# ***  guess_low-abundance_barcodes  ***
# ======================================


def parser_guess_barcodes(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "in_barcodes", help="The barcode counts file produced by common_barcodes."
    )
    parser.add_argument(
        "in_picard_metrics", help="The demultiplexing read metrics produced by Picard."
    )
    parser.add_argument(
        "out_summary_tsv",
        help="""Path to the summary file (.tsv format). It includes several columns: 
                                            (sample_name, expected_barcode_1, expected_barcode_2, 
                                            expected_barcode_1_name, expected_barcode_2_name, 
                                            expected_barcodes_read_count, guessed_barcode_1, 
                                            guessed_barcode_2, guessed_barcode_1_name, 
                                            guessed_barcode_2_name, guessed_barcodes_read_count, 
                                            match_type), 
                                            where the expected values are those used by Picard during demultiplexing
                                            and the guessed values are based on the barcodes seen among the data.""",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--readcount_threshold",
        default=None,
        type=int,
        help="""If specified, guess barcodes for samples with fewer than this many reads.""",
    )
    group.add_argument(
        "--sample_names",
        nargs="*",
        help="If specified, only guess barcodes for these sample names.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--outlier_threshold",
        help="threshold of how far from unbalanced a sample must be to be considered an outlier.",
        type=float,
        default=0.775,
    )
    parser.add_argument(
        "--expected_assigned_fraction",
        help="The fraction of reads expected to be assigned. An exception is raised if fewer than this fraction are assigned.",
        type=float,
        default=0.7,
    )
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument(
        "--number_of_negative_controls",
        help="If specified, the number of negative controls in the pool, for calculating expected number of reads in the rest of the pool.",
        type=int,
    )
    group2.add_argument(
        "--neg_control_prefixes",
        nargs="+",
        help="If specified, the sample name prefixes assumed for counting negative controls. Case-insensitive.",
        type=str,
        default=["neg", "water", "NTC", "H2O"],
    )
    parser.add_argument(
        "--rows_limit",
        default=1000,
        type=int,
        help="""The number of rows to use from the in_barcodes.""",
    )

    util.cmd.common_args(
        parser, (("loglevel", None), ("version", None), ("tmp_dir", None))
    )
    util.cmd.attach_main(parser, main_guess_barcodes, split_args=True)
    return parser


def main_guess_barcodes(
    in_barcodes,
    in_picard_metrics,
    out_summary_tsv,
    sample_names,
    outlier_threshold,
    expected_assigned_fraction,
    number_of_negative_controls,
    readcount_threshold,
    rows_limit,
    neg_control_prefixes,
):
    """
    Guess the barcode value for a sample name,
        based on the following:
         - a list is made of novel barcode pairs seen in the data, but not in the picard metrics
         - for the sample in question, get the most abundant novel barcode pair where one of the
           barcodes seen in the data matches one of the barcodes in the picard metrics (partial match)
         - if there are no partial matches, get the most abundant novel barcode pair

        Limitations:
         - If multiple samples share a barcode with multiple novel barcodes, disentangling them
           is difficult or impossible

    The names of samples to guess are selected:
      - explicitly by name, passed via argument, OR
      - explicitly by read count threshold, OR
      - automatically (if names or count threshold are omitted)
        based on basic outlier detection of deviation from an assumed-balanced pool with
        some number of negative controls
    """

    bh = util.illumina_indices.IlluminaBarcodeHelper(
        in_barcodes, in_picard_metrics, rows_limit
    )
    guessed_barcodes = bh.find_uncertain_barcodes(
        sample_names=sample_names,
        outlier_threshold=outlier_threshold,
        expected_assigned_fraction=expected_assigned_fraction,
        number_of_negative_controls=number_of_negative_controls,
        readcount_threshold=readcount_threshold,
        neg_control_prefixes=neg_control_prefixes,
    )
    bh.write_guessed_barcodes(out_summary_tsv, guessed_barcodes)


__commands__.append(("guess_barcodes", parser_guess_barcodes))


# ============================
# ***  IlluminaDirectory   ***
# ============================


class IlluminaDirectory(object):
    """A class that handles Illumina data directories"""

    def __init__(self, uri):
        self.uri = uri
        self.path = None
        self.tempDir = None
        self.runinfo = None
        self.samplesheet = None

    def __enter__(self):
        self.load()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return 0

    def load(self):
        if self.path is None:
            if "://" in self.uri:
                raise NotImplementedError("boto s3 download here uri -> tarball")
                # tarball = util.file.mkstempfname('.tar.gz')
                # # TODO: download here, uri -> tarball
                # self._extract_tarball(tarball)
                # os.unlink(tarball)
            else:
                if os.path.isdir(self.uri):
                    self.path = self.uri
                else:
                    self._extract_tarball(self.uri)
            self._fix_path()

    def _fix_path(self):
        assert self.path is not None
        # this is not the correct root-level directory
        # sometimes this points to one level up
        while True:
            if os.path.isdir(
                os.path.join(self.path, "Data", "Intensities", "BaseCalls")
            ):
                # found it! self.path is correct
                break
            else:
                subdirs = list(
                    os.path.join(self.path, x)
                    for x in os.listdir(self.path)
                    if os.path.isdir(os.path.join(self.path, x))
                )
                if len(subdirs) == 1:
                    # follow the rabbit hole
                    self.path = subdirs[0]
                else:
                    # don't know where to go now!
                    raise Exception(
                        "cannot find Data/Intensities/BaseCalls/ inside %s (%s)"
                        % (self.uri, self.path)
                    )

    def _extract_tarball(self, tarfile):
        self.tempDir = tempfile.mkdtemp(prefix="IlluminaDirectory-")
        self.path = self.tempDir
        util.file.extract_tarball(tarfile, self.tempDir)

    def close(self):
        if self.tempDir:
            shutil.rmtree(self.tempDir)
            self.tempDir = None

    def get_RunInfo(self):
        if self.runinfo is None:
            runinfo_file = os.path.join(self.path, "RunInfo.xml")
            util.file.check_paths(runinfo_file)
            self.runinfo = RunInfo(runinfo_file)
        return self.runinfo

    def get_SampleSheet(self, only_lane=None, append_run_id=None, **kwargs):
        if self.samplesheet is None:
            samplesheet_file = os.path.join(self.path, "SampleSheet.csv")
            util.file.check_paths(samplesheet_file)
            self.samplesheet = SampleSheet(
                samplesheet_file,
                only_lane=only_lane,
                append_run_id=append_run_id,
                **kwargs
            )
        return self.samplesheet

    def get_intensities_dir(self):
        return os.path.join(self.path, "Data", "Intensities")

    def get_BCLdir(self):
        return os.path.join(self.get_intensities_dir(), "BaseCalls")


# ==================
# ***  RunInfo   ***
# ==================


class RunInfo(object):
    """A class that reads the RunInfo.xml file emitted by Illumina
    MiSeq and HiSeq machines.
    """

    def __init__(self, xml_fname):
        self.fname = xml_fname
        self.root = xml.etree.ElementTree.parse(xml_fname).getroot()

    def get_fname(self):
        return self.fname

    def get_run_id(self):
        return self.root[0].attrib["Id"]

    def get_flowcell_raw(self):
        return self.root[0].find("Flowcell").text

    def get_flowcell(self):
        fc = self.get_flowcell_raw()
        # slice in the case where the ID has a prefix of zeros
        if re.match(r"^0+-", fc):
            if "-" in fc:
                # miseq often adds a bunch of leading zeros and a dash in front
                fc = "-".join(fc.split("-")[1:])
        # >=5 to avoid an exception here: https://github.com/broadinstitute/picard/blob/2.17.6/src/main/java/picard/illumina/IlluminaBasecallsToSam.java#L510
        # <= 15 to limit the bytes added to each bam record
        assert len(fc) >= 5, "The flowcell ID must be five or more characters in length"
        if len(fc) > 15:
            log.warning(
                "The provided flowcell ID is longer than 15 characters. Is that correct?"
            )
        return fc

    def _get_rundate_obj(self):
        """
        Access the text of the <Date> node in the RunInfo.xml file
        and returns an arrow date object.
        """
        rundate = self.root[0].find("Date").text
        # possible formats found in RunInfo.xml:
        #   "170712" (YYMMDD)
        #   "20170712" (YYYYMMDD)
        #   "6/27/2018 4:59:20 PM" (M/D/YYYY h:mm:ss A)
        #   "2021-04-21T20:48:39Z" (YYYY-MM-DDTHH:mm:ssZ) [seen on NextSeq 2000]
        datestring_formats = [
            "YYMMDD",
            "YYYYMMDD",
            "M/D/YYYY h:mm:ss A",
            "YYYY-MM-DDTHH:mm:ssZ",
        ]
        for datestring_format in datestring_formats:
            try:
                date_parsed = arrow.get(rundate, datestring_format)
                return date_parsed
            except arrow.parser.ParserError:
                pass
        raise arrow.parser.ParserError(
            "The date string seen in RunInfo.xml ('%s') did not match known Illumina formats: %s"
            % (rundate, datestring_formats)
        )

    def get_rundate_american(self):
        return str(self._get_rundate_obj().format("MM/DD/YYYY"))

    def get_rundate_iso(self):
        return str(self._get_rundate_obj().format("YYYY-MM-DD"))

    def get_machine(self):
        return self.root[0].find("Instrument").text

    def get_read_structure(self):
        reads = []
        for x in self.root[0].find("Reads").findall("Read"):
            order = int(x.attrib["Number"])
            read = x.attrib["NumCycles"] + (
                x.attrib["IsIndexedRead"] == "Y" and "B" or "T"
            )
            reads.append((order, read))
        return "".join([r for _, r in sorted(reads)])

    def num_reads(self):
        return sum(
            1
            for x in self.root[0].find("Reads").findall("Read")
            if x.attrib["IsIndexedRead"] == "N"
        )

    def get_lane_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["LaneCount"])

    def get_surface_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["SurfaceCount"])

    def get_swath_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["SwathCount"])

    def get_tile_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["TileCount"])

    def get_section_count(self):
        layout = self.root[0].find("FlowcellLayout")
        # not ever flowcell type has sections but some do (ex. NextSeq 550 does)
        # return 1 in the event it's not listed in the RunInfo.xml file
        return int(layout.attrib.get("SectionPerLane", 1))

    def tile_count(self):
        lane_count = self.get_lane_count()
        surface_count = self.get_surface_count()
        swath_count = self.get_swath_count()
        tile_count = self.get_tile_count()
        section_count = self.get_section_count()

        total_tile_count = (
            lane_count * surface_count * swath_count * tile_count * section_count
        )
        return total_tile_count

    def machine_model_from_tile_count(self):
        """
        Return machine name and lane count based on tile count
        Machine names aim to conform to the NCBI SRA controlled
        vocabulary for Illumina sequencers available here:
          https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.common.xsd?view=co&content-type=text%2Fplain
        """
        tc = self.tile_count()

        machine = None
        if tc == 2:
            log.info("Detected %s tiles, interpreting as MiSeq nano run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 8:
            log.info("Detected %s tiles, interpreting as MiSeq micro run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 16:
            log.info("Detected %s tiles, interpreting as iSeq run.", tc)
            machine = {"machine": "Illumina iSeq 100", "lane_count": 1}
        elif tc == 28:
            log.info("Detected %s tiles, interpreting as MiSeq run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 38:
            log.info("Detected %s tiles, interpreting as MiSeq run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 128:
            log.info("Detected %s tiles, interpreting as HiSeq2k run.", tc)
            machine = {"machine": "Illumina HiSeq 2500", "lane_count": 2}
        elif tc == 132:
            # NextSeq P2 kit can be used on either NextSeq 1000 or 2000
            # so we cannot know which from the tile count alone
            log.info("Detected %s tiles, interpreting as NextSeq 1000/2000 P2 run.", tc)
            machine = {"machine": "NextSeq 1000/2000", "lane_count": 1}
        elif tc == 264:
            log.info("Detected %s tiles, interpreting as NextSeq 2000 P3 run.", tc)
            machine = {"machine": "NextSeq 2000", "lane_count": 2}
        elif tc == 288:
            # NextSeq 550 is a NextSeq 500 that can also read arrays.
            # Since we cannot tell them apart based on tile count, we call it the 550
            log.info(
                "Detected %s tiles, interpreting as NextSeq 550 (mid-output) run.", tc
            )
            machine = {"machine": "NextSeq 550", "lane_count": 4}
        elif tc == 624:
            log.info(
                "Detected %s tiles, interpreting as Illumina NovaSeq 6000 run.", tc
            )
            machine = {"machine": "Illumina NovaSeq 6000", "lane_count": 2}
        elif tc == 768:
            # HiSeq 2000 and 2500 have the same number of tiles
            # Defaulting to the newer HiSeq 2500
            log.info("Detected %s tiles, interpreting as HiSeq2500 run.", tc)
            machine = {"machine": "Illumina HiSeq 2500", "lane_count": 8}
        elif tc == 864:
            # NextSeq 550 is a NextSeq 500 that can also read arrays.
            # Since we cannot tell them apart based on tile count, we call it the 550
            log.info(
                "Detected %s tiles, interpreting as NextSeq 550 (high-output) run.", tc
            )
            machine = {"machine": "NextSeq 550", "lane_count": 4}
        elif tc == 896:
            log.info("Detected %s tiles, interpreting as HiSeq4k run.", tc)
            machine = {"machine": "Illumina HiSeq 4000", "lane_count": 8}
        elif tc == 1408:
            log.info(
                "Detected %s tiles, interpreting as Illumina NovaSeq 6000 run.", tc
            )
            machine = {"machine": "Illumina NovaSeq 6000", "lane_count": 2}
        elif tc == 3744:
            log.info(
                "Detected %s tiles, interpreting as Illumina NovaSeq 6000 run.", tc
            )
            machine = {"machine": "Illumina NovaSeq 6000", "lane_count": 4}
        elif tc > 3744:
            log.info("Tile count: %s tiles (unknown instrument type).", tc)
        return machine

    def get_flowcell_chemistry(self):
        guessed_sequencer = self.infer_sequencer_model()
        return guessed_sequencer["chemistry"]

    def get_flowcell_lane_count(self):
        guessed_sequencer = self.infer_sequencer_model()
        try:
            return self.get_lane_count()
        except Exception as e:
            return guessed_sequencer.get("lane_count", None)

    def get_machine_model(self):
        guessed_sequencer = self.infer_sequencer_model()
        return guessed_sequencer["machine"]

    @classmethod
    def get_machines_for_flowcell_id(cls, fcid):
        sequencer_by_fcid = []
        for key in cls.flowcell_to_machine_model_and_chemistry:
            if re.search(key, fcid):
                sequencer_by_fcid.append(
                    cls.flowcell_to_machine_model_and_chemistry[key]
                )
        return sequencer_by_fcid

    def infer_sequencer_model(self):
        fcid = self.get_flowcell_raw()
        sequencer_by_tile_count = self.machine_model_from_tile_count()
        sequencers_by_fcid = self.get_machines_for_flowcell_id(fcid)

        if len(sequencers_by_fcid) > 1:
            raise LookupError("Multiple sequencers possible: %s", fcid)

        log.debug("self.tile_count(): %s", self.tile_count())

        # always return sequencer model based on flowcell ID, if we can
        if len(sequencers_by_fcid) > 0:
            if (
                sequencer_by_tile_count is not None
                and sequencers_by_fcid[0]["machine"]
                != sequencer_by_tile_count["machine"]
            ):
                log.warning(
                    "Sequencer type inferred from flowcell ID: %s does not match sequencer inferred from tile count: %s; is this a new machine type?"
                    % (
                        sequencers_by_fcid[0]["machine"],
                        sequencer_by_tile_count["machine"],
                    )
                )
            return sequencers_by_fcid[0]
        # otherwise return based on tile count if we can
        elif sequencer_by_tile_count is not None:
            log.warning(
                "Sequencer type unknown flowcell ID: %s, yet sequencer type was inferred for tile count: %s; is this a new flowcell ID pattern?"
                % (fcid, self.tile_count())
            )
            return sequencer_by_tile_count
        # otherwise we do not know
        else:
            log.warning(
                "Tile count: %s and flowcell ID: %s are both novel; is this a new machine type?"
                % (self.tile_count(), fcid)
            )
            return {"machine": "UNKNOWN", "lane_count": self.get_lane_count()}

    # Machine names aim to conform to the NCBI SRA controlled
    # vocabulary for Illumina sequencers available here:
    #   https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.common.xsd?view=co&content-type=text%2Fplain
    flowcell_to_machine_model_and_chemistry = {
        r"[A-Z,0-9]{5}AAXX": {
            "machine": "Illumina Genome Analyzer IIx",
            "chemistry": "All",
            "lane_count": 8,
            "note": "",
        },
        r"[A-Z,0-9]{5}ABXX": {
            "machine": "Illumina HiSeq 2000",
            "chemistry": "V2 Chemistry",
            "lane_count": 8,
            "note": "",
        },
        r"[A-Z,0-9]{5}ACXX": {
            "machine": "Illumina HiSeq 2000",
            "chemistry": "V3 Chemistry",
            "lane_count": 8,
            "note": "Also used on transient 2000E",
        },
        r"[A-Z,0-9]{5}(?:ANXX|AN\w\w)": {
            "machine": "Illumina HiSeq 2500",
            "chemistry": "V4 Chemistry",
            "lane_count": 8,
            "note": "High output",
        },
        r"[A-Z,0-9]{5}(?:ADXX|AD\w\w)": {
            "machine": "Illumina HiSeq 2500",
            "chemistry": "V1 Chemistry",
            "lane_count": 2,
            "note": "Rapid run",
        },
        r"[A-Z,0-9]{5}AMXX": {
            "machine": "Illumina HiSeq 2500",
            "chemistry": "V2 Chemistry (beta)",
            "lane_count": 2,
            "note": "Rapid run",
        },
        r"[A-Z,0-9]{5}(?:BCXX|BC\w\w)": {
            "machine": "Illumina HiSeq 2500",
            "chemistry": "V2 Chemistry",
            "lane_count": 2,
            "note": "Rapid run",
        },
        # NextSeq 550 is a NextSeq 500 that can also read arrays.
        # Since we cannot tell them apart based on tile count, we call it the 550
        r"[A-Z,0-9]{5}AFX\w": {
            "machine": "NextSeq 550",
            "chemistry": "Mid-Output NextSeq",
            "lane_count": 4,
            "note": "",
        },
        # NextSeq 550 is a NextSeq 500 that can also read arrays.
        # Since we cannot tell them apart based on tile count, we call it the 550
        r"[A-Z,0-9]{5}AGXX": {
            "machine": "NextSeq 550",
            "chemistry": "V1 Chemistry",
            "lane_count": 4,
            "note": "High-output",
        },
        # NextSeq 550 is a NextSeq 500 that can also read arrays.
        # Since we cannot tell them apart based on tile count, we call it the 550
        r"[A-Z,0-9]{5}(?:BGXX|BG\w\w)": {
            "machine": "NextSeq 550",
            "chemistry": "V2/V2.5 Chemistry",
            "lane_count": 4,
            "note": "High-output",
        },
        # r'[A-Z,0-9]{5}(?:AAAC|AAA\w)':{ # suffix not confirmed
        #     "machine":    "NextSeq 1000/2000",
        #     "chemistry":  "P2 Chemistry",
        #     "lane_count":  1,
        #     "note":       "Mid-output"
        # },
        # r'[A-Z,0-9]{5}(?:AAAC|AAA\w)':{ # suffix not confirmed
        #     "machine":    "NextSeq 2000",
        #     "chemistry":  "P3 Chemistry",
        #     "lane_count":  2,
        #     "note":       "High-output"
        # },
        r"[A-Z,0-9]{5}(?:BBXX|BB\w\w)": {
            "machine": "Illumina HiSeq 4000",
            "chemistry": "Illumina HiSeq 4000",
            "lane_count": 8,
            "note": "",
        },
        r"[A-Z,0-9]{5}(?:ALXX:AL\w\w)": {
            "machine": "HiSeq X Ten",
            "chemistry": "V1/V2.5 Chemistry",
            "lane_count": 8,
            "note": "",
        },
        r"[A-Z,0-9]{5}(?:CCXX:CC\w\w)": {
            "machine": "HiSeq X Ten",
            "chemistry": "V2/V2.5 Chemistry",
            "lane_count": 8,
            "note": "",
        },
        r"[A-Z,0-9]{5}DR\w\w": {
            "machine": "Illumina NovaSeq 6000",
            "chemistry": "V1 Chemistry",
            "lane_count": 2,
            "note": "S1/SP",
        },
        r"[A-Z,0-9]{5}DM\w\w": {
            "machine": "Illumina NovaSeq 6000",
            "chemistry": "V1 Chemistry",
            "lane_count": 2,
            "note": "S2",
        },
        r"[A-Z,0-9]{5}DS\w\w": {
            "machine": "Illumina NovaSeq 6000",
            "chemistry": "V1 Chemistry",
            "lane_count": 4,
            "note": "S4",
        },
        r"BNS417.*": {
            "machine": "Illumina iSeq 100",
            "chemistry": "V1",
            "lane_count": 1,
            "note": "AKA Firefly",
        },
        r"[0-9]{9}-\w{5}": {
            "machine": "Illumina MiSeq",
            "chemistry": "V1/V2/V3 Chemistry",
            "lane_count": 1,
            "note": "",
        },
    }


# ======================
# ***  SampleSheet   ***
# ======================


class SampleSheetError(Exception):
    def __init__(self, message, fname):
        super(SampleSheetError, self).__init__(
            "Failed to load {} ({})".format(fname, message)
        )


class SampleSheet(object):
    """A class that reads an Illumina SampleSheet.csv or alternative/simplified
    tab-delimited versions as well.
    """

    def __init__(
        self,
        infile,
        use_sample_name  = True,
        only_lane        = None,
        allow_non_unique = False,
        append_run_id    = None,
        collapse_duplicates = False,
        rev_comp_barcode_2        = False,
        barcode_columns_to_revcomp = None # list of (additional) column names to reverse complement
    ):
        self.fname            = infile
        self.use_sample_name  = use_sample_name
        if only_lane is not None:
            only_lane         = str(only_lane)
        self.only_lane        = only_lane
        self.allow_non_unique = allow_non_unique
        self.append_run_id    = append_run_id
        self.rows             = []
        self._rowsOriginal    = []
        # state-holding class attributes
        self.duplicate_rows_collapsed             = False
        self.barcodes_revcomped_relative_to_input = False # idempotent
        self.barcodes_revcomped_column_names      = set()

        self._detect_and_load_sheet(infile)

        barcode_columns_to_revcomp = barcode_columns_to_revcomp or []

        # see rev_comp_barcode_2 is set, or barcode_columns_to_revcomp is: present, iterable, and not empty
        if rev_comp_barcode_2 or (barcode_columns_to_revcomp and isinstance(barcode_columns_to_revcomp, (list, tuple, set, str))):
            columns_to_revcomp = ['barcode_2'] if rev_comp_barcode_2 else []
            columns_to_revcomp += barcode_columns_to_revcomp

            self.rev_comp_barcode_values( barcode_columns_to_revcomp=columns_to_revcomp, inplace=True)

        if self.can_be_collapsed:
            if not allow_non_unique:
                raise SampleSheetError("Duplicate indices found in sample sheet", infile)
            else:
                if collapse_duplicates:
                    self.collapse_sample_index_duplicates()
            

    def _detect_and_load_sheet(self, infile):
        if infile.endswith((".csv", ".csv.gz")):
            # one of a few possible CSV formats (watch out for line endings from other OSes)
            with util.file.open_or_gzopen(infile, "rU") as inf:
                header = None
                miseq_skip = False
                row_num = 0
                for line_no, line in enumerate(inf):
                    if line_no == 0:
                        # remove BOM from first line, if present
                        #   see: https://en.wikipedia.org/wiki/Byte_order_mark
                        line = line.replace("\ufeff", "")

                    # if this is a blank line, skip parsing and continue to the next line...
                    if len(line.rstrip("\r\n").strip()) == 0:
                        continue
                    csv.register_dialect(
                        "samplesheet", quoting=csv.QUOTE_MINIMAL, escapechar="\\"
                    )
                    row = next(
                        csv.reader([line.strip().rstrip("\n")], dialect="samplesheet")
                    )
                    row = [
                        item.strip() for item in row
                    ]  # remove leading/trailing whitespace from each item
                    if miseq_skip:
                        if line.startswith("[Data]"):
                            # start paying attention *after* this line
                            miseq_skip = False
                        # otherwise, skip all the miseq headers
                    elif line.startswith("["):
                        # miseq: ignore all lines until we see "[Data]"
                        miseq_skip = True
                    elif header is None:
                        header = row
                        if all(x in header for x in ["Sample_ID", "Index"]):
                            # this is a Broad Platform MiSeq-generated SampleSheet.csv
                            keymapper = {
                                "Sample_ID": "sample",
                                "Index": "barcode_1",
                                "Index2": "barcode_2",
                                "Sample_Name": "sample_name",
                            }
                            header = list(map(keymapper.get, header))
                        elif "Sample_ID" in header:
                            # this is a MiSeq-generated SampleSheet.csv
                            keymapper = {
                                "Sample_ID": "sample",
                                "index": "barcode_1",
                                "index2": "barcode_2",
                                "Sample_Name": "sample_name",
                            }
                            header = list(map(keymapper.get, header))
                        elif "SampleID" in header:
                            # this is a Broad Platform HiSeq-generated SampleSheet.csv
                            keymapper = {
                                "SampleID": "sample",
                                "Index": "barcode_1",
                                "Index2": "barcode_2",
                                "libraryName": "library_id_per_sample",
                                "FCID": "flowcell",
                                "Lane": "lane",
                            }
                            header = list(map(keymapper.get, header))
                        elif len(row) == 3:
                            # hopefully this is a Broad walk-up submission sheet (_web_iww_htdocs_seq...)
                            header = ["sample", "barcode_1", "barcode_2"]
                            if "sample" not in row[0].lower():
                                # this is an actual data row! (no header exists in this file)
                                row_num += 1
                                self.rows.append(
                                    {
                                        "sample": row[0],
                                        "barcode_1": row[1],
                                        "barcode_2": row[2],
                                        "row_num": str(row_num),
                                    }
                                )
                        else:
                            raise SampleSheetError("unrecognized filetype", infile)
                        for h in ("sample", "barcode_1"):
                            assert h in header
                    else:
                        # data rows
                        row_num += 1

                        # pad the row with null strings if it is shorter than the header list
                        # sometimes a MiSeq produces an out-of-spec CSV file that lacks trailing commas,
                        # removing null values that should be present to ensure a length match with the header
                        while len(row) < len(header):
                            row.append("")

                        assert len(header) == len(row)
                        row = dict((k, v) for k, v in zip(header, row) if k and v)
                        row["row_num"] = str(row_num)
                        if (
                            self.only_lane is not None
                            and row.get("lane")
                            and self.only_lane != row["lane"]
                        ):
                            continue
                        if ("sample" in row and row["sample"]) and (
                            "barcode_1" in row and row["barcode_1"]
                        ):
                            self.rows.append(row)
            # go back and re-shuffle miseq columns if use_sample_name applies
            if (
                self.use_sample_name
                and "sample_name" in header
                and all(row.get("sample_name") for row in self.rows)
            ):
                for row in self.rows:
                    row["library_id_per_sample"] = row["sample"]
                    row["sample"] = row["sample_name"]
            for row in self.rows:
                if "sample_name" in row:
                    del row["sample_name"]

        elif infile.endswith((".txt", ".txt.gz", ".tsv")):
            # ===================================================================================================
            # our custom tab file format: sample, barcode_1, barcode_2, library_id_per_sample, [barcode_3], [...]
            # Required tsv column headings are:
            #    'sample' (must correspond to a biological sample; must be unique within a samplesheet)
            #    'library_id_per_sample' (must be unique within a samplesheet)
            #    'barcode_1'
            #    'barcode_2' (if paired reads, omit if single-end)
            #
            # Additional columns are optional and may include:
            #    'barcode_3'
            #
            # Note:
            #   ('sample' x 'library_id_per_sample') must be unique within a samplesheet and correspond to 
            #   independent libraries from the same original sample. 
            #
            #   External software or pipelines using this format may require additional columns to be present, 
            #   such as:
            #     'library_strategy'
            #     'library_source'
            #     'library_selection'
            #     'design_description'
            #
            #     Controlled vocabulary following the strict ontology on the third tab of:
            #       https://www.ncbi.nlm.nih.gov/core/assets/sra/files/SRA_metadata_acc_example.xlsx
            # ===================================================================================================
            self.rows = []
            row_num = 0
            for row in util.file.read_tabfile_dict(infile):
                assert row.get("sample") and row.get("barcode_1")
                row_num += 1
                row["row_num"] = str(row_num)
                self.rows.append(row)
        else:
            raise SampleSheetError("unrecognized filetype", infile)

        if not self.rows:
            raise SampleSheetError("empty file", infile)

        # populate library IDs, run IDs (ie BAM filenames)
        for row in self.rows:
            row["library"] = row["sample"]
            if row.get("library_id_per_sample"):
                row["library"] += ".l" + row["library_id_per_sample"]
            row["run"] = row["library"]
        if len(set(row["run"] for row in self.rows)) != len(self.rows):
            if self.allow_non_unique:
                log.warning("non-unique library IDs in this lane")
                unique_count = {}
                for row in self.rows:
                    unique_count.setdefault(row["library"], 0)
                    unique_count[row["library"]] += 1
                    row["run"] += ".r" + str(unique_count[row["library"]])
            else:
                raise SampleSheetError("non-unique library IDs in this lane", infile)
        if self.append_run_id:
            for row in self.rows:
                row["run"] += "." + self.append_run_id

        # escape sample, run, and library IDs to be filename-compatible
        for row in self.rows:
            row["sample_original"] = row["sample"]
            row["sample"]  = util.file.string_to_file_name(row["sample"])
            row["library"] = util.file.string_to_file_name(row["library"])
            row["run"]     = util.file.string_to_file_name(row["run"])

        # are we single or double indexed?
        if all(row.get("barcode_2") for row in self.rows):
            self.indexes = 2
        elif any(row.get("barcode_2") for row in self.rows):
            raise SampleSheetError(
                "inconsistent single/double barcoding in sample sheet", infile
            )
        else:
            self.indexes = 1

    def collapse_sample_index_duplicates(self, output_tsv=None, overwrite_instance_data=True):
        """
        Read in a sample sheet TSV, detect repeated grouping by (barcode_1[, barcode_2]),
        and collapse those rows into a single row if a group has >=2 rows.
        Un-duplicated rows remain as-is.

        The grouping logic:
          * If barcode_2 is present, group by (barcode_1, barcode_2).
          * Otherwise, group only by (barcode_1).

        For each group of >=2:
          - If barcode_2 does not exist, sample => <barcode_1>
          - If barcode_2 exists:
              * If it is non-empty for that group, sample => <barcode_1>-<barcode_2>
              * If it is empty for that group, sample => <barcode_1>
          - Other columns => collapsed with an MD5-based logic if multiple distinct values exist.
        """
        assert len(self.rows) > 0, "No sample sheet rows to collapse"

        hash_if_longer_than=32

        # convert [dict()] -> dataframe
        #   pd.json_normalize() is used to allow for conversion of nested dicts
        #     see: https://stackoverflow.com/a/53831756
        #   for simple lists of dicts, pd.DataFrame() is sufficient
        df = pd.json_normalize(self.rows).fillna("")

        # Read data, preserve column order
        # df = pd.read_csv(input_tsv, sep="\t", dtype=str).fillna("")
        original_cols = df.columns.tolist()

        # Determine which columns to group by
        grouping_cols = ["barcode_1"]
        has_barcode2 = "barcode_2" in df.columns
        if has_barcode2:
            grouping_cols.append("barcode_2")

        # Identify rows that belong to a duplicated group
        duplicated_mask = df.duplicated(subset=grouping_cols, keep=False)

        # Rows not in duplicates: pass unchanged
        df_unique = df[~duplicated_mask].copy()

        # Rows that are duplicates
        df_duplicates = df[duplicated_mask].copy()

        # Group them
        grouped_dups = df_duplicates.groupby(grouping_cols, sort=False)

        collapsed_rows = []
        for group_keys, group_df in grouped_dups:
            """
            group_keys is either:
              - A single value if grouping_cols == ['barcode_1']
              - A tuple (b1, b2) if grouping_cols == ['barcode_1','barcode_2']
            """
            if not isinstance(group_keys, tuple):
                group_keys = (group_keys,)  # unify type

            b1 = group_keys[0]
            b2 = group_keys[1] if has_barcode2 else None

            if len(group_df) == 1:
                # Single row in this group -> not truly duplicated, pass unmodified
                collapsed_rows.append(group_df.iloc[0])
            else:
                # Actual duplicates => build one collapsed row
                row_dict = {}

                # If "barcode_1" is a column, set it
                if "barcode_1" in df.columns:
                    row_dict["barcode_1"] = b1

                # If "barcode_2" is a column, set it
                if has_barcode2:
                    row_dict["barcode_2"] = b2

                # Build the sample string
                # If no barcode_2 column => sample = <b1>
                # If b2 is not empty => sample = <b1>-<b2>
                # If b2 is empty => sample = <b1>
                if not has_barcode2:
                    row_dict["sample"] = f"{b1}"
                else:
                    if b2 and b2.strip() != "":
                        row_dict["sample"] = f"{b1}-{b2}"
                    else:
                        row_dict["sample"] = f"{b1}"

                # For other columns, apply the collapse function
                for col in original_cols:
                    if col in ("barcode_1", "barcode_2", "sample"):
                        continue
                    col_values = group_df[col].tolist()
                    row_dict[col] = util.misc.collapse_dup_strs_to_str_or_md5(
                        col_values,
                        suffix="_muxed",
                        hash_if_longer_than=hash_if_longer_than
                    )

                # set the new run value; this becomes the bam file basename
                row_dict["run"] = f"{row_dict['sample']}.l{row_dict['library_id_per_sample']}" 

                if self.append_run_id:
                    row_dict["run"] += "." + self.append_run_id

                """
                File naming schemes seen in the wild:

                NovaSeq XP:
                  <SampleName>_S<SampleNumber>_L00<LaneNumber>_R<ReadNumber>_001.fastq.gz
                  ex.
                    NA10831_S1_L001_R1_001.fastq.gz

                NovaSeq standard:
                  <SampleName>_S<SampleNumber>_R<ReadNumber>_001.fastq.gz
                  ex.
                    NA10831_S1_R1_001.fastq.gz

                Broad walkup NovaSeq:
                  <Sample_ID>_S<SampleNumber>_L00<LaneNumber>_R<ReadNumber>_001.fastq.gz
                    Where Sample_ID is constructed as follows (empirical):
                      <fcid>_<lane>_<Sample_Name>
                        The <Sample_Name> is the original sample name, prefixed by 
                        a numeric string corresponding to the 
                        "Investigator Name" and/or "Experiment Name"
                        in the [header] block of the Illumina CSV 
                        returned with flowcell data. (ex. '0420593812' in the examples below)
                  ex.
                    22J5GLLT4_6_0420593812_B13Pool1a_S1_L006_R1_001.fastq.gz
                    22J5GLLT4_6_0420593812_VGG_21760_S49_L006_R2_001.fastq.gz

                Custom formats we have used:
                  <sample_name>.l<library_id>.<run_id>.<lane>.bam
                  ex.
                    USA-MA-Broad_MGH-22906-2024.lNDM_B1.HLWLWDRX5.1.bam
                    USA-MA-Broad_MGH-21721-2024.lDN_B20_C1.HJC2FDRX5.1.bam
                """

                collapsed_rows.append(pd.Series(row_dict))

        # Create a DataFrame from collapsed rows
        df_collapsed = pd.DataFrame(collapsed_rows) if collapsed_rows else pd.DataFrame(columns=original_cols)

        # Combine the collapsed duplicates with the unchanged unique rows
        out_df = pd.concat([df_collapsed, df_unique], ignore_index=True)
        out_df = out_df[original_cols]  # restore original column ordering


        # save collapsed tsv based on location of input tsv
        #if output_tsv is None:
        #    base = os.path.splitext(os.path.realpath(os.path.basename(self.fname)))[0]
        #    output_tsv = f"{base}_pools-collapsed.tsv"

        rows_collapsed = out_df.to_dict(orient="records")

        if len(rows_collapsed) < len(self.rows):
            # if we collapsed rows, log the change
            log.info("%s: %i rows collapsed ----> %i rows", os.path.basename(self.fname), len(self.rows), len(rows_collapsed))
            self.duplicate_rows_collapsed = True
        else:
            log.info("%s: ZERO rows collapsed (i.e. no duplicate barcodes or barcode pairs were present)", os.path.basename(self.fname), len(self.rows), len(rows_collapsed))

        excluded_tsv_output_columns=('row_num','library','run','sample_original')
        filtered_cols = [col for col in original_cols if col not in excluded_tsv_output_columns]
        if output_tsv is not None:
            log.info("Saving collapsed sample sheet to: %s", os.path.realpath(output_tsv))
            out_df[filtered_cols].to_csv(output_tsv, sep="\t", index=False)

        if overwrite_instance_data:
            self._rowsOriginal = self.rows
            self.rows = rows_collapsed

        return rows_collapsed

    @property
    def can_be_collapsed(self) -> bool:
        """
        Return True if, within the specified TSV, at least one group of rows
        shares the same (barcode_1[, barcode_2]) -- i.e., if collapsing would occur.
        Otherwise, return False.

        This function does NOT perform collapsing; it only checks the possibility.
        """
        assert len(self.rows) > 0, "No sample sheet rows to collapse"

        # convert [dict()] -> dataframe
        #   pd.json_normalize() is used to allow for conversion of nested dicts
        df = pd.json_normalize(self.rows).fillna("")

        grouping_cols = ["barcode_1"]
        if "barcode_2" in df.columns:
            grouping_cols.append("barcode_2")

        # Check if any row is duplicated based on grouping_cols
        duplicated_mask = df.duplicated(subset=grouping_cols, keep=False)
        # Return True if there's at least one group of 2 or more
        return duplicated_mask.any()

    def inner_demux_mapper(self):
        """
        Build a DataFrame mapping each (barcode_1,[barcode_2]) group to all 'barcode_3' values,
        preserving the original 'sample', assigning (or reusing) 'Inline_Index_ID',
        and computing two columns for each row:
          1) 'run'         -> "<sample>.l<library_id>[.<append_run_id>]"
          2) 'muxed_run'   -> A "collapsed sample" style run string:
                              "<muxed_sample>.l<library_id>[.<append_run_id>]"
             where muxed_sample is "<barcode_1>-<barcode_2>" if barcode_2 is non-empty,
             or just "<barcode_1>" otherwise.

        The final DataFrame will have 'sample' as its index.

        Each dict from self.rows corresponds to one row of data with keys as column headers.
        self.append_run_id : str, optional
            If provided, this string is appended to both 'run' and 'muxed_run'
            as ".{append_run_id}".

        Returns
        -------
        out_df : pd.DataFrame
            Columns:
              - barcode_1
              - barcode_2 (only if present)
              - barcode_3
              - sample             (becomes the index)
              - Inline_Index_ID
              - run
              - muxed_run
        """

        # Load into a DataFrame (flatten if nested)
        df = pd.json_normalize(self.rows)

        # Detect presence of barcode_2
        has_barcode2 = ("barcode_2" in df.columns)

        columns = df.columns

        # Grouping columns: always barcode_1; add barcode_2 if it exists
        grouping_cols = ["barcode_1"]
        if has_barcode2:
            grouping_cols.append("barcode_2")

        def fill_inline_index_ids(group):
            """
            If Inline_Index_ID doesn't exist, create it by enumerating rows in the group.
            If Inline_Index_ID exists but is empty, fill those rows with enumerated IDs.
            """
            if "Inline_Index_ID" not in group.columns:
                group["Inline_Index_ID"] = [str(i + 1) for i in range(len(group))]
            else:
                mask = (group["Inline_Index_ID"].isna()) | (group["Inline_Index_ID"] == "")
                next_ids = (str(i + 1) for i in range(len(group)))
                group.loc[mask, "Inline_Index_ID"] = [next(next_ids) for _ in range(sum(mask))]
            return group

        # Fill missing Inline_Index_ID values on a per-group basis
        df = df.groupby(grouping_cols, group_keys=False, as_index=False)
        df = df[columns].apply(fill_inline_index_ids, include_groups=False)

        # Build "run" column = "<sample>.l<library_id_per_sample>[.<append_run_id>]"
        def build_run_string(row):
            sample_val = row.get("sample", "")
            lib_val    = row.get("library_id_per_sample", "")
            run_str = f"{sample_val}.l{lib_val}"
            if self.append_run_id:
                run_str += f".{self.append_run_id}"
            return run_str

        df["run"] = df.apply(build_run_string, axis=1)

        # Build "muxed_run" column = 
        #   if barcode_2 is present & non-empty => "<barcode_1>-<barcode_2>"
        #   else => "<barcode_1>"
        #   then ".l<library_id_per_sample>[.<append_run_id>]"
        def build_muxed_run_string(row):
            b1 = row.get("barcode_1", "")
            b2 = row.get("barcode_2", "") if has_barcode2 else None
            if b2 and b2.strip():
                muxed_sample = f"{b1}-{b2}"
            else:
                muxed_sample = b1

            lib_val = row.get("library_id_per_sample", "")
            muxed_run_str = f"{muxed_sample}.l{lib_val}"
            if self.append_run_id:
                muxed_run_str += f".{self.append_run_id}"
            return muxed_run_str

        df["muxed_run"] = df.apply(build_muxed_run_string, axis=1)

        # Construct final output columns
        out_cols = ["barcode_1"]
        if has_barcode2:
            out_cols.append("barcode_2")
        out_cols.extend(["barcode_3", "sample", "Inline_Index_ID", "run", "muxed_run"])

        # Keep only columns that exist in df
        existing_cols = [c for c in out_cols if c in df.columns]
        out_df = df[existing_cols].copy()

        # Set 'sample' as the index
        out_df.set_index("sample", inplace=True)

        return out_df

    def make_barcodes_file(self, outFile):
        """Create input file for Picard ExtractBarcodes"""
        if self.num_indexes() == 2:
            header = [
                "barcode_name",
                "library_name",
                "barcode_sequence_1",
                "barcode_sequence_2",
            ]
        else:
            header = ["barcode_name", "library_name", "barcode_sequence_1"]
        with open(outFile, "wt") as outf:
            outf.write("\t".join(header) + "\n")
            for row in self.rows:
                out = {
                    "barcode_sequence_1": row["barcode_1"],
                    "barcode_sequence_2": row.get("barcode_2", ""),
                    "barcode_name": row["sample"],
                    "library_name": row["library"],
                }
                outf.write("\t".join(out[h] for h in header) + "\n")

    def write_tsv(self, outFile, force=False):
        """Write sample sheet to a tab-delimited file using csv dictwriter"""
        
        # before writing to outFile, check if it already
        # exists and raise an error unless force=True
        if os.path.exists(outFile) and not force:
            raise FileExistsError(f"Output file {outFile} already exists. Use force=True to overwrite.")

        with open(outFile, "wt") as outf:
            writer = csv.DictWriter(outf, self.rows[0].keys(), delimiter="\t")
            writer.writeheader()
            for row in self.rows:
                writer.writerow(row)

    def rev_comp_barcode_values(self, barcode_columns_to_revcomp=None, inplace=True):
        """
        Reverse-complement all barcode values in the sample sheet
        for the specified column(s). If no column is specified,
        only the 'barcode_2' column is reverse-complemented.
        If inplace=True, modify the instance data; otherwise, return a new SampleSheet object.
        """
        barcode_columns_to_revcomp = barcode_columns_to_revcomp or ["barcode_2"]

        # barcodes_revcomped_relative_to_input barcodes_revcomped_column_names

        if type(barcode_columns_to_revcomp) is str:
            barcode_columns_to_revcomp = [barcode_columns_to_revcomp]

        if inplace:
            self._rowsOriginal = self.rows.copy()
            for row_idx, row in enumerate(self.rows):
                for column_name in barcode_columns_to_revcomp:
                    if column_name in row:
                        try:
                            row[column_name] = util.misc.reverse_complement(row[column_name])
                        # If the barcode is not a valid DNA sequence, ignore it
                        except Exception as e:
                            log.warning("Failed to reverse-complement barcode value on line %s of %s: '%s'", row_idx+1, self.fname, row[column_name] )
                            row[column_name] = row[column_name]
                            pass
                        self.barcodes_revcomped_relative_to_input = True
                        self.barcodes_revcomped_column_names.add(column_name)
            return self
        else:
            new_sheet_fp = util.file.mkstempfname(f'{os.path.basename(self.fname)}_rev-comped-{"-".join(barcode_columns_to_revcomp)}.txt')
            log.debug("Creating a new SampleSheet object with reverse-complemented barcodes: %s", new_sheet_fp)
            shutil.copyfile(self.fname, new_sheet_fp)

            new_ss = SampleSheet(
                                    new_sheet_fp,
                                    barcode_columns_to_revcomp = barcode_columns_to_revcomp,
                                    allow_non_unique          = self.allow_non_unique
                                )
            return new_ss

    def make_params_file(self, bamDir, outFile):
        """Create input file for Picard IlluminaBasecallsToXXX"""
        if self.num_indexes() == 2:
            header = [
                "OUTPUT",
                "SAMPLE_ALIAS",
                "LIBRARY_NAME",
                "BARCODE_1",
                "BARCODE_2",
            ]
        else:
            header = ["OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_1"]
        with open(outFile, "wt") as outf:
            outf.write("\t".join(header) + "\n")
            # add one catchall entry at the end called Unmatched
            rows = self.rows + [
                {
                    "barcode_1": "N",
                    "barcode_2": "N",
                    "sample": "Unmatched",
                    "library": "Unmatched",
                    "run": "Unmatched",
                }
            ]
            for row in rows:
                out = {
                    "BARCODE_1": row["barcode_1"],
                    "BARCODE_2": row.get("barcode_2", ""),
                    "SAMPLE_ALIAS": row["sample"],
                    "LIBRARY_NAME": row["library"],
                }
                out["OUTPUT"] = os.path.join(bamDir, row["run"] + ".bam")
                outf.write("\t".join(out[h] for h in header) + "\n")

    def get_fname(self):
        return self.fname

    def get_rows(self):
        return self.rows

    def print_rows(self, row_indices=None):
        """ Print rows, optionally only a subset of them 
            specified via row_indices """
        rows_selected = [self.rows[i] for i in row_indices] if row_indices else self.rows

        for r in rows_selected:
            longest_key_len = len(max(r.keys(), key=len))
            print(f"Sample: {r['sample']}")
            for k,v in r.items():
                if k!='sample':
                    print(f"\t{k:<{longest_key_len+1}}: {v}")
            print("")

    def num_indexes(self):
        """
            Return 1 or 2 depending on whether pools are single or double indexed.
            Note that this field refers only to to Illumina adapter indices;
            it does not include nested inner barcodes (e.g., SwiftSeq barcode_3).
        """
        return self.indexes

    @property
    def num_samples(self):
        return len(self.rows)

    def fetch_by_index(self, idx):
        idx = str(idx)
        for row in self.rows:
            if idx == row["row_num"]:
                return row
        return None


# =============================
# ***  miseq_fastq_to_bam   ***
# =============================


def miseq_fastq_to_bam(
    outBam,
    sampleSheet,
    inFastq1,
    inFastq2=None,
    runInfo=None,
    sequencing_center=None,
    JVMmemory=tools.picard.FastqToSamTool.jvmMemDefault,
):
    """Convert fastq read files to a single bam file. Fastq file names must conform
    to patterns emitted by Miseq machines. Sample metadata must be provided
    in a SampleSheet.csv that corresponds to the fastq filename. Specifically,
    the _S##_ index in the fastq file name will be used to find the corresponding
    row in the SampleSheet
    """

    # match miseq based on fastq filenames
    mo = re.match(r"^\S+_S(\d+)_L001_R(\d)_001.fastq(?:.gz|)$", inFastq1)
    assert mo, (
        "fastq filename %s does not match the patterns used by an Illumina Miseq machine"
        % inFastq1
    )
    assert (
        mo.group(2) == "1"
    ), "fastq1 must correspond to read 1, not read %s" % mo.group(2)
    sample_num = mo.group(1)
    if inFastq2:
        mo = re.match(r"^\S+_S(\d+)_L001_R(\d)_001.fastq(?:.gz|)$", inFastq2)
        assert mo, (
            "fastq filename %s does not match the patterns used by an Illumina Miseq machine"
            % inFastq2
        )
        assert (
            mo.group(2) == "2"
        ), "fastq2 must correspond to read 2, not read %s" % mo.group(2)
        assert (
            mo.group(1) == sample_num
        ), "fastq1 (%s) and fastq2 (%s) must have the same sample number" % (
            sample_num,
            mo.group(1),
        )

    # load metadata
    samples = SampleSheet(sampleSheet, allow_non_unique=True)
    sample_info = samples.fetch_by_index(sample_num)
    assert sample_info, "sample %s not found in %s" % (sample_num, sampleSheet)
    sampleName = sample_info["sample"]
    log.info("Using sample name: %s", sampleName)
    if sample_info.get("barcode_2"):
        barcode = "-".join((sample_info["barcode_1"], sample_info["barcode_2"]))
    else:
        barcode = sample_info["barcode_1"]
    picardOpts = {
        "LIBRARY_NAME": sample_info["library"],
        "PLATFORM": "illumina",
        "VERBOSITY": "WARNING",
        "QUIET": "TRUE",
    }
    if runInfo:
        runInfo = RunInfo(runInfo)
        flowcell = runInfo.get_flowcell()
        picardOpts["RUN_DATE"] = runInfo.get_rundate_iso()
        if inFastq2:
            assert (
                runInfo.num_reads() == 2
            ), "paired fastqs given for a single-end RunInfo.xml"
        else:
            assert (
                runInfo.num_reads() == 1
            ), "second fastq missing for a paired-end RunInfo.xml"
    else:
        flowcell = "A"
    if sequencing_center is None and runInfo:
        sequencing_center = runInfo.get_machine()
    if sequencing_center:
        picardOpts["SEQUENCING_CENTER"] = util.file.string_to_file_name(
            sequencing_center
        )
    picardOpts["PLATFORM_UNIT"] = ".".join((flowcell, "1", barcode))
    if len(flowcell) > 5:
        flowcell = flowcell[:5]
    picardOpts["READ_GROUP_NAME"] = flowcell

    # run Picard
    picard = tools.picard.FastqToSamTool()
    picard.execute(
        inFastq1,
        inFastq2,
        sampleName,
        outBam,
        picardOptions=picard.dict_to_picard_opts(picardOpts),
        JVMmemory=JVMmemory,
    )
    return 0


def parser_miseq_fastq_to_bam(parser=argparse.ArgumentParser()):
    parser.add_argument("outBam", help="Output BAM file.")
    parser.add_argument("sampleSheet", help="Input SampleSheet.csv file.")
    parser.add_argument(
        "inFastq1", help="Input fastq file; 1st end of paired-end reads if paired."
    )
    parser.add_argument(
        "--inFastq2",
        help="Input fastq file; 2nd end of paired-end reads.",
        default=None,
    )
    parser.add_argument("--runInfo", help="Input RunInfo.xml file.", default=None)
    parser.add_argument(
        "--sequencing_center",
        default=None,
        help="Name of your sequencing center (default is the sequencing machine ID from the RunInfo.xml)",
    )
    parser.add_argument(
        "--JVMmemory",
        default=tools.picard.FastqToSamTool.jvmMemDefault,
        help="JVM virtual memory size (default: %(default)s)",
    )
    util.cmd.common_args(
        parser, (("loglevel", None), ("version", None), ("tmp_dir", None))
    )
    util.cmd.attach_main(parser, miseq_fastq_to_bam, split_args=True)
    return parser


__commands__.append(("miseq_fastq_to_bam", parser_miseq_fastq_to_bam))


# ==============================
# ***  extract_fc_metadata   ***
# ==============================
def extract_fc_metadata(flowcell, outRunInfo, outSampleSheet):
    """Extract RunInfo.xml and SampleSheet.csv from the provided Illumina directory"""
    illumina = IlluminaDirectory(flowcell)
    illumina.load()
    shutil.copy(illumina.get_RunInfo().get_fname(), outRunInfo)
    shutil.copy(illumina.get_SampleSheet().get_fname(), outSampleSheet)
    return 0


def parser_extract_fc_metadata(parser=argparse.ArgumentParser()):
    parser.add_argument("flowcell", help="Illumina directory (possibly tarball)")
    parser.add_argument("outRunInfo", help="Output RunInfo.xml file.")
    parser.add_argument("outSampleSheet", help="Output SampleSheet.csv file.")
    util.cmd.common_args(
        parser, (("loglevel", None), ("version", None), ("tmp_dir", None))
    )
    util.cmd.attach_main(parser, extract_fc_metadata, split_args=True)
    return parser


__commands__.append(("extract_fc_metadata", parser_extract_fc_metadata))


# ==================
# ***  Swiftseq/cDNAtag-Seq demux   ***
# ==================


def read_json(filepath):
    with open(filepath, "r") as file:
        data = json.load(file)
    return data


def create_lut(sample_sheet, outDir, unmatched_name, pool_ids=None, append_run_id=None, check_sample_sheet_consistency=False):
    pool_ids = pool_ids or []

    # Create look-up table that shows number of reads for each sample and maps i7/i5/inline barcodes

    # Load sample_sheet
    barcodes_df = pd.read_csv(sample_sheet, sep="\t")
    
    df_out = f"{outDir}/bc2sample_lut.csv"

    i7_barcodes = list(set(barcodes_df["barcode_1"].values))
    barcodes = sorted(list(set(barcodes_df["barcode_3"].values)))

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

    def duplication_check(df, primary_cols, secondary_col, error_message_header=None, error_message=None):
        default_error_message_header = "Error: More than one '{column_to_check_for_duplicates}' value present for distinct combinations of the columns {affected_column_names}:"
        default_error_message        = "'{duplicated_values}' appears {duplicate_count} times, for {affected_values}."

        error_message_header = error_message_header or default_error_message_header
        error_message        = error_message        or default_error_message

        # Check which placeholders are actually present in the template
        err_header_has_dup_check_col         = '{column_to_check_for_duplicates}' in error_message_header
        err_header_has_affected_column_names = '{affected_column_names}' in error_message_header

        err_msg_has_affected_values       = '{affected_values}' in error_message
        err_msg_has_dup_val               = '{duplicated_values}' in error_message
        err_msg_has_count                 = '{duplicate_count}' in error_message
        

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
            print(problem_header)
            for problem in problems:
                print(f"\t{problem}")
            problem_found = True
        if problem_found:
            raise ValueError("Problem(s) found in sample sheet; see above for details.")

    pool_dfs      = []
    unmatched_dfs = []

    for pool in barcodes_df["muxed_pool"].unique():
        # Get and load splitcode stats report json
        splitcode_summary_file = glob.glob(f"{outDir}/{pool}_summary.json")[0]
        splitcode_summary = read_json(splitcode_summary_file)
        
        samplesheet_rows_for_pool_df = barcodes_df[barcodes_df["muxed_pool"] == pool]

        splitcode_summary_df = pd.DataFrame.from_records(splitcode_summary["tag_qc"])

        splitcode_summary_df['run'] = splitcode_summary_df['tag'].copy()
        splitcode_summary_df['run'] = splitcode_summary_df['run'].str.removesuffix('_R1')
        
        splitcode_summary_df_h0_df = splitcode_summary_df[splitcode_summary_df["distance"] == 0]
        splitcode_summary_df_h1_df = splitcode_summary_df[splitcode_summary_df["distance"] == 1].copy()

        splitcode_summary_df_h1_df = splitcode_summary_df_h1_df.rename(columns={"count": "count_h1"})

        samplesheet_rows_for_pool_hx_df = samplesheet_rows_for_pool_df.join(
                                            splitcode_summary_df_h0_df.set_index('run'), 
                                            on='run')

        samplesheet_rows_for_pool_hx_df = pd.merge(samplesheet_rows_for_pool_hx_df,
                                            splitcode_summary_df_h1_df[['run','count_h1']].rename(columns={'run':'run_h1'}),
                                             left_on='run',
                                             right_on='run_h1',
                                             how='left')
        samplesheet_rows_for_pool_hx_df = samplesheet_rows_for_pool_hx_df.drop(columns=['run_h1','distance','tag']) # dropping distance since we've added a col with different distance (as indicated by _h1 suffix)
        # fil NA values in 'count_h1' and cast to int
        samplesheet_rows_for_pool_hx_df["count_h1"] = samplesheet_rows_for_pool_hx_df["count_h1"].fillna(0).astype(int)

        pool_dfs.append(samplesheet_rows_for_pool_hx_df)

        unmatched_dict = {
            "sample": f"{unmatched_name}.{pool}",
            "library_id_per_sample": list(set(samplesheet_rows_for_pool_hx_df["library_id_per_sample"]))[0],
            "run": f"{unmatched_name}.{pool}",
            "muxed_pool": pool,
            "count": splitcode_summary["n_processed"] - splitcode_summary["n_assigned"],
            "count_h1": 0,
            "barcode_1": list(samplesheet_rows_for_pool_hx_df["barcode_1"])[0],
            "barcode_2": list(samplesheet_rows_for_pool_hx_df["barcode_2"])[0],
            "barcode_3": "N" * len(list(samplesheet_rows_for_pool_hx_df["barcode_3"])[0]),
        }
        unmatched_df = pd.DataFrame.from_dict([unmatched_dict])
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

    df_lut.to_csv(df_out, index=False)

    return df_out

def plot_read_counts(df_lut_path, outDir):
    df_lut = pd.read_csv(df_lut_path)
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

    axs[1].set_yscale("log")

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


def plot_sorted_curve(df_lut_path, outDir, unmatched_name):
    df_lut = pd.read_csv(df_lut_path)

    fig, axs = plt.subplots(figsize=(10, 10), nrows=4, sharex=True)
    fontsize = 14

    # Define colors
    unique_library_ids = df_lut["library_id"].nunique()
    tab20_colors = plt.cm.tab20.colors
    pool_colors = (tab20_colors * (unique_library_ids // 20 + 1))[:unique_library_ids]

    reads_per_bc = pd.DataFrame()
    for i, pool in enumerate(sorted(df_lut["library_id"].unique())):
        num_reads = df_lut[
            (df_lut["library_id"] == pool)
            & (df_lut["inline_barcode"] != unmatched_name)
        ]["num_reads_total"]
        num_reads = sorted(num_reads, reverse=True)
        for ax in axs[:2]:
            ax.scatter(
                np.arange(len(num_reads)),
                num_reads,
                label=f'{pool.split("_")[-1]}',
                color=pool_colors[i],
            )
            ax.plot(np.arange(len(num_reads)), num_reads, color=pool_colors[i])

        # Calculate total and fractions
        total_reads = sum(num_reads)
        fractions = np.array([read / total_reads for read in num_reads])

        for ax in axs[2:]:
            ax.scatter(np.arange(len(fractions)), fractions * 100, color=pool_colors[i])
            ax.plot(np.arange(len(fractions)), fractions * 100, color=pool_colors[i])
        show_bottom_values = 5
        reads_per_bc[f"{pool}_reads"] = num_reads
        reads_per_bc[f"{pool}_reads_%"] = fractions * 100

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

    # Save data frame that contains the number of reads and fraction for each barcode in a sorted table
    reads_per_bc.to_csv(f"{outDir}/reads_per_bc.csv", index=False)

    plt.tight_layout()

    fig.savefig(
        f"{outDir}/reads_per_pool_sorted_curve.pdf", bbox_inches="tight", dpi=300
    )
    fig.savefig(
        f"{outDir}/reads_per_pool_sorted_curve.png", bbox_inches="tight", dpi=300
    )


def generate_bams_from_fastqs(
    sample,
    sampleSheet,
    bc_idxs,
    unmatched_name,
    outDir,
    platform_name,
    flowcell,
    lane,
    run_date,
    readgroup_name,
    sequencing_center,
):
    # Load sample sheet
    barcodes_df = pd.read_csv(sampleSheet, sep="\t")

    for bc_idx in bc_idxs + [unmatched_name]:
        r1_fastq = f"{outDir}/{sample}_demuxed/Barcode_{bc_idx}_0.fastq"
        r2_fastq = f"{outDir}/{sample}_demuxed/Barcode_{bc_idx}_1.fastq"

        # Get library/pool name (TO-DO: More stable way to get library name)
        library_name = sample.split("_")[-3]

        # Get i7 and i5 barcode sequences
        i7s = barcodes_df[barcodes_df["library_id_per_sample"] == library_name][
            "barcode_1"
        ].values
        if len(set(i7s)) > 1:
            raise ValueError(
                f"More than one i7 barcode sequence found for library {library_name}."
            )
        i7_barcode = i7s[0]

        i5s = barcodes_df[barcodes_df["library_id_per_sample"] == library_name][
            "barcode_2"
        ].values
        if len(set(i5s)) > 1:
            raise ValueError(
                f"More than one i5 barcode sequence found for library {library_name}."
            )
        i5_barcode = i5s[0]

        # Get sample name and inline barcode sequence
        if bc_idx != unmatched_name:
            barcodes_df_temp = barcodes_df[
                (barcodes_df["library_id_per_sample"] == library_name)
                & (barcodes_df["Inline_Index_ID"].astype(str) == str(bc_idx))
            ]

            # Get sample name
            sample_names = barcodes_df_temp["sample"].values
            if len(set(sample_names)) > 1:
                raise ValueError(
                    f"More than one sample name found for library {library_name} and barcode {bc_idx}."
                )
            sample_name = sample_names[0]

            # Get inline barcode sequence
            inline_seqs = barcodes_df_temp["barcode_3"].values
            if len(set(inline_seqs)) > 1:
                raise ValueError(
                    f"More than one inline barcode sequence found for library {library_name} and barcode {bc_idx}."
                )
            inline_barcode = inline_seqs[0]

        # Handle files where inline barcode was not matched
        else:
            sample_name = unmatched_name + "_" + library_name
            inline_barcode = unmatched_name + "_" + library_name

        platform_unit = f"{flowcell}.{lane}.{i7_barcode}-{i5_barcode}-{inline_barcode}"
        bam_outfolder = f"{outDir}/demuxed_bams"
        
        


        #tools.samtools.getReadGroups(inBam)

        picardOptions = [
            f"LIBRARY_NAME={library_name}",
            f"PLATFORM={platform_name}",
            f"PLATFORM_UNIT={platform_unit}",
            f"SEQUENCING_CENTER={sequencing_center}",
            f"READ_GROUP_NAME={readgroup_name}",
            f"RUN_DATE={run_date}",
        ]

        # output bam file named in the form <sample_name>.l<library_id>.<run_id>.<lane>.bam
        out_bam = f"{bam_outfolder}/{sample_name}.l{library_name}.{flowcell}.{lane}.bam"

        tools.picard.FastqToSamTool().execute( fastq1,
                                                fastqs2,
                                                sampleName,
                                                out_bam,
                                                picardOptions=picardOptions,
                                                JVMmemory=jvm_memory)

        # generate a new header to apply to the bam file
        #
        #tools.samtools.dumpHeader(inBam, outHeader)
        # read_utils.fastq_to_bam(
        #     inFastq1,
        #     inFastq2,
        #     outBam,
        #     sampleName=None,
        #     header=newHeader,
        #     JVMmemory=tools.picard.FastqToSamTool.jvmMemDefault,
        #     picardOptions=None
        # )

        # Change name of generated bam file
        # generated_bam = glob.glob(
        #     f"{bam_outfolder}/*_fastq_to_ubam/call-FastqToUBAM/work/{sample_name}.bam"
        # )[0]
        # destination = f"{bam_outfolder}/{sample_name}.{flowcell}.{lane}.bam"
        # shutil.move(generated_bam, destination)


def splitcode_demux(
    inDir,
    lane,
    outDir,

    max_hamming_dist=1,
    unmatched_name=None,

    illumina_run_directory=None,
    sampleSheet=None,
    runinfo=None,
    
    read_structure=None,
    platform_name=None,
    flowcell=None,
    run_id=None,
    run_date=None,
    #readgroup_name=None,
    sequencing_center=None,
    sequencer_model=None,

    rev_comp_barcodes_before_demux=None,

    threads=None,
    jvm_memory=None
):
    """
    Args:
    inDir                (str) File path to folder containing gzipped bam or fastq files from demux using outer barcodes
    outDir               (str) Output directory for BAM files and other output files
    sampleSheet          (str) Overwrite file path to csv samplesheet. Input tab or CSV file with columns:
                               sample, library_id_per_sample, I7_Index_ID, barcode_1, I5_Index_ID, barcode_2, Inline_Index_ID, barcode_3.
                               Default is to look for a SampleSheet.csv in the inDir.
    unmatched_name       (str) ID for reads that don't match an inline barcode. Default: "unmatched"
    max_hamming_dist     (int) Max allowed Hamming distance for inline barcode matching. Default: 1
    threads              (int) Threads used by splitcode demultiplexing. Default: 32
    runinfo              (str) Override RunInfo. Input xml file. Default: Look for a RunInfo.xml file in the inDir
    read_structure       (str) Read structure. Default: read from RunInfo.xml
    platform_name        (str) Platform name (used to populate BAM header). Default: read from RunInfo.xml
    flowcell             (str) Flowcell ID (used to populate BAM header). Default: read from RunInfo.xml
    lane                 (str) Lane number (used to populate BAM header). Default: read from RunInfo.xml
    run_date             (str) Run date (used to populate BAM header). Default: read from RunInfo.xml
    readgroup_name       (str) Readgroup name (used to populate BAM header). Default: read from RunInfo.xml
    sequencing_center    (str) Sequencing center (used to populate BAM header). Default: read from RunInfo.xml
    sequencer_model      (str) Sequencer model (used to populate BAM header). Default: inferred from RunInfo.xml
    
    rev_comp_barcodes_before_demux (list(str)) List of barcode columns to reverse complement before demux.
    """

    threads = threads or util.misc.available_cpu_count()
    unmatched_name = unmatched_name or "Unmatched"

    illumina_dir=None
    if illumina_run_directory:
        illumina_dir = IlluminaDirectory(illumina_run_directory)
        illumina_dir.load()

    # ToDo: if illumina_run_directory is not provided, but we lack provided runinfo, 
    #       inspect the inDir to see if it looks like an illumina run directory
    #       and try to load it as such (then interpreting Analysis/, fastq/ and other sub-dirs)
    #       as potential locations for the fastq files (maybe recursive find them)

    if runinfo:
        runinfo = RunInfo(runinfo)
    else:
        assert illumina_dir, "No Illumina Run Directory provided, and runinfo not provided."
        runinfo = illumina_dir.get_RunInfo()
    # data from RunInfo.xml is used in the absense of user input
    # so assert that we have runinfo available
        
    # if RunInfo.xml is needed to obtain any values not provided, 
    # so unless all values are provided, assert that runinfo is available
    if not all([flowcell, run_date, read_structure, run_id, sequencing_center]):
        assert runinfo, "No RunInfo.xml provided, nor found in a specified Illumina Run Directory."

    flowcell          = flowcell          or runinfo.get_flowcell()
    run_date          = run_date          or runinfo.get_rundate_iso()
    read_structure    = read_structure    or runinfo.get_read_structure()
    run_id            = run_id            or f"{runinfo.get_flowcell()}.{lane}" # runinfo.get_run_id()
    sequencing_center = sequencing_center or runinfo.get_machine() #f"{runinfo.get_machine_model()}.{runinfo.get_machine()}"
    sequencing_center = util.file.string_to_file_name(sequencing_center)
    platform_model    = sequencer_model or runinfo.get_machine_model()
    platform_model    = util.file.string_to_file_name(platform_model)

    if not platform_name:
        log.warning("No platform name provided. Defaulting to 'ILLUMINA'.")
    platform_name     = platform_name or "ILLUMINA"

    if sampleSheet:
        samples = SampleSheet(
            sampleSheet,
            only_lane=lane,
            append_run_id=run_id,
            allow_non_unique=True,
            # barcode_columns_to_revcomp:
            #   For --rev_comp_barcodes_before_demux: 
            #    1) None if not passed
            #    2) 'barcode_2' if passed without value
            #    3) barcodes specified if values passed with --rev_comp_barcodes_before_demux
            barcode_columns_to_revcomp = rev_comp_barcodes_before_demux
        )
    else:
        samples = illumina_dir.get_SampleSheet(
            only_lane=lane,
            append_run_id=run_id,
            allow_non_unique=True,
            # barcode_columns_to_revcomp:
            #   For --rev_comp_barcodes_before_demux: 
            #    1) None if not passed
            #    2) 'barcode_2' if passed without value
            #    3) barcodes specified if values passed with --rev_comp_barcodes_before_demux
            barcode_columns_to_revcomp = rev_comp_barcodes_before_demux
        )

    # debugging / perhaps useful later
    lane_count = runinfo.get_lane_count()

    # Create outDir
    os.makedirs(outDir, exist_ok=True)
    out_demux_dir_path = f"{outDir}"
    os.makedirs(out_demux_dir_path, exist_ok=True)

    # Load samplesheet into dataframe
    barcodes_df = pd.json_normalize(samples.get_rows()).fillna("")

    inner_demux_barcode_map_df = None
    if samples.can_be_collapsed:
        #samples.collapse_sample_index_duplicates()
        #collapsed_barcodes_df = pd.json_normalize(samples.get_rows()).fillna("")
        inner_demux_barcode_map_df = samples.inner_demux_mapper()
    else:
        log.error("The outer (barcode_1,barcode_2) sequences in the sample sheet do not appear to be collapsible.")

    # TODO: guardrails around missing barcode_3 values or a mixture of rows with it present/absent

    # Load inline barcodes and inline barcode indices
    inline_barcodes = sorted(list(set(barcodes_df["barcode_3"].values)))

    # Instantiate the splitcode and samtools tools
    splitcode_tool = SplitCodeTool()
    samtools = tools.samtools.SamtoolsTool()

    demux_bams_to_sample_library_id_map         = defaultdict(list)
    pool_id_to_sample_library_id_map            = defaultdict(list)
    sample_to_demux_bam_map                     = {}
    pool_id_to_splitcode_keepfile               = {}
    pool_id_to_splitcode_configfile             = {}
    pool_id_to_pool_bam                         = {}
    sample_library_id_to_fastqs                 = {}
    pool_ids_successfully_demuxed_via_splitcode = []

    # Iterate over rows
    for sample_name, sample_row in inner_demux_barcode_map_df.iterrows():
        #b1 = sample_row["barcode_1"]
        #b3 = sample_row["barcode_3"]
        #inline_index_id = sample_row["Inline_Index_ID"]
        #run_str = sample_row["run"]
        muxed_pool_str = sample_row["muxed_run"]

        bam_to_glob_for = f"{sample_row['muxed_run']}*.bam"
        found_bam_files = glob.glob(f"{inDir}/{bam_to_glob_for}".replace("//","/"))
        found_bam_file = found_bam_files[0] if found_bam_files else None
        if found_bam_file:
            if muxed_pool_str not in pool_id_to_pool_bam:
                pool_id_to_pool_bam[muxed_pool_str] = found_bam_file

            demux_bams_to_sample_library_id_map[found_bam_file].append(sample_row["run"])
            sample_to_demux_bam_map[sample_name] = found_bam_file
            pool_id_to_sample_library_id_map[sample_row["muxed_run"]].append(sample_row["run"])
        else:
            raise FileNotFoundError(f"No bam file found: for {bam_to_glob_for}")

    demux_bams_to_sample_library_id_map = dict(demux_bams_to_sample_library_id_map)
    pool_id_to_sample_library_id_map = dict(pool_id_to_sample_library_id_map)

    for pool_id, sample_libraries in pool_id_to_sample_library_id_map.items():
        # ======== create splitcode config file ========
        splitcode_config = util.file.mkstempfname(f'splitcode_{pool_id}_config.txt')
        with open(splitcode_config, "w") as config_fh:
            config_tsv_writer = csv.writer(config_fh, delimiter="\t")

            # Write header columns
            config_header = [
                                "tag",
                                "id",
                                "location",
                                "distance",
                                "left",
                                "right"
                            ]

            config_tsv_writer.writerow(config_header)

            # Create splitcode config files
            for sample_name, sample_row in inner_demux_barcode_map_df.iterrows():
                barcode_sequence  = sample_row["barcode_3"]
                barcode_len       = len(barcode_sequence)
                inline_index_id   = sample_row["Inline_Index_ID"]
                sample_library_id = sample_row["run"]

                if sample_row["muxed_run"] == pool_id:
                    # Columns: Barcode, tag_name, location, max allowed Hamming distance, trim from left on/off, trim from right on/off
                    # R1 barcode (using the "left" (-> remove from the left) column to remove the inline barcodes from R1 reads)
                    config_line_r1 = [
                                        barcode_sequence, 
                                        f"{sample_library_id}_R1", 
                                        f"0:0:{barcode_len}", 
                                        str(max_hamming_dist), 
                                        "1", 
                                        "0"
                                     ]
                    config_tsv_writer.writerow(config_line_r1)

                    # === retain for now ===
                    # R2 barcode (using the "right" (-> remove from the right) column to remove the inline barcodes from R2 reads)
                    # config_line_r2 = [
                    #                     barcode_sequence, 
                    #                     f"{inline_index_id}_R2", 
                    #                     f"1:-{barcode_len}:0", 
                    #                     str(max_hamming_dist), 
                    #                     "0", 
                    #                     "1"
                    #                  ]
                    # #config_f.write("\t".join(config_line_r2)+"\n")
                    # config_tsv_writer.writerow(config_line_r2)
                    # =======================
        pool_id_to_splitcode_configfile[pool_id] = splitcode_config


        # ======== create splitcode keep file ==========
    
        # ToDo: consider... rather than create a single splitcode keep file for each pool
        #       provide the option to create 1..N splitcode keep files for each pool
        #       where N is the number of samples in each pool
        #       to allow for increased parallelization of the splitcode demuxing,
        #       one job per library (at the expense of more I/O)

        splitcode_sample_keep_file = util.file.mkstempfname(f"splitcode_{pool_id}_keepfile.txt")
        with open(splitcode_sample_keep_file, "w") as splitcode_sample_keep_fh:
            splitcode_sample_keepfile_tsv_writer = csv.writer(splitcode_sample_keep_fh, delimiter="\t")

            for sample_library_id in sample_libraries:
                # ToDo: write fastq intermediate files to tmp dir once conversion to bam is in place downstream
                sample_output_prefix = f"{out_demux_dir_path}/{sample_library_id}"
                splitcode_sample_keep_row = [
                                                f"{sample_library_id}_R1",
                                                sample_output_prefix,
                                                #f"Barcode_{inline_index_id}" + "\n"
                                                #f"{sample_name}" + "\n"
                                                #f"{sampled_demuxed_basename}"
                                            ]
                splitcode_sample_keepfile_tsv_writer.writerow(splitcode_sample_keep_row)

                assert (sample_library_id not in sample_library_id_to_fastqs), "%s detected twice when it should be present only once" % sample_library_id
                sample_library_id_to_fastqs[sample_library_id] = [f"{sample_output_prefix}_R1.fastq",f"{sample_output_prefix}_R2.fastq"]
                  
        pool_id_to_splitcode_keepfile[pool_id] = splitcode_sample_keep_file
        unmatched_sample_name = f"{unmatched_name}.{pool_id}"
        unmatched_output_prefix = f"{out_demux_dir_path}/{unmatched_sample_name}"
        sample_library_id_to_fastqs[unmatched_sample_name] = [f"{unmatched_output_prefix}_R1.fastq",f"{unmatched_output_prefix}_R2.fastq"]

            # ToDo: glob the input files for either fastq or bam patterns
            # if the input is a bam file, convert to fastq files
            # 
            #with tools.samtools.SamtoolsTool().bam2fq_tmp(inBam) as (fqin1, fqin2), \
            #     util.file.tmp_dir('_splitcode') as t_dir:
            #     pass

            # Unzip fastq files using pigz
            #fastq_files = " ".join(f"'{f}'" for f in glob.glob(fastqs))
            #subprocess.run(["pigz", "-d", *fastq_files], check=True)            



    # TO-DO: Enable parallel processing of pools/samples
    # scatter tile-specific barcode files among workers to store barcode counts in SQLite
    
    # workers = util.misc.sanitize_thread_count(threads)
    # with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
    #     futures = [
    #         executor.submit(
    #             util.file.count_occurrences_in_tsv_sqlite_backed,
    #             tf,
    #             bf,
    #             include_noise=includeNoise,
    #         )
    #         for bf, tf in barcodefile_tempfile_tuples
    #     ]
    #     for future in concurrent.futures.as_completed(futures):
    #         tmp_db, barcode_file = future.result()
    #         log.debug(
    #             "done reading barcodes from %s; adding to total...", barcode_file
    #         )
    #         # gather and reduce counts from separate SQLite databases into one
    #         reduce_db.add_counts_from_other_db(tmp_db)
    #         os.unlink(tmp_db)

    #tmp dir to store intermediate fastq files from splitcode
    #with util.file.tmp_dir('_splitcode') as t_dir:

    # ======== run splitcode ==========

    for pool_id in sorted(pool_id_to_sample_library_id_map.keys()):
        #break # TODO: remove this
        splitcode_config   = pool_id_to_splitcode_configfile[pool_id]
        splitcode_keepfile = pool_id_to_splitcode_keepfile[pool_id]

        pool_bam_file = pool_id_to_pool_bam[pool_id]

        with tools.samtools.SamtoolsTool().bam2fq_tmp(pool_bam_file) as (fqin1, fqin2):
            n_fastqs = 2

            # output args
            unmapped_r1   = f"{out_demux_dir_path}/{unmatched_name}.{pool_id}_R1.fastq"
            unmapped_r2   = f"{out_demux_dir_path}/{unmatched_name}.{pool_id}_R2.fastq"
            summary_stats = f"{out_demux_dir_path}/{pool_id}_summary.json"

            # Optional: Pass 'predemux_r2_trim_3prime_num_bp':8 for '--trim-3 0,8' 
            # to remove 8 bases from end of R2 read 
            # (note: trimming happens before demuxing) 
            splitcode_kwargs={
                "n_fastqs": n_fastqs,
                "threads": threads,
                "config_file": splitcode_config,
                "keep_file": splitcode_keepfile,
                "unassigned_r1": unmapped_r1,
                "unassigned_r2": unmapped_r2,
                "summary_stats": summary_stats,
                "r1": fqin1,
                "r2": fqin2,
                "splitcode_opts": ["--no-output", "--no-outb"]
            }
            splitcode_tool.execute(**splitcode_kwargs)


            # maps sample_library_id to fastqs
            # ex.
            #   'VGG_19883.lPool_3.HHJYWDRX5.1' -> ["VGG_19883.lPool_3.HHJYWDRX5.1_R1.fastq","VGG_19883.lPool_3.HHJYWDRX5.1_R2.fastq"]
            #print(f"sample_library_id_to_fastqs {sample_library_id_to_fastqs}")
            pool_ids_successfully_demuxed_via_splitcode.append(pool_id)

            # ToDo: convert the splitcode output file fastqs to bams
            #       and reheader them with the appropriate readgroups and other attributes

    
    # Create look-up table that shows number of reads for each sample and maps i7/i5/inline barcodes
    df_out = create_lut(sampleSheet, outDir, unmatched_name, pool_ids_successfully_demuxed_via_splitcode, append_run_id=run_id)

    picard_style_splitcode_metrics_path = os.path.splitext(df_out)[0] + "_picard-style.txt"
    convert_splitcode_demux_metrics_to_picard_style(df_out, picard_style_splitcode_metrics_path)

    # Plot number and fraction of reads per inline barcode per pool
    plot_read_counts(df_out, outDir)

    # Plot a sorted curve per pool and save csv file with sorted read numbers and fractions for QC
    plot_sorted_curve(df_out, outDir, unmatched_name)

    
    df = pd.read_csv(df_out)
    df = df.rename(columns={"inline_barcode": "barcode_3"})

    # change index of df to 'run' column
    df = df.set_index('run')

    # Generate bam files
    for sample,(fq1,fq2) in sample_library_id_to_fastqs.items():
        # print the row of df with index matching sample
        if sample not in df.index:
            log.warning("Sample %s not found in the samplesheet data.", sample)
            continue

        samplesheet_row_for_sample = df.loc[sample]
        samplesheet_row_idx = df.index.get_loc(sample)
        samplesheet_row_idx = samplesheet_row_idx if type(samplesheet_row_idx)==int else None

        barcodes_for_sample_joined = "-".join([
                                        samplesheet_row_for_sample["barcode_1"], 
                                        samplesheet_row_for_sample["barcode_2"], 
                                        samplesheet_row_for_sample["barcode_3"]])
        platform_unit = ".".join([
                                    flowcell,
                                    str(lane), 
                                    barcodes_for_sample_joined
                                ])

        sample = samplesheet_row_for_sample.name #long-form name in the form {sample}.l{library_name}.{flowcell}.{lane}
        sample_name = samplesheet_row_for_sample["sample"]
        library_name = samplesheet_row_for_sample["library_id"]
        readgroup_name = ".".join([
                                    flowcell[:4],
                                    str(lane),
                                    str(samplesheet_row_idx+1) if samplesheet_row_idx is not None else util.misc.md5_digest(barcodes_for_sample_joined,last_n_chr=4)
                                    ])
        picardOptions = [
            f"LIBRARY_NAME={sample}",
            f"PLATFORM={platform_name}",
            f"PLATFORM_UNIT={platform_unit}",
            f"PLATFORM_MODEL={platform_model}",
            f"SEQUENCING_CENTER={sequencing_center}",
            f"RUN_DATE={run_date}",
            f"READ_GROUP_NAME={readgroup_name}",
        ]

        # output bam file named in the form <sample_name>.l<library_id>.<run_id>.<lane>.bam
        #out_bam = f"{out_demux_dir_path}/{sample}.l{library_name}.{flowcell}.{lane}.bam"
        out_bam = f"{out_demux_dir_path}/{sample}.bam"

        tools.picard.FastqToSamTool().execute(  fq1,
                                                fq2,
                                                sample_name,
                                                out_bam,
                                                picardOptions=picardOptions,
                                                JVMmemory=jvm_memory)

        # generate a new header to apply to the bam file
        #
        #tools.samtools.dumpHeader(inBam, outHeader)
        # read_utils.fastq_to_bam(
        #     inFastq1,
        #     inFastq2,
        #     outBam,
        #     sampleName=None,
        #     header=newHeader,
        #     JVMmemory=tools.picard.FastqToSamTool.jvmMemDefault,
        #     picardOptions=None
        # )

        # Change name of generated bam file
        # generated_bam = glob.glob(
        #     f"{bam_outfolder}/*_fastq_to_ubam/call-FastqToUBAM/work/{sample_name}.bam"
        # )[0]
        # destination = f"{bam_outfolder}/{sample_name}.{flowcell}.{lane}.bam"
        # shutil.move(generated_bam, destination)


        
        # generate_bams_from_fastqs(
        #     sample,
        #     sampleSheet,
        #     bc_idxs,
        #     unmatched_name,
        #     outDir,
        #     platform_name,
        #     flowcell,
        #     lane,
        #     run_date,
        #     readgroup_name,
        #     sequencing_center,
        # )


def main_splitcode_demux(args):
    """
    Main function to call splitcode_demux from the command line.
    """
    splitcode_demux(
        inDir=args.inDir,
        lane=args.lane,
        outDir=args.outDir,
        sampleSheet=args.sampleSheet,
        illumina_run_directory=args.illumina_run_directory,
        unmatched_name=args.unmatched_name,
        max_hamming_dist=args.max_hamming_dist,
        threads=args.threads,
        runinfo=args.runinfo,
        platform_name=args.platform_name,
        flowcell=args.flowcell,
        run_date=args.run_date,
        #readgroup_name=args.readgroup_name,
        sequencing_center=args.sequencing_center,
        rev_comp_barcodes_before_demux=args.rev_comp_barcodes_before_demux,
        jvm_memory=args.jvm_memory
    )


def parser_splitcode_demux(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="Demultiplex FASTQ files using inline barcodes with splitcode."
        )

    # Positional arguments
    parser.add_argument(
        "inDir", help="File path to folder containing gzipped FASTQ files."
    )
    parser.add_argument(
        "lane",
        default=None,
        help="Lane number (used to populate BAM header). Default: read from RunInfo.xml.",
    )
    parser.add_argument(
        "outDir", help="Output directory for BAM files and other output files."
    )

    # Optional arguments
    parser.add_argument(
        "--sampleSheet",
        default=None,
        help="""Override SampleSheet. Input tab or CSV file with columns:
                sample, library_id_per_sample, I7_Index_ID, barcode_1, I5_Index_ID, barcode_2, Inline_Index_ID, barcode_3..
                Default is to look for a SampleSheet.csv in the inDir.""",
    )
    parser.add_argument(
        "--illuminaRunDirectory",
        dest="illumina_run_directory",
        default=None,
        help="""Path to Illumina Run Directory. Default is to look for a RunInfo.xml file in the inDir.""",
    )

    # Optional arguments specific to splitcode_demux
    parser.add_argument(
        "--unmatched_name",
        default="unmatched",
        help="ID for reads that don't match an inline barcode. Default: 'unmatched'.",
    )
    parser.add_argument(
        "--max_hamming_dist",
        type=int,
        default=1,
        help="Max allowed Hamming distance for inline barcode matching. Default: 1.",
    )

    # Additional BAM header fields
    parser.add_argument(
        "--runInfo",
        default=None,
        dest="runinfo",
        help="""Override RunInfo. Input xml file.
                Default is to look for a RunInfo.xml file in the inDir.""",
    )
    parser.add_argument(
        "--platform_name",
        default=None,
        help="Platform name (used to populate BAM header). Default: read from RunInfo.xml.",
    )
    parser.add_argument(
        "--flowcell",
        default=None,
        help="Flowcell ID (used to populate BAM header). Default: read from RunInfo.xml.",
    )
    parser.add_argument(
        "--run_date",
        default=None,
        help="Run date (used to populate BAM header). Default: read from RunInfo.xml.",
    )
    # parser.add_argument(
    #     "--readgroup_name",
    #     default=None,
    #     help="Readgroup name (used to populate BAM header). Default: read from RunInfo.xml.",
    # )
    parser.add_argument(
        "--sequencing_center",
        default=None,
        help="Sequencing center (used to populate BAM header). Default: read from RunInfo.xml.",
    )
    parser.add_argument(
        "--rev_comp_barcodes_before_demux",
        help="""Reverse complement barcodes before demultiplexing.
                If specified without setting a value, 
                    the "barcode_2" column will be reverse-complemented.
                If one or more values are specified, 
                    the columns with those names will be reverse-complemented.
                (and if not specified, barcodes will not be reverse-complemented)""",
        nargs='*',
        action=util.cmd.storeMultiArgsOrFallBackToConst,
        type=str,
        const=["barcode_2"],
    )
    parser.add_argument(
        "--JVMmemory",
        dest="jvm_memory",
        help="JVM virtual memory size (default: %(default)s)",
        default=tools.picard.FastqToSamTool.jvmMemDefault,
    )

    # Attach common arguments. Adjust the default for threads if needed.
    util.cmd.common_args(
        parser,
        (("threads", 32), ("loglevel", None), ("version", None), ("tmp_dir", None)),
    )

    # Attach main function so that the parser can execute the main when run from CLI.
    util.cmd.attach_main(parser, main_splitcode_demux)
    return parser


__commands__.append(("splitcode_demux", parser_splitcode_demux))

def convert_splitcode_demux_metrics_to_picard_style(
    in_splitcode_csv_path,
    out_splitcode_tsv_metrics_path,
    demux_function="viral-core.SplitcodeMetrics",
    catchall_name="unmatched",
    combine_innerbarcode_unmatched=False,
    report_within_pools=True
):
    """
    Convert a custom Splitcode demux CSV file into a Picard-style
    ExtractIlluminaBarcodes.BarcodeMetric TSV file.

    :param in_splitcode_csv_path: Path to the Splitcode-format CSV input
    :param out_splitcode_tsv_metrics_path: Path to the desired Picard-style metrics TSV output
    :param demux_function: The string to insert after 'ExtractIlluminaBarcodes$' in
                           the '## METRICS CLASS' comment line
    :param catchall_name: A placeholder name for collapsed 'N' rows (not strictly used here)
    :param combine_innerbarcode_unmatched: If True, we collapse “all-N” rows into one single row.
                                           If False, we leave them as-is (no collapsing).
    :param report_within_pools: If True, the ratio/percentage columns are computed for each unique
                                (barcode_1, barcode_2) group. If False, they are computed globally.
                                Also, if True, we check only the *last* barcode segment for 'N's.
                                If False, we check the entire combined barcode for 'N's.
    """

    #
    # Required columns for the first set of metrics
    #
    required_cols = [
        "sample",                     # -> BARCODE_NAME
        "barcode_1",                 # -> part of BARCODE
        "barcode_2",                 # -> part of BARCODE
        "num_reads_hdistance0",      # -> PERFECT_MATCHES, PF_PERFECT_MATCHES
        "num_reads_hdistance1",      # -> ONE_MISMATCH_MATCHES, PF_ONE_MISMATCH_MATCHES
        "num_reads_total"            # -> READS, PF_READS
    ]

    # Read the CSV
    df = pd.read_csv(in_splitcode_csv_path)

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

    # If 'inline_barcode' is present, treat it as a potential "third" barcode
    barcode_3_col = "inline_barcode" if "inline_barcode" in df.columns else None

    #
    # Combine up to 3 barcodes into 'BARCODE' (with '-') and 'BARCODE_WITHOUT_DELIMITER' (concatenated)
    #
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

    #
    # Assign 'BARCODE_NAME' and 'LIBRARY_NAME'
    #
    df["BARCODE_NAME"] = df["sample"]
    if "run" in df.columns and not df["run"].dropna().empty:
        df["LIBRARY_NAME"] = df["run"]
    elif "library_id" in df.columns and not df["library_id"].dropna().empty:
        df["LIBRARY_NAME"] = df["library_id"]
    else:
        df["LIBRARY_NAME"] = df["BARCODE_NAME"]

    #
    # Fill the first set of columns (READS, PF_READS, etc.)
    #
    df["READS"] = df["num_reads_total"]
    df["PF_READS"] = df["num_reads_total"]
    df["PERFECT_MATCHES"] = df["num_reads_hdistance0"]
    df["PF_PERFECT_MATCHES"] = df["num_reads_hdistance0"]
    df["ONE_MISMATCH_MATCHES"] = df["num_reads_hdistance1"]
    df["PF_ONE_MISMATCH_MATCHES"] = df["num_reads_hdistance1"]

    #
    # Define the "all N" checkers
    #
    def is_all_N(barcode_str: str) -> bool:
        """Return True if the entire (combined) barcode is all 'N' and non-empty."""
        if not isinstance(barcode_str, str):
            return False
        return len(barcode_str) > 0 and all(ch == 'N' for ch in barcode_str)

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

    #
    # Row Collapsing for "all-N" if combine_innerbarcode_unmatched == True
    #
    if combine_innerbarcode_unmatched:
        # The definition of "all-N" changes depending on report_within_pools:
        #  - If report_within_pools=True => check ONLY the last barcode piece
        #  - Otherwise => check entire combined barcode
        check_all_n_func = last_barcode_is_all_N if report_within_pools else is_all_N

        # Identify any rows that are "all N" by that definition
        df["IS_ALL_N"] = df["BARCODE"].apply(check_all_n_func)

        # Separate them out
        df_all_n = df[df["IS_ALL_N"]]
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
            collapsed["BARCODE"] = "N"
            collapsed["BARCODE_WITHOUT_DELIMITER"] = "N"
            collapsed["BARCODE_NAME"] = ""
            collapsed["LIBRARY_NAME"] = ""
            collapsed["READS"] = sum_reads
            collapsed["PF_READS"] = sum_pf_reads
            collapsed["PERFECT_MATCHES"] = sum_perfect_matches
            collapsed["PF_PERFECT_MATCHES"] = sum_pf_perfect_matches
            collapsed["ONE_MISMATCH_MATCHES"] = sum_one_mismatch
            collapsed["PF_ONE_MISMATCH_MATCHES"] = sum_pf_one_mismatch

            # Recombine the "not all-N" rows plus this one collapsed row
            df = pd.concat([df_not_all_n, pd.DataFrame([collapsed])], ignore_index=True)

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
        """Compute sums, maxima, etc. *within* the group, fill ratio/percentage columns."""
        sum_of_reads = group["READS"].sum()
        max_of_reads = group["READS"].max() if len(group) > 0 else 1
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

    #
    # Perform the group-based or global stats
    #
    if report_within_pools:
        # Group by (barcode_1, barcode_2)
        grouped = df.groupby(["barcode_1", "barcode_2"], group_keys=False)
        df = grouped[df.columns].apply(compute_stats_per_group)
    else:
        # Global
        df = compute_stats_per_group(df)

    #
    # Prepare the output
    #
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

# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == "__main__":
    util.cmd.main_argparse(__commands__, __doc__)
