#!/usr/bin/env python3
"""
Utilities for demultiplexing Illumina data.
"""

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import os
import os.path
import re
import gc
import csv
import sqlite3, itertools
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree
from collections import defaultdict
import concurrent.futures

import arrow

import util.cmd
import util.file
import util.misc
import tools.picard
from util.illumina_indices import IlluminaIndexReference, IlluminaBarcodeHelper

log = logging.getLogger(__name__)

# =========================
# ***  illumina_demux   ***
# =========================


def parser_illumina_demux(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('outDir', help='Output directory for BAM files.')

    parser.add_argument('--outMetrics',
                        help='Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file.',
                        default=None)
    parser.add_argument('--commonBarcodes',
                        help='''Write a TSV report of all barcode counts, in descending order. 
                                Only applicable for read structures containing "B"''',
                        default=None)
    parser.add_argument('--max_barcodes',
                        help='''Cap the commonBarcodes report length to this size (default: %(default)s)''',
                        default=10000, type=int)
    parser.add_argument('--sampleSheet',
                        default=None,
                        help='''Override SampleSheet. Input tab or CSV file w/header and four named columns:
                                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2.
                                Default is to look for a SampleSheet.csv in the inDir.''')
    parser.add_argument('--runInfo',
                        default=None,
                        dest="runinfo",
                        help='''Override RunInfo. Input xml file.
                                Default is to look for a RunInfo.xml file in the inDir.''')
    parser.add_argument('--flowcell', help='Override flowcell ID (default: read from RunInfo.xml).', default=None)
    parser.add_argument('--read_structure',
                        help='Override read structure (default: read from RunInfo.xml).',
                        default=None)
    parser.add_argument('--append_run_id',
                        help='If specified, output filenames will include the flowcell ID and lane number.',
                        action='store_true')

    for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list:
        if opt not in ('read_structure', 'num_processors'):
            parser.add_argument('--' + opt,
                                help='Picard ExtractIlluminaBarcodes ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt))
    for opt in tools.picard.IlluminaBasecallsToSamTool.option_list:
        if opt == 'adapters_to_check':
            parser.add_argument('--' + opt,
                                nargs='*',
                                help='Picard IlluminaBasecallsToSam ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))
        elif opt in ('read_structure', 'num_processors'):
            pass
        else:
            parser.add_argument('--' + opt,
                                help='Picard IlluminaBasecallsToSam ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))

    parser.add_argument('--JVMmemory',
                        help='JVM virtual memory size (default: %(default)s)',
                        default=tools.picard.IlluminaBasecallsToSamTool.jvmMemDefault)
    util.cmd.common_args(parser, (('threads', tools.picard.IlluminaBasecallsToSamTool.defaults['num_processors']), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_illumina_demux)
    return parser


def main_illumina_demux(args):
    ''' Read Illumina runs & produce BAM files, demultiplexing to one bam per sample, or 
        for simplex runs, a single bam will be produced bearing the flowcell ID.
        Wraps together Picard's ExtractBarcodes (for multiplexed samples) and IlluminaBasecallsToSam
        while handling the various required input formats. Also can
        read Illumina BCL directories, tar.gz BCL directories.
    '''

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
        samples = SampleSheet(args.sampleSheet, only_lane=args.lane, append_run_id=run_id)
    else:
        samples = illumina.get_SampleSheet(only_lane=args.lane, append_run_id=run_id)


    link_locs=False
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
            illumina.get_BCLdir(),
            args.lane,
            read_structure,
            link_locs=link_locs
        )
    except subprocess.CalledProcessError as e:
        log.warning("CheckIlluminaDirectory failed for %s", illumina.get_BCLdir())
        if os.path.exists(os.path.join(illumina.get_intensities_dir(), "s.locs")):
            # recurse to remove broken links in directory
            log.info("This run has an 's.locs' file; checking for and removing broken per-tile symlinks...")
            broken_links = util.file.find_broken_symlinks(illumina.get_intensities_dir())
            if len(broken_links):
                for lpath in broken_links:
                    log.info("Removing broken symlink: %s", lpath)
                    os.unlink(lpath)

            # call CheckIlluminaDirectory with LINK_LOCS=true
            link_locs=True

            log.info("Checking run directory with Picard...")
            tools.picard.CheckIlluminaDirectoryTool().execute(
                illumina.get_BCLdir(),
                args.lane,
                read_structure,
                link_locs=link_locs
            )
        else:
            log.error("CheckIlluminaDirectory failed for %s", illumina.get_BCLdir())

    multiplexed_samples = True if 'B' in read_structure else False            
    
    if multiplexed_samples:
        assert samples is not None, "This looks like a multiplexed run since 'B' is in the read_structure: a SampleSheet must be given."
    else:
        assert samples==None, "A SampleSheet may not be provided unless 'B' is present in the read_structure"
        if args.commonBarcodes:
            log.warning("--commonBarcodes was set but 'B' is not present in the read_structure; emitting an empty file.")
            util.file.touch(args.commonBarcodes)

    # B in read structure indicates barcoded multiplexed samples
    if multiplexed_samples:
        # Picard ExtractIlluminaBarcodes
        extract_input = util.file.mkstempfname('.txt', prefix='.'.join(['barcodeData', flowcell, str(args.lane)]))
        barcodes_tmpdir = tempfile.mkdtemp(prefix='extracted_barcodes-')
        samples.make_barcodes_file(extract_input)
        out_metrics = (args.outMetrics is None) and util.file.mkstempfname('.metrics.txt') or args.outMetrics
        picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
                          if hasattr(args, opt) and getattr(args, opt) != None)
        picardOpts['read_structure'] = read_structure
        tools.picard.ExtractIlluminaBarcodesTool().execute(
            illumina.get_BCLdir(),
            args.lane,
            extract_input,
            barcodes_tmpdir,
            out_metrics,
            picardOptions=picardOpts,
            JVMmemory=args.JVMmemory)

        if args.commonBarcodes:
            # this step can take > 2 hours on a large high-output flowcell
            # so kick it to the background while we demux
            #count_and_sort_barcodes(barcodes_tmpdir, args.commonBarcodes)
            executor = concurrent.futures.ProcessPoolExecutor()
            executor.submit(count_and_sort_barcodes, barcodes_tmpdir, args.commonBarcodes, truncateToLength=args.max_barcodes, threads=util.misc.sanitize_thread_count(args.threads))

        # Picard IlluminaBasecallsToSam
        basecalls_input = util.file.mkstempfname('.txt', prefix='.'.join(['library_params', flowcell, str(args.lane)]))
        samples.make_params_file(args.outDir, basecalls_input)

    picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.IlluminaBasecallsToSamTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    picardOpts['run_start_date'] = run_date
    picardOpts['read_structure'] = read_structure
    if args.threads:
        picardOpts['num_processors'] = args.threads
    if not picardOpts.get('sequencing_center') and runinfo:
        picardOpts['sequencing_center'] = runinfo.get_machine()

    if picardOpts.get('sequencing_center'):
        picardOpts["sequencing_center"] = util.file.string_to_file_name(picardOpts["sequencing_center"])

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
            JVMmemory=args.JVMmemory)
    else:
        tools.picard.IlluminaBasecallsToSamTool().execute_single_sample(
            illumina.get_BCLdir(),
            os.path.join(args.outDir,flowcell+".bam"),
            flowcell,
            args.lane,
            flowcell,
            picardOptions=picardOpts,
            JVMmemory=args.JVMmemory)

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


__commands__.append(('illumina_demux', parser_illumina_demux))


# ==========================
# ***  lane_metrics   ***
# ==========================

def parser_lane_metrics(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.')
    parser.add_argument('outPrefix', help='''Prefix path to the *.illumina_lane_metrics and *.illumina_phasing_metrics files.''')
    parser.add_argument('--read_structure',
                        help='Override read structure (default: read from RunInfo.xml).',
                        default=None)
    parser.add_argument('--JVMmemory',
                        help='JVM virtual memory size (default: %(default)s)',
                        default=tools.picard.ExtractIlluminaBarcodesTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_lane_metrics)
    return parser

def main_lane_metrics(args):
    '''
        Write out lane metrics to a tsv file.
    '''
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

    picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.CollectIlluminaLaneMetricsTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    picardOpts['read_structure'] = read_structure
    tools.picard.CollectIlluminaLaneMetricsTool().execute(
        illumina.path,
        output_dir,
        output_prefix,
        picardOptions=picardOpts,
        JVMmemory=args.JVMmemory)

    illumina.close()
    return 0

__commands__.append(('lane_metrics', parser_lane_metrics))


# ==========================
# ***  common_barcodes   ***
# ==========================

def parser_common_barcodes(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Illumina BCL directory (or tar.gz of BCL directory). This is the top-level run directory.')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('outSummary', help='''Path to the summary file (.tsv format). It includes several columns: 
                                            (barcode1, likely_index_name1, barcode2, likely_index_name2, count), 
                                            where likely index names are either the exact match index name for the barcode 
                                            sequence, or those Hamming distance of 1 away.''')

    parser.add_argument('--truncateToLength',
                        help='If specified, only this number of barcodes will be returned. Useful if you only want the top N barcodes.',
                        type=int,
                        default=None)
    parser.add_argument('--omitHeader', 
                        help='If specified, a header will not be added to the outSummary tsv file.',
                        action='store_true')
    parser.add_argument('--includeNoise', 
                        help='If specified, barcodes with periods (".") will be included.',
                        action='store_true')
    parser.add_argument('--outMetrics',
                        help='Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file.',
                        default=None)
    parser.add_argument('--sampleSheet',
                        default=None,
                        help='''Override SampleSheet. Input tab or CSV file w/header and four named columns:
                                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2.
                                Default is to look for a SampleSheet.csv in the inDir.''')
    parser.add_argument('--flowcell', help='Override flowcell ID (default: read from RunInfo.xml).', default=None)
    parser.add_argument('--read_structure',
                        help='Override read structure (default: read from RunInfo.xml).',
                        default=None)

    for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list:
        if opt not in ('read_structure', 'num_processors'):
            parser.add_argument('--' + opt,
                                help='Picard ExtractIlluminaBarcodes ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt))

    parser.add_argument('--JVMmemory',
                        help='JVM virtual memory size (default: %(default)s)',
                        default=tools.picard.ExtractIlluminaBarcodesTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_common_barcodes)
    return parser

def main_common_barcodes(args):
    ''' 
        Extract Illumina barcodes for a run and write a TSV report 
        of the barcode counts in descending order
    '''

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
    barcode_file = util.file.mkstempfname('.txt', prefix='.'.join(['barcodeData', flowcell, str(args.lane)]))
    barcodes_tmpdir = tempfile.mkdtemp(prefix='extracted_barcodes-')
    samples.make_barcodes_file(barcode_file)
    out_metrics = (args.outMetrics is None) and util.file.mkstempfname('.metrics.txt') or args.outMetrics
    picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    picardOpts['read_structure'] = read_structure
    tools.picard.ExtractIlluminaBarcodesTool().execute(
        illumina.get_BCLdir(),
        args.lane,
        barcode_file,
        barcodes_tmpdir,
        out_metrics,
        picardOptions=picardOpts,
        JVMmemory=args.JVMmemory)

    count_and_sort_barcodes(barcodes_tmpdir, args.outSummary, args.truncateToLength, args.includeNoise, args.omitHeader)

    # clean up
    os.unlink(barcode_file)
    shutil.rmtree(barcodes_tmpdir)
    illumina.close()
    return 0

__commands__.append(('common_barcodes', parser_common_barcodes))

def count_and_sort_barcodes(barcodes_dir, outSummary, truncateToLength=None, includeNoise=False, omitHeader=False, threads=None):
    # collect the barcode file paths for all tiles
    tile_barcode_files = [os.path.join(barcodes_dir, filename) for filename in os.listdir(barcodes_dir)]

    # count all of the barcodes present in the tile files
    log.info("reading barcodes in all tile files")

    with util.file.CountDB() as reduce_db:
        barcodefile_tempfile_tuples = [(tile_barcode_file,util.file.mkstempfname('sqlite_.db')) for tile_barcode_file in tile_barcode_files]

        # scatter tile-specific barcode files among workers to store barcode counts in SQLite
        workers = util.misc.sanitize_thread_count(threads)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(util.file.count_occurrences_in_tsv_sqlite_backed, tf, bf, include_noise=includeNoise) for bf,tf in barcodefile_tempfile_tuples]
            for future in concurrent.futures.as_completed(futures):
                tmp_db, barcode_file = future.result()
                log.debug("done reading barcodes from %s; adding to total...",barcode_file)
                # gather and reduce counts from separate SQLite databases into one 
                reduce_db.add_counts_from_other_db(tmp_db)
                os.unlink(tmp_db)

        illumina_reference = IlluminaIndexReference()

        log.info("Number of barcodes seen %s",reduce_db.get_num_IDS())

        # write the barcodes and their corresponding counts
        with open(outSummary, 'w') as tsvfile:
            log.info("sorting counts...")
            log.info("writing output...")
            writer = csv.writer(tsvfile, delimiter='\t')
            # write the header unless the user has specified not to do so
            if not omitHeader:
                writer.writerow(("Barcode1", "Likely_Index_Names1", "Barcode2", "Likely_Index_Names2", "Count"))

            for num_processed,row in enumerate(reduce_db.get_counts_descending()):

                if truncateToLength and num_processed>truncateToLength:
                    break

                barcode,count = row

                writer.writerow((barcode[:8], ",".join([x for x in illumina_reference.guess_index(barcode[:8], distance=1)] or ["Unknown"]), 
                            barcode[8:], ",".join([x for x in illumina_reference.guess_index(barcode[8:], distance=1)] or ["Unknown"]), 
                            count))

                if num_processed%50000==0:
                    log.debug("written %s barcode summaries to output file",num_processed)

    log.info("done")

# ======================================
# ***  guess_low-abundance_barcodes  ***
# ======================================

def parser_guess_barcodes(parser=argparse.ArgumentParser()):
    parser.add_argument('in_barcodes', help='The barcode counts file produced by common_barcodes.')
    parser.add_argument('in_picard_metrics', help='The demultiplexing read metrics produced by Picard.')
    parser.add_argument('out_summary_tsv', help='''Path to the summary file (.tsv format). It includes several columns: 
                                            (sample_name, expected_barcode_1, expected_barcode_2, 
                                            expected_barcode_1_name, expected_barcode_2_name, 
                                            expected_barcodes_read_count, guessed_barcode_1, 
                                            guessed_barcode_2, guessed_barcode_1_name, 
                                            guessed_barcode_2_name, guessed_barcodes_read_count, 
                                            match_type), 
                                            where the expected values are those used by Picard during demultiplexing
                                            and the guessed values are based on the barcodes seen among the data.''')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--readcount_threshold',
                        default=None,
                        type=int,
                        help='''If specified, guess barcodes for samples with fewer than this many reads.''')
    group.add_argument('--sample_names',
                        nargs='*',
                        help='If specified, only guess barcodes for these sample names.',
                        type=str,
                        default=None)
    parser.add_argument('--outlier_threshold', 
                        help='threshold of how far from unbalanced a sample must be to be considered an outlier.',
                        type=float,
                        default=0.675)
    parser.add_argument('--expected_assigned_fraction', 
                        help='The fraction of reads expected to be assigned. An exception is raised if fewer than this fraction are assigned.',
                        type=float,
                        default=0.7)
    parser.add_argument('--number_of_negative_controls',
                        help='The number of negative controls in the pool, for calculating expected number of reads in the rest of the pool.',
                        type=int,
                        default=1)

    parser.add_argument('--rows_limit',
                        default=1000,
                        type=int,
                        help='''The number of rows to use from the in_barcodes.''')

    
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_guess_barcodes, split_args=True)
    return parser

def main_guess_barcodes(in_barcodes, 
                        in_picard_metrics, 
                        out_summary_tsv, 
                        sample_names, 
                        outlier_threshold, 
                        expected_assigned_fraction, 
                        number_of_negative_controls, 
                        readcount_threshold, 
                        rows_limit):
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

    bh = util.illumina_indices.IlluminaBarcodeHelper(in_barcodes, in_picard_metrics, rows_limit)
    guessed_barcodes = bh.find_uncertain_barcodes(sample_names=sample_names, 
                                                    outlier_threshold=outlier_threshold, 
                                                    expected_assigned_fraction=expected_assigned_fraction, 
                                                    number_of_negative_controls=number_of_negative_controls, 
                                                    readcount_threshold=readcount_threshold)
    bh.write_guessed_barcodes(out_summary_tsv, guessed_barcodes)

__commands__.append(('guess_barcodes', parser_guess_barcodes))


# ============================
# ***  IlluminaDirectory   ***
# ============================


class IlluminaDirectory(object):
    ''' A class that handles Illumina data directories
    '''

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
            if '://' in self.uri:
                raise NotImplementedError('boto s3 download here uri -> tarball')
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
            if os.path.isdir(os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')):
                # found it! self.path is correct
                break
            else:
                subdirs = list(os.path.join(self.path, x) for x in os.listdir(self.path)
                               if os.path.isdir(os.path.join(self.path, x)))
                if len(subdirs) == 1:
                    # follow the rabbit hole
                    self.path = subdirs[0]
                else:
                    # don't know where to go now!
                    raise Exception('cannot find Data/Intensities/BaseCalls/ inside %s (%s)' % (self.uri, self.path))

    def _extract_tarball(self, tarfile):
        self.tempDir = tempfile.mkdtemp(prefix='IlluminaDirectory-')
        self.path = self.tempDir
        util.file.extract_tarball(tarfile, self.tempDir)

    def close(self):
        if self.tempDir:
            shutil.rmtree(self.tempDir)
            self.tempDir = None

    def get_RunInfo(self):
        if self.runinfo is None:
            runinfo_file = os.path.join(self.path, 'RunInfo.xml')
            util.file.check_paths(runinfo_file)
            self.runinfo = RunInfo(runinfo_file)
        return self.runinfo

    def get_SampleSheet(self, only_lane=None, append_run_id=None):
        if self.samplesheet is None:
            samplesheet_file = os.path.join(self.path, 'SampleSheet.csv')
            util.file.check_paths(samplesheet_file)
            self.samplesheet = SampleSheet(samplesheet_file, only_lane=only_lane, append_run_id=append_run_id)
        return self.samplesheet

    def get_intensities_dir(self):
        return os.path.join(self.path, 'Data', 'Intensities')

    def get_BCLdir(self):
        return os.path.join(self.get_intensities_dir(), 'BaseCalls')


# ==================
# ***  RunInfo   ***
# ==================


class RunInfo(object):
    ''' A class that reads the RunInfo.xml file emitted by Illumina
        MiSeq and HiSeq machines.
    '''

    def __init__(self, xml_fname):
        self.fname = xml_fname
        self.root = xml.etree.ElementTree.parse(xml_fname).getroot()

    def get_fname(self):
        return self.fname

    def get_flowcell(self):
        fc = self.root[0].find('Flowcell').text
        # slice in the case where the ID has a prefix of zeros
        if re.match(r"^0+-", fc):
            if '-' in fc:
                # miseq often adds a bunch of leading zeros and a dash in front
                fc = "-".join(fc.split('-')[1:])
        # >=5 to avoid an exception here: https://github.com/broadinstitute/picard/blob/2.17.6/src/main/java/picard/illumina/IlluminaBasecallsToSam.java#L510
        # <= 15 to limit the bytes added to each bam record
        assert len(fc) >= 5,"The flowcell ID must be five or more characters in length"
        if len(fc) > 15:
            log.warning("The provided flowcell ID is longer than 15 characters. Is that correct?")
        return fc

    def _get_rundate_obj(self):
        """
            Access the text of the <Date> node in the RunInfo.xml file
            and returns an arrow date object.
        """
        rundate = self.root[0].find('Date').text
        # possible formats found in RunInfo.xml:
        #   "170712" (YYMMDD)
        #   "20170712" (YYYYMMDD)
        #   "6/27/2018 4:59:20 PM" (M/D/YYYY h:mm:ss A)
        datestring_formats = [
            "YYMMDD",
            "YYYYMMDD",
            "M/D/YYYY h:mm:ss A"
        ]
        for datestring_format in datestring_formats:
            try:
                date_parsed = arrow.get(rundate, datestring_format)
                return date_parsed
            except arrow.parser.ParserError:
                pass
        raise arrow.parser.ParserError("The date string seen in RunInfo.xml ('%s') did not match known Illumina formats: %s" % (rundate,datestring_formats) )

    def get_rundate_american(self):
        return str(self._get_rundate_obj().format("MM/DD/YYYY"))

    def get_rundate_iso(self):
        return str(self._get_rundate_obj().format("YYYY-MM-DD"))

    def get_machine(self):
        return self.root[0].find('Instrument').text

    def get_read_structure(self):
        reads = []
        for x in self.root[0].find('Reads').findall('Read'):
            order = int(x.attrib['Number'])
            read = x.attrib['NumCycles'] + (x.attrib['IsIndexedRead'] == 'Y' and 'B' or 'T')
            reads.append((order, read))
        return ''.join([r for _, r in sorted(reads)])

    def num_reads(self):
        return sum(1 for x in self.root[0].find('Reads').findall('Read') if x.attrib['IsIndexedRead'] == 'N')

# ======================
# ***  SampleSheet   ***
# ======================

class SampleSheetError(Exception):
    def __init__(self, message, fname):
        super(SampleSheetError, self).__init__(
            'Failed to read SampleSheet {}. {}'.format(
                fname, message))

class SampleSheet(object):
    ''' A class that reads an Illumina SampleSheet.csv or alternative/simplified
        tab-delimited versions as well.
    '''

    def __init__(self, infile, use_sample_name=True, only_lane=None, allow_non_unique=False, append_run_id=None):
        self.fname = infile
        self.use_sample_name = use_sample_name
        if only_lane is not None:
            only_lane = str(only_lane)
        self.only_lane = only_lane
        self.allow_non_unique = allow_non_unique
        self.append_run_id = append_run_id
        self.rows = []
        self._detect_and_load_sheet(infile)

    def _detect_and_load_sheet(self, infile):
        if infile.endswith(('.csv','.csv.gz')):
            # one of a few possible CSV formats (watch out for line endings from other OSes)
            with util.file.open_or_gzopen(infile, 'rU') as inf:
                header = None
                miseq_skip = False
                row_num = 0
                for line_no, line in enumerate(inf):
                    if line_no==0:
                        # remove BOM, if present
                        line = line.replace('\ufeff','')

                    # if this is a blank line, skip parsing and continue to the next line...
                    if len(line.rstrip('\r\n').strip()) == 0:
                        continue
                    csv.register_dialect('samplesheet', quoting=csv.QUOTE_MINIMAL, escapechar='\\')
                    row = next(csv.reader([line.strip().rstrip('\n')], dialect="samplesheet"))
                    row = [item.strip() for item in row] # remove leading/trailing whitespace from each item
                    if miseq_skip:
                        if line.startswith('[Data]'):
                            # start paying attention *after* this line
                            miseq_skip = False
                        # otherwise, skip all the miseq headers
                    elif line.startswith('['):
                        # miseq: ignore all lines until we see "[Data]"
                        miseq_skip = True
                    elif header is None:
                        header = row
                        if all(x in header for x in ['Sample_ID','Index']):
                            # this is a Broad Platform MiSeq-generated SampleSheet.csv
                            keymapper = {
                                'Sample_ID': 'sample',
                                'Index': 'barcode_1',
                                'Index2': 'barcode_2',
                                'Sample_Name': 'sample_name'
                            }
                            header = list(map(keymapper.get, header))
                        elif 'Sample_ID' in header:
                            # this is a MiSeq-generated SampleSheet.csv
                            keymapper = {
                                'Sample_ID': 'sample',
                                'index': 'barcode_1',
                                'index2': 'barcode_2',
                                'Sample_Name': 'sample_name'
                            }
                            header = list(map(keymapper.get, header))
                        elif 'SampleID' in header:
                            # this is a Broad Platform HiSeq-generated SampleSheet.csv
                            keymapper = {
                                'SampleID': 'sample',
                                'Index': 'barcode_1',
                                'Index2': 'barcode_2',
                                'libraryName': 'library_id_per_sample',
                                'FCID': 'flowcell',
                                'Lane': 'lane'
                            }
                            header = list(map(keymapper.get, header))
                        elif len(row) == 3:
                            # hopefully this is a Broad walk-up submission sheet (_web_iww_htdocs_seq...)
                            header = ['sample', 'barcode_1', 'barcode_2']
                            if 'sample' not in row[0].lower():
                                # this is an actual data row! (no header exists in this file)
                                row_num += 1
                                self.rows.append({
                                    'sample': row[0],
                                    'barcode_1': row[1],
                                    'barcode_2': row[2],
                                    'row_num': str(row_num)
                                })
                        else:
                            raise SampleSheetError('unrecognized filetype', infile)
                        for h in ('sample', 'barcode_1'):
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
                        row['row_num'] = str(row_num)
                        if (self.only_lane is not None and row.get('lane') and self.only_lane != row['lane']):
                            continue
                        if ('sample' in row and row['sample']) and ('barcode_1' in row and row['barcode_1']):
                            self.rows.append(row)
            # go back and re-shuffle miseq columns if use_sample_name applies
            if (self.use_sample_name and 'sample_name' in header and all(row.get('sample_name') for row in self.rows)):
                for row in self.rows:
                    row['library_id_per_sample'] = row['sample']
                    row['sample'] = row['sample_name']
            for row in self.rows:
                if 'sample_name' in row:
                    del row['sample_name']
        elif infile.endswith(('.txt','.txt.gz')):
            # our custom tab file format: sample, barcode_1, barcode_2, library_id_per_sample
            self.rows = []
            row_num = 0
            for row in util.file.read_tabfile_dict(infile):
                assert row.get('sample') and row.get('barcode_1')
                row_num += 1
                row['row_num'] = str(row_num)
                self.rows.append(row)
        else:
            raise SampleSheetError('unrecognized filetype', infile)

        if not self.rows:
            raise SampleSheetError('empty file', infile)

        # populate library IDs, run IDs (ie BAM filenames)
        for row in self.rows:
            row['library'] = row['sample']
            if row.get('library_id_per_sample'):
                row['library'] += '.l' + row['library_id_per_sample']
            row['run'] = row['library']
        if len(set(row['run'] for row in self.rows)) != len(self.rows):
            if self.allow_non_unique:
                log.warning("non-unique library IDs in this lane")
                unique_count = {}
                for row in self.rows:
                    unique_count.setdefault(row['library'], 0)
                    unique_count[row['library']] += 1
                    row['run'] += '.r' + str(unique_count[row['library']])
            else:
                raise SampleSheetError('non-unique library IDs in this lane', infile)
        if self.append_run_id:
            for row in self.rows:
                row['run'] += '.' + self.append_run_id

        # escape sample, run, and library IDs to be filename-compatible
        for row in self.rows:
            row['sample'] = util.file.string_to_file_name(row['sample'])
            row['library'] = util.file.string_to_file_name(row['library'])
            row['run'] = util.file.string_to_file_name(row['run'])

        # are we single or double indexed?
        if all(row.get('barcode_2') for row in self.rows):
            self.indexes = 2
        elif any(row.get('barcode_2') for row in self.rows):
            raise SampleSheetError('inconsistent single/double barcoding in sample sheet', infile)
        else:
            self.indexes = 1

    def make_barcodes_file(self, outFile):
        ''' Create input file for Picard ExtractBarcodes '''
        if self.num_indexes() == 2:
            header = ['barcode_name', 'library_name', 'barcode_sequence_1', 'barcode_sequence_2']
        else:
            header = ['barcode_name', 'library_name', 'barcode_sequence_1']
        with open(outFile, 'wt') as outf:
            outf.write('\t'.join(header) + '\n')
            for row in self.rows:
                out = {
                    'barcode_sequence_1': row['barcode_1'],
                    'barcode_sequence_2': row.get('barcode_2', ''),
                    'barcode_name': row['sample'],
                    'library_name': row['library']
                }
                outf.write('\t'.join(out[h] for h in header) + '\n')

    def make_params_file(self, bamDir, outFile):
        ''' Create input file for Picard IlluminaBasecallsToXXX '''
        if self.num_indexes() == 2:
            header = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1', 'BARCODE_2']
        else:
            header = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1']
        with open(outFile, 'wt') as outf:
            outf.write('\t'.join(header) + '\n')
            # add one catchall entry at the end called Unmatched
            rows = self.rows + [{
                'barcode_1': 'N',
                'barcode_2': 'N',
                'sample': 'Unmatched',
                'library': 'Unmatched',
                'run': 'Unmatched'
            }]
            for row in rows:
                out = {
                    'BARCODE_1': row['barcode_1'],
                    'BARCODE_2': row.get('barcode_2', ''),
                    'SAMPLE_ALIAS': row['sample'],
                    'LIBRARY_NAME': row['library']
                }
                out['OUTPUT'] = os.path.join(bamDir, row['run'] + ".bam")
                outf.write('\t'.join(out[h] for h in header) + '\n')

    def get_fname(self):
        return self.fname

    def get_rows(self):
        return self.rows

    def num_indexes(self):
        ''' Return 1 or 2 depending on whether pools are single or double indexed '''
        return self.indexes

    def fetch_by_index(self, idx):
        idx = str(idx)
        for row in self.rows:
            if idx == row['row_num']:
                return row
        return None

# =============================
# ***  miseq_fastq_to_bam   ***
# =============================


def miseq_fastq_to_bam(outBam, sampleSheet, inFastq1, inFastq2=None, runInfo=None,
                       sequencing_center=None,
                       JVMmemory=tools.picard.FastqToSamTool.jvmMemDefault):
    ''' Convert fastq read files to a single bam file. Fastq file names must conform
        to patterns emitted by Miseq machines. Sample metadata must be provided
        in a SampleSheet.csv that corresponds to the fastq filename. Specifically,
        the _S##_ index in the fastq file name will be used to find the corresponding
        row in the SampleSheet
    '''

    # match miseq based on fastq filenames
    mo = re.match(r"^\S+_S(\d+)_L001_R(\d)_001.fastq(?:.gz|)$", inFastq1)
    assert mo, "fastq filename %s does not match the patterns used by an Illumina Miseq machine" % inFastq1
    assert mo.group(2) == '1', "fastq1 must correspond to read 1, not read %s" % mo.group(2)
    sample_num = mo.group(1)
    if inFastq2:
        mo = re.match(r"^\S+_S(\d+)_L001_R(\d)_001.fastq(?:.gz|)$", inFastq2)
        assert mo, "fastq filename %s does not match the patterns used by an Illumina Miseq machine" % inFastq2
        assert mo.group(2) == '2', "fastq2 must correspond to read 2, not read %s" % mo.group(2)
        assert mo.group(1) == sample_num, "fastq1 (%s) and fastq2 (%s) must have the same sample number" % (
            sample_num, mo.group(1))

    # load metadata
    samples = SampleSheet(sampleSheet, allow_non_unique=True)
    sample_info = samples.fetch_by_index(sample_num)
    assert sample_info, "sample %s not found in %s" % (sample_num, sampleSheet)
    sampleName = sample_info['sample']
    log.info("Using sample name: %s", sampleName)
    if sample_info.get('barcode_2'):
        barcode = '-'.join((sample_info['barcode_1'], sample_info['barcode_2']))
    else:
        barcode = sample_info['barcode_1']
    picardOpts = {
        'LIBRARY_NAME': sample_info['library'],
        'PLATFORM': 'illumina',
        'VERBOSITY': 'WARNING',
        'QUIET': 'TRUE',
    }
    if runInfo:
        runInfo = RunInfo(runInfo)
        flowcell = runInfo.get_flowcell()
        picardOpts['RUN_DATE'] = runInfo.get_rundate_iso()
        if inFastq2:
            assert runInfo.num_reads() == 2, "paired fastqs given for a single-end RunInfo.xml"
        else:
            assert runInfo.num_reads() == 1, "second fastq missing for a paired-end RunInfo.xml"
    else:
        flowcell = 'A'
    if sequencing_center is None and runInfo:
        sequencing_center = runInfo.get_machine()
    if sequencing_center:
        picardOpts['SEQUENCING_CENTER'] = util.file.string_to_file_name(sequencing_center)
    picardOpts['PLATFORM_UNIT'] = '.'.join((flowcell, '1', barcode))
    if len(flowcell) > 5:
        flowcell = flowcell[:5]
    picardOpts['READ_GROUP_NAME'] = flowcell

    # run Picard
    picard = tools.picard.FastqToSamTool()
    picard.execute(inFastq1,
                   inFastq2,
                   sampleName,
                   outBam,
                   picardOptions=picard.dict_to_picard_opts(picardOpts),
                   JVMmemory=JVMmemory)
    return 0


def parser_miseq_fastq_to_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('outBam', help='Output BAM file.')
    parser.add_argument('sampleSheet', help='Input SampleSheet.csv file.')
    parser.add_argument('inFastq1', help='Input fastq file; 1st end of paired-end reads if paired.')
    parser.add_argument('--inFastq2', help='Input fastq file; 2nd end of paired-end reads.', default=None)
    parser.add_argument('--runInfo', help='Input RunInfo.xml file.', default=None)
    parser.add_argument(
        '--sequencing_center',
        default=None,
        help='Name of your sequencing center (default is the sequencing machine ID from the RunInfo.xml)')
    parser.add_argument('--JVMmemory',
                        default=tools.picard.FastqToSamTool.jvmMemDefault,
                        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, miseq_fastq_to_bam, split_args=True)
    return parser


__commands__.append(('miseq_fastq_to_bam', parser_miseq_fastq_to_bam))



# ==============================
# ***  extract_fc_metadata   ***
# ==============================
def extract_fc_metadata(flowcell, outRunInfo, outSampleSheet):
    ''' Extract RunInfo.xml and SampleSheet.csv from the provided Illumina directory
    '''
    illumina = IlluminaDirectory(flowcell)
    illumina.load()
    shutil.copy(illumina.get_RunInfo().get_fname(), outRunInfo)
    shutil.copy(illumina.get_SampleSheet().get_fname(), outSampleSheet)
    return 0
def parser_extract_fc_metadata(parser=argparse.ArgumentParser()):
    parser.add_argument('flowcell', help='Illumina directory (possibly tarball)')
    parser.add_argument('outRunInfo', help='Output RunInfo.xml file.')
    parser.add_argument('outSampleSheet', help='Output SampleSheet.csv file.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, extract_fc_metadata, split_args=True)
    return parser
__commands__.append(('extract_fc_metadata', parser_extract_fc_metadata))


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
