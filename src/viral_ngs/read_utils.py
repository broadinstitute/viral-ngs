#!/usr/bin/env python3
"""
Utilities for working with sequence reads, such as converting formats and
fixing mate pairs.
"""

__author__ = "irwin@broadinstitute.org, dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import math
import os
import tempfile
import shutil
import sys
import concurrent.futures
from contextlib import contextmanager
import functools

from Bio import SeqIO
import pysam

import util.cmd
import util.file
import util.misc
from util.file import mkstempfname
import tools.bwa
import tools.cdhit
import tools.picard
import tools.samtools
import tools.mvicuna
import tools.prinseq
import tools.novoalign
import tools.gatk

log = logging.getLogger(__name__)

# =======================
# ***  purge_unmated  ***
# =======================


def purge_unmated(inFastq1, inFastq2, outFastq1, outFastq2, regex=r'^@(\S+)/[1|2]$'):
    '''Use mergeShuffledFastqSeqs to purge unmated reads, and
       put corresponding reads in the same order.
       Corresponding sequences must have sequence identifiers
       of the form SEQID/1 and SEQID/2.
    '''
    tempOutput = mkstempfname()
    mergeShuffledFastqSeqsPath = os.path.join(util.file.get_scripts_path(), 'mergeShuffledFastqSeqs.pl')
    cmdline = [mergeShuffledFastqSeqsPath, '-t', '-r', regex, '-f1', inFastq1, '-f2', inFastq2, '-o', tempOutput]
    log.debug(' '.join(cmdline))
    util.misc.run_and_print(cmdline, check=True)
    shutil.move(tempOutput + '.1.fastq', outFastq1)
    shutil.move(tempOutput + '.2.fastq', outFastq2)
    return 0


def parser_purge_unmated(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq1', help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2', help='Input fastq file; 2nd end of paired-end reads.')
    parser.add_argument('outFastq1', help='Output fastq file; 1st end of paired-end reads.')
    parser.add_argument('outFastq2', help='Output fastq file; 2nd end of paired-end reads.')
    parser.add_argument(
        "--regex",
        help="Perl regular expression to parse paired read IDs (default: %(default)s)",
        default=r'^@(\S+)/[1|2]$'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, purge_unmated, split_args=True)
    return parser


__commands__.append(('purge_unmated', parser_purge_unmated))


# ===============================
# ***  index_fasta_samtools   ***
# ===============================


def parser_index_fasta_samtools(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta', help='Reference genome, FASTA format.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_index_fasta_samtools)
    return parser


def main_index_fasta_samtools(args):
    '''Index a reference genome for Samtools.'''
    tools.samtools.SamtoolsTool().faidx(args.inFasta, overwrite=True)
    return 0


__commands__.append(('index_fasta_samtools', parser_index_fasta_samtools))

# =============================
# ***  index_fasta_picard   ***
# =============================


def parser_index_fasta_picard(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta', help='Input reference genome, FASTA format.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.CreateSequenceDictionaryTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s CreateSequenceDictionary, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_index_fasta_picard)
    return parser


def main_index_fasta_picard(args):
    '''Create an index file for a reference genome suitable for Picard/GATK.'''
    tools.picard.CreateSequenceDictionaryTool().execute(
        args.inFasta, overwrite=True,
        picardOptions=args.picardOptions,
        JVMmemory=args.JVMmemory
    )
    return 0


__commands__.append(('index_fasta_picard', parser_index_fasta_picard))

# =============================
# ***  mkdup_picard   ***
# =============================


def parser_mkdup_picard(parser=argparse.ArgumentParser()):
    parser.add_argument('inBams', help='Input reads, BAM format.', nargs='+')
    parser.add_argument('outBam', help='Output reads, BAM format.')
    parser.add_argument('--outMetrics', help='Output metrics file. Default is to dump to a temp file.', default=None)
    parser.add_argument(
        "--remove",
        help="Instead of marking duplicates, remove them entirely (default: %(default)s)",
        default=False,
        action="store_true",
        dest="remove"
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.MarkDuplicatesTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s MarkDuplicates, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_mkdup_picard)
    return parser


def main_mkdup_picard(args):
    '''Mark or remove duplicate reads from BAM file.'''
    opts = list(args.picardOptions)
    if args.remove:
        opts = ['REMOVE_DUPLICATES=true'] + opts
    tools.picard.MarkDuplicatesTool().execute(
        args.inBams, args.outBam,
        args.outMetrics, picardOptions=opts,
        JVMmemory=args.JVMmemory
    )
    return 0


__commands__.append(('mkdup_picard', parser_mkdup_picard))

# =============================
# ***  revert_bam_picard   ***
# =============================

def parser_revert_sam_common(parser=argparse.ArgumentParser()):
    parser.add_argument('--clearTags', dest='clear_tags', default=False, action='store_true', 
                        help='When supplying an aligned input file, clear the per-read attribute tags')
    parser.add_argument("--tagsToClear", type=str, nargs='+', dest="tags_to_clear", default=["XT", "X0", "X1", "XA", 
                        "AM", "SM", "BQ", "CT", "XN", "OC", "OP"], 
                        help='A space-separated list of tags to remove from all reads in the input bam file (default: %(default)s)')
    parser.add_argument(
        '--doNotSanitize',
        dest="do_not_sanitize",
        action='store_true',
        help="""When being reverted, picard's SANITIZE=true
                is set unless --doNotSanitize is given. 
                Sanitization is a destructive operation that removes reads so the bam file is consistent.
                From the picard documentation:
                    'Reads discarded include (but are not limited to) paired reads with missing mates, duplicated records, 
                    records with mismatches in length of bases and qualities.'
                For more information see:
                    https://broadinstitute.github.io/picard/command-line-overview.html#RevertSam
                """
    )
    try:
        parser.add_argument(
            '--JVMmemory',
            default=tools.picard.RevertSamTool.jvmMemDefault,
            help='JVM virtual memory size (default: %(default)s)'
        )
    except argparse.ArgumentError:
        # if --JVMmemory is already present in the parser, continue on...
        pass
    return parser

def parser_revert_bam_picard(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('outBam', help='Output reads, BAM format.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.RevertSamTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s RevertSam, OPTIONNAME=value ...'
    )
    parser = parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_revert_bam_picard, split_args=True)
    return parser


def main_revert_bam_picard(inBam, outBam, clear_tags=False, tags_to_clear=None, picardOptions=None, do_not_sanitize=False, JVMmemory=None):
    '''Revert BAM to raw reads'''
    picardOptions = picardOptions or []
    tags_to_clear = tags_to_clear or []

    if clear_tags:
        for tag in tags_to_clear:
            picardOptions.append("ATTRIBUTE_TO_CLEAR={}".format(tag))
            
            if not do_not_sanitize and not any(opt.startswith("SANITIZE") for opt in picardOptions):
                picardOptions.append('SANITIZE=true')

    tools.picard.RevertSamTool().execute(inBam, outBam, picardOptions=picardOptions, JVMmemory=JVMmemory)
    return 0
__commands__.append(('revert_bam_picard', parser_revert_bam_picard))

@contextmanager
def revert_bam_if_aligned(inBam, revert_bam=None, clear_tags=True, tags_to_clear=None, picardOptions=None, JVMmemory="4g", sanitize=False):
    '''
        revert the bam file if aligned. 
        If a revertBam file path is specified, it is used, otherwise a temp file is created.
    '''

    revertBamOut = revert_bam if revert_bam else util.file.mkstempfname('.bam')
    picardOptions = picardOptions or []
    tags_to_clear = tags_to_clear or []

    outBam = inBam
    with pysam.AlignmentFile(inBam, 'rb', check_sq=False) as bam:
        # if it looks like the bam is aligned, revert it
        if 'SQ' in bam.header and len(bam.header['SQ'])>0:
            if clear_tags:
                for tag in tags_to_clear:
                    picardOptions.append("ATTRIBUTE_TO_CLEAR={}".format(tag))

            if not any(opt.startswith("SORT_ORDER") for opt in picardOptions):
                picardOptions.append('SORT_ORDER=queryname')

            if sanitize and not any(opt.startswith("SANITIZE") for opt in picardOptions):
                picardOptions.append('SANITIZE=true')

            tools.picard.RevertSamTool().execute(
                inBam, revertBamOut, picardOptions=picardOptions 
            )
            outBam = revertBamOut
        else:
            # if we don't need to produce a revertBam file
            # but the user has specified one anyway
            # simply touch the output
            if revert_bam:
                log.warning("An output was specified for 'revertBam', but the input is unaligned, so RevertSam was not needed. Touching the output.")
                util.file.touch(revertBamOut)
                # TODO: error out? run RevertSam anyway?

    yield outBam

    if revert_bam == None:
        os.unlink(revertBamOut)

# =========================
# ***  generic picard   ***
# =========================


def parser_picard(parser=argparse.ArgumentParser()):
    parser.add_argument('command', help='picard command')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.PicardTools.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[], nargs='*',
        help='Optional arguments to Picard, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_picard)
    return parser


def main_picard(args):
    '''Generic Picard runner.'''
    tools.picard.PicardTools().execute(args.command, picardOptions=args.picardOptions, JVMmemory=args.JVMmemory)
    return 0


__commands__.append(('picard', parser_picard))

# ===================
# ***  sort_bam   ***
# ===================


def parser_sort_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input bam file.')
    parser.add_argument('outBam', help='Output bam file, sorted.')
    parser.add_argument(
        'sortOrder',
        help='How to sort the reads. [default: %(default)s]',
        choices=tools.picard.SortSamTool.valid_sort_orders,
        default=tools.picard.SortSamTool.default_sort_order
    )
    parser.add_argument(
        "--index",
        help="Index outBam (default: %(default)s)",
        default=False,
        action="store_true",
        dest="index"
    )
    parser.add_argument(
        "--md5",
        help="MD5 checksum outBam (default: %(default)s)",
        default=False,
        action="store_true",
        dest="md5"
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.SortSamTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s SortSam, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_sort_bam)
    return parser


def main_sort_bam(args):
    '''Sort BAM file'''
    opts = list(args.picardOptions)
    if args.index:
        opts = ['CREATE_INDEX=true'] + opts
    if args.md5:
        opts = ['CREATE_MD5_FILE=true'] + opts
    tools.picard.SortSamTool().execute(
        args.inBam, args.outBam, args.sortOrder,
        picardOptions=opts, JVMmemory=args.JVMmemory
    )
    return 0


__commands__.append(('sort_bam', parser_sort_bam))

# ====================
# ***  downsample_bams  ***
# ====================

def parser_downsample_bams(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bams', help='Input bam files.', nargs='+')
    parser.add_argument('--outPath', dest="out_path", type=str, help="""Output path. If not provided, 
                        downsampled bam files will be written to the same paths as each 
                        source bam file""")
    parser.add_argument('--readCount', dest="specified_read_count", type=int, help='The number of reads to downsample to.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--deduplicateBefore', action='store_true', dest="deduplicate_before", help='de-duplicate reads before downsampling.')
    group.add_argument('--deduplicateAfter', action='store_true', dest="deduplicate_after", help='de-duplicate reads after downsampling.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.DownsampleSamTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s DownsampleSam, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_downsample_bams, split_args=True)
    return parser


def main_downsample_bams(in_bams, out_path, specified_read_count=None, deduplicate_before=False, deduplicate_after=False, picardOptions=None, threads=None, JVMmemory=None):
    '''Downsample multiple bam files to the smallest read count in common, or to the specified count.'''
    if picardOptions is None:
        picardOptions = []

    opts = list(picardOptions) + []

    def get_read_counts(bams):
        samtools = tools.samtools.SamtoolsTool()
        # get read counts for bam files provided
        read_counts = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=util.misc.sanitize_thread_count(threads)) as executor:
            for x, result in zip(bams, executor.map(samtools.count, bams)):
                read_counts[x] = int(result)

        return read_counts

    def get_downsample_target_count(bams, readcount_requested):
        read_counts = get_read_counts(bams)
        downsample_target = 0
        if readcount_requested is not None:
            downsample_target = readcount_requested
            if downsample_target > min(read_counts.values()):
                raise ValueError("The smallest file has %s reads, which is less than the downsample target specified, %s. Please reduce the target count or omit it to use the read count of the smallest input file." % (min(read_counts.values()), downsample_target))
        else:
            downsample_target = min(read_counts.values())
        return downsample_target

    def downsample_bams(data_pairs, downsample_target, threads=None, JVMmemory=None):
        downsamplesam = tools.picard.DownsampleSamTool()
        workers = util.misc.sanitize_thread_count(threads)
        JVMmemory = JVMmemory if JVMmemory else tools.picard.DownsampleSamTool.jvmMemDefault
        
        jvm_worker_memory_mb = str(int(util.misc.convert_size_str(JVMmemory,"m")[:-1])//workers)+"m"
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_file = {executor.submit(downsamplesam.downsample_to_approx_count, *(list(fp)+[downsample_target]), JVMmemory=jvm_worker_memory_mb): fp[0] for fp in data_pairs}
            for future in concurrent.futures.as_completed(future_to_file):
                f = future_to_file[future]
                try:
                    data = future.result()
                except Exception as exc:
                    raise

    def dedup_bams(data_pairs, threads=None, JVMmemory=None):
        workers = util.misc.sanitize_thread_count(threads)
        JVMmemory = JVMmemory if JVMmemory else tools.picard.DownsampleSamTool.jvmMemDefault
        jvm_worker_memory_mb = str(int(util.misc.convert_size_str(JVMmemory,"m")[:-1])//workers)+"m"
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_file = {executor.submit(rmdup_mvicuna_bam, *fp, JVMmemory=jvm_worker_memory_mb): fp[0] for fp in data_pairs}
            for future in concurrent.futures.as_completed(future_to_file):
                f = future_to_file[future]
                try:
                    data = future.result()
                except Exception as exc:
                    raise

    data_pairs = []
    if deduplicate_before:
        with util.file.tempfnames(suffixes=[ '{}.dedup.bam'.format(os.path.splitext(os.path.basename(x))[0]) for x in in_bams]) as deduped_bams:
            data_pairs = list(zip(in_bams, deduped_bams))
            dedup_bams(data_pairs, threads=threads, JVMmemory=JVMmemory)
            downsample_target = get_downsample_target_count(deduped_bams, specified_read_count)
            
            if out_path:
                util.file.mkdir_p(out_path)
                data_pairs = list(zip(deduped_bams, [os.path.join(out_path, os.path.splitext(os.path.basename(x))[0]+".dedup.downsampled-{}.bam".format(downsample_target)) for x in in_bams]))
            else:
                data_pairs = list(zip(deduped_bams, [os.path.splitext(x)[0]+".dedup.downsampled-{}.bam".format(downsample_target) for x in in_bams]))
            downsample_bams(data_pairs, downsample_target, threads=threads, JVMmemory=JVMmemory)
    else:
        downsample_target = get_downsample_target_count(in_bams, specified_read_count)
        if deduplicate_after:
            log.warning("--deduplicateAfter has been specified. Read count of output files is not guaranteed.")
            with util.file.tempfnames(suffixes=[ '{}.downsampled-{}.bam'.format(os.path.splitext(os.path.basename(x))[0], downsample_target) for x in in_bams]) as downsampled_tmp_bams:
                data_pairs = list(zip(in_bams, downsampled_tmp_bams))
                downsample_bams(data_pairs, downsample_target, threads=threads, JVMmemory=JVMmemory)
                
                if out_path:
                    util.file.mkdir_p(out_path)
                    data_pairs = list(zip(downsampled_tmp_bams, [os.path.join(out_path, os.path.splitext(os.path.basename(x))[0]+".downsampled-{}.dedup.bam").format(downsample_target) for x in in_bams]))
                else:
                    data_pairs = list(zip(downsampled_tmp_bams, [os.path.splitext(x)[0]+".downsampled-{}.dedup.bam".format(downsample_target) for x in in_bams]))
                dedup_bams(data_pairs, threads=threads, JVMmemory=JVMmemory)
        else:
            if out_path:
                util.file.mkdir_p(out_path)
                data_pairs = list(zip(in_bams, [os.path.join(out_path, os.path.splitext(os.path.basename(x))[0]+".downsampled-{}.bam".format(downsample_target)) for x in in_bams]))
            else:
                data_pairs = list(zip(in_bams, [os.path.splitext(x)[0]+".downsampled-{}.bam".format(downsample_target) for x in in_bams]))
            downsample_bams(data_pairs, downsample_target, threads=threads, JVMmemory=JVMmemory)
    return 0

__commands__.append(('downsample_bams', parser_downsample_bams))

# ====================
# ***  merge_bams  ***
# ====================


def parser_merge_bams(parser=argparse.ArgumentParser()):
    parser.add_argument('inBams', help='Input bam files.', nargs='+')
    parser.add_argument('outBam', help='Output bam file.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.MergeSamFilesTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s MergeSamFiles, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_merge_bams)
    return parser


def main_merge_bams(args):
    '''Merge multiple BAMs into one'''
    opts = list(args.picardOptions) + ['USE_THREADING=true']
    tools.picard.MergeSamFilesTool().execute(args.inBams, args.outBam, picardOptions=opts, JVMmemory=args.JVMmemory)
    return 0


__commands__.append(('merge_bams', parser_merge_bams))

# ====================
# ***  filter_bam  ***
# ====================


def parser_filter_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input bam file.')
    parser.add_argument('readList', help='Input file of read IDs.')
    parser.add_argument('outBam', help='Output bam file.')
    parser.add_argument(
        "--exclude",
        help="""If specified, readList is a list of reads to remove from input.
            Default behavior is to treat readList as an inclusion list (all unnamed
            reads are removed).""",
        default=False,
        action="store_true",
        dest="exclude"
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s FilterSamReads, OPTIONNAME=value ...'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_filter_bam)
    return parser


def main_filter_bam(args):
    '''Filter BAM file by read name'''
    tools.picard.FilterSamReadsTool().execute(
        args.inBam,
        args.exclude,
        args.readList,
        args.outBam,
        picardOptions=args.picardOptions,
        JVMmemory=args.JVMmemory
    )
    return 0


__commands__.append(('filter_bam', parser_filter_bam))



# =======================
# ***  fastq_to_bam   ***
# =======================


def fastq_to_bam(
    inFastq1,
    inFastq2,
    outBam,
    sampleName=None,
    header=None,
    JVMmemory=tools.picard.FastqToSamTool.jvmMemDefault,
    picardOptions=None
):
    ''' Convert a pair of fastq paired-end read files and optional text header
        to a single bam file.
    '''
    picardOptions = picardOptions or []

    if header:
        fastqToSamOut = mkstempfname('.bam')
    else:
        fastqToSamOut = outBam
    if sampleName is None:
        sampleName = 'Dummy'    # Will get overwritten by rehead command
    if header:
        # With the header option, rehead will be called after FastqToSam.
        # This will invalidate any md5 file, which would be a slow to construct
        # on our own, so just disallow and let the caller run md5sum if desired.
        if any(opt.lower() == 'CREATE_MD5_FILE=True'.lower() for opt in picardOptions):
            raise Exception("""CREATE_MD5_FILE is not allowed with '--header.'""")
    tools.picard.FastqToSamTool().execute(
        inFastq1, inFastq2,
        sampleName, fastqToSamOut,
        picardOptions=picardOptions,
        JVMmemory=JVMmemory
    )

    if header:
        tools.samtools.SamtoolsTool().reheader(fastqToSamOut, header, outBam)

    return 0


def parser_fastq_to_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq1', help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2', help='Input fastq file; 2nd end of paired-end reads.')
    parser.add_argument('outBam', help='Output bam file.')
    headerGroup = parser.add_mutually_exclusive_group(required=True)
    headerGroup.add_argument('--sampleName', help='Sample name to insert into the read group header.')
    headerGroup.add_argument('--header', help='Optional text file containing header.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FastqToSamTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='''Optional arguments to Picard\'s FastqToSam,
                OPTIONNAME=value ...  Note that header-related options will be
                overwritten by HEADER if present.'''
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, fastq_to_bam, split_args=True)
    return parser


__commands__.append(('fastq_to_bam', parser_fastq_to_bam))


def join_paired_fastq(
        inFastqs,
        output,
        outFormat
):
    ''' Join paired fastq reads into single reads with Ns
    '''
    inFastqs = list(inFastqs)
    if output == '-':
        output = sys.stdout
    SeqIO.write(util.file.join_paired_fastq(inFastqs, output_format=outFormat), output, outFormat)
    return 0


def parser_join_paired_fastq(parser=argparse.ArgumentParser()):
    parser.add_argument('output', help='Output file.')
    parser.add_argument('inFastqs', nargs='+', help='Input fastq file (2 if paired, 1 if interleaved)')
    parser.add_argument('--outFormat', default='fastq', help='Output file format.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, join_paired_fastq, split_args=True)
    return parser

__commands__.append(('join_paired_fastq', parser_join_paired_fastq))


# ======================
# ***  split_reads   ***
# ======================
defaultIndexLen = 2
defaultMaxReads = 1000
defaultFormat = 'fastq'


def split_bam(inBam, outBams):
    '''Split BAM file equally into several output BAM files. '''
    samtools = tools.samtools.SamtoolsTool()
    picard = tools.picard.PicardTools()

    # get totalReadCount and maxReads
    # maxReads = totalReadCount / num files, but round up to the nearest
    # even number in order to keep read pairs together (assuming the input
    # is sorted in query order and has no unmated reads, which can be
    # accomplished by Picard RevertSam with SANITIZE=true)
    totalReadCount = samtools.count(inBam)
    maxReads = int(math.ceil(float(totalReadCount) / len(outBams) / 2) * 2)
    log.info("splitting %d reads into %d files of %d reads each", totalReadCount, len(outBams), maxReads)

    # load BAM header into memory
    header = samtools.getHeader(inBam)
    if 'SO:queryname' not in header[0]:
        raise Exception('Input BAM file must be sorted in queryame order')

    # dump to bigsam
    bigsam = mkstempfname('.sam')
    samtools.view(['-@', '3'], inBam, bigsam)

    # split bigsam into little ones
    with util.file.open_or_gzopen(bigsam, 'rt') as inf:
        for outBam in outBams:
            log.info("preparing file " + outBam)
            tmp_sam_reads = mkstempfname('.sam')
            with open(tmp_sam_reads, 'wt') as outf:
                for row in header:
                    outf.write('\t'.join(row) + '\n')
                for _ in range(maxReads):
                    line = inf.readline()
                    if not line:
                        break
                    outf.write(line)
                if outBam == outBams[-1]:
                    for line in inf:
                        outf.write(line)
            picard.execute(
                "SamFormatConverter", [
                    'INPUT=' + tmp_sam_reads, 'OUTPUT=' + outBam, 'VERBOSITY=WARNING'
                ],
                JVMmemory='512m'
            )
            os.unlink(tmp_sam_reads)
    os.unlink(bigsam)


def parser_split_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument('outBams', nargs='+', help='Output BAM files')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, split_bam, split_args=True)
    return parser


__commands__.append(('split_bam', parser_split_bam))

# =======================
# ***  reheader_bam   ***
# =======================


def parser_reheader_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('rgMap', help='Tabular file containing three columns: field, old, new.')
    parser.add_argument('outBam', help='Output reads, BAM format.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_reheader_bam)
    return parser


def main_reheader_bam(args):
    ''' Copy a BAM file (inBam to outBam) while renaming elements of the BAM header.
        The mapping file specifies which (key, old value, new value) mappings. For
        example:
            LB	lib1	lib_one
            SM	sample1	Sample_1
            SM	sample2	Sample_2
            SM	sample3	Sample_3
            CN	broad	BI
    '''
    # read mapping file
    mapper = dict((a + ':' + b, a + ':' + c) for a, b, c in util.file.read_tabfile(args.rgMap))
    # read and convert bam header
    header_file = mkstempfname('.sam')
    with open(header_file, 'wt') as outf:
        for row in tools.samtools.SamtoolsTool().getHeader(args.inBam):
            if row[0] == '@RG':
                row = [mapper.get(x, x) for x in row]
            outf.write('\t'.join(row) + '\n')
    # write new bam with new header
    tools.samtools.SamtoolsTool().reheader(args.inBam, header_file, args.outBam)
    os.unlink(header_file)
    return 0


__commands__.append(('reheader_bam', parser_reheader_bam))


def parser_reheader_bams(parser=argparse.ArgumentParser()):
    parser.add_argument('rgMap', help='Tabular file containing three columns: field, old, new.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_reheader_bams)
    return parser


def main_reheader_bams(args):
    ''' Copy BAM files while renaming elements of the BAM header.
        The mapping file specifies which (key, old value, new value) mappings. For
        example:
            LB  lib1  lib_one
            SM  sample1 Sample_1
            SM  sample2 Sample_2
            SM  sample3 Sample_3
            CN  broad   BI
            FN  in1.bam out1.bam
            FN  in2.bam out2.bam
    '''
    # read mapping file
    mapper = dict((a + ':' + b, a + ':' + c) for a, b, c in util.file.read_tabfile(args.rgMap) if a != 'FN')
    files = list((b, c) for a, b, c in util.file.read_tabfile(args.rgMap) if a == 'FN')
    header_file = mkstempfname('.sam')
    # read and convert bam headers
    for inBam, outBam in files:
        if os.path.isfile(inBam):
            with open(header_file, 'wt') as outf:
                for row in tools.samtools.SamtoolsTool().getHeader(inBam):
                    if row[0] == '@RG':
                        row = [mapper.get(x, x) for x in row]
                    outf.write('\t'.join(row) + '\n')
            # write new bam with new header
            tools.samtools.SamtoolsTool().reheader(inBam, header_file, outBam)
    os.unlink(header_file)
    return 0


__commands__.append(('reheader_bams', parser_reheader_bams))

# ============================
# ***  dup_remove_mvicuna  ***
# ============================


def mvicuna_fastqs_to_readlist(inFastq1, inFastq2, readList):
    # Run M-Vicuna on FASTQ files
    outFastq1 = mkstempfname('.1.fastq')
    outFastq2 = mkstempfname('.2.fastq')
    if inFastq2 is None or os.path.getsize(inFastq2) < 10:
        tools.mvicuna.MvicunaTool().rmdup_single(inFastq1, outFastq1)
    else:
        tools.mvicuna.MvicunaTool().rmdup((inFastq1, inFastq2), (outFastq1, outFastq2), None)

    # Make a list of reads to keep
    with open(readList, 'at') as outf:
        for fq in (outFastq1, outFastq2):
            with util.file.open_or_gzopen(fq, 'rt') as inf:
                line_num = 0
                for line in inf:
                    if (line_num % 4) == 0:
                        idVal = line.rstrip('\n')[1:]
                        if idVal.endswith('/1'):
                            outf.write(idVal[:-2] + '\n')
                        # single-end reads do not have /1 /2 mate suffix
                        # so pass through their IDs
                        if not (idVal.endswith('/1') or idVal.endswith('/2')):
                            outf.write(idVal + '\n')
                    line_num += 1
    os.unlink(outFastq1)
    os.unlink(outFastq2)


def rmdup_cdhit_bam(inBam, outBam, max_mismatches=None, jvm_memory=None):
    ''' Remove duplicate reads from BAM file using cd-hit-dup.
    '''
    max_mismatches = max_mismatches or 4
    tmp_dir = tempfile.mkdtemp()

    tools.picard.SplitSamByLibraryTool().execute(inBam, tmp_dir)

    s2fq_tool = tools.picard.SamToFastqTool()
    cdhit = tools.cdhit.CdHit()
    out_bams = []
    for f in os.listdir(tmp_dir):
        out_bam = mkstempfname('.bam')
        out_bams.append(out_bam)
        library_sam = os.path.join(tmp_dir, f)

        in_fastqs = mkstempfname('.1.fastq'), mkstempfname('.2.fastq')

        s2fq_tool.execute(library_sam, in_fastqs[0], in_fastqs[1])
        if not os.path.getsize(in_fastqs[0]) > 0 and not os.path.getsize(in_fastqs[1]) > 0:
            continue

        out_fastqs = mkstempfname('.1.fastq'), mkstempfname('.2.fastq')
        options = {
            '-e': max_mismatches,
        }
        if in_fastqs[1] is not None and os.path.getsize(in_fastqs[1]) > 10:
            options['-i2'] = in_fastqs[1]
            options['-o2'] = out_fastqs[1]

        log.info("executing cd-hit-est on library " + library_sam)
        # cd-hit-dup cannot operate on piped fastq input because it reads twice
        cdhit.execute('cd-hit-dup', in_fastqs[0], out_fastqs[0], options=options, background=True)

        tools.picard.FastqToSamTool().execute(out_fastqs[0], out_fastqs[1], f, out_bam, JVMmemory=jvm_memory)
        for fn in in_fastqs:
            os.unlink(fn)

    with util.file.fifo(name='merged.sam') as merged_bam:
        merge_opts = ['SORT_ORDER=queryname']
        tools.picard.MergeSamFilesTool().execute(out_bams, merged_bam, picardOptions=merge_opts, JVMmemory=jvm_memory, background=True)
        tools.picard.ReplaceSamHeaderTool().execute(merged_bam, inBam, outBam, JVMmemory=jvm_memory)


def parser_rmdup_cdhit_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('outBam', help='Output reads, BAM format.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)',
        dest='jvm_memory'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, rmdup_cdhit_bam, split_args=True)
    return parser


__commands__.append(('rmdup_cdhit_bam', parser_rmdup_cdhit_bam))

def _merge_fastqs_and_mvicuna(lb, files):
    readList = mkstempfname('.keep_reads.txt')
    log.info("executing M-Vicuna DupRm on library " + lb)

    # create merged FASTQs per library
    infastqs = (mkstempfname('.1.fastq'), mkstempfname('.2.fastq'))
    for d in range(2):
        with open(infastqs[d], 'wt') as outf:
            for fprefix in files:
                fn = '%s_%d.fastq' % (fprefix, d + 1)

                if os.path.isfile(fn):
                    with open(fn, 'rt') as inf:
                        for line in inf:
                            outf.write(line)
                    os.unlink(fn)
                else:
                    log.warning(
                        """no reads found in %s,
                                assuming that's because there's no reads in that read group""", fn
                    )

    # M-Vicuna DupRm to see what we should keep (append IDs to running file)
    if os.path.getsize(infastqs[0]) > 0 or os.path.getsize(infastqs[1]) > 0:
        mvicuna_fastqs_to_readlist(infastqs[0], infastqs[1], readList)
    for fn in infastqs:
        os.unlink(fn)

    return readList

def rmdup_mvicuna_bam(inBam, outBam, JVMmemory=None):
    ''' Remove duplicate reads from BAM file using M-Vicuna. The
        primary advantage to this approach over Picard's MarkDuplicates tool
        is that Picard requires that input reads are aligned to a reference,
        and M-Vicuna can operate on unaligned reads.
    '''

    # Convert BAM -> FASTQ pairs per read group and load all read groups
    tempDir = tempfile.mkdtemp()
    tools.picard.SamToFastqTool().per_read_group(inBam, tempDir, picardOptions=['VALIDATION_STRINGENCY=LENIENT'], JVMmemory=JVMmemory)
    read_groups = [x[1:] for x in tools.samtools.SamtoolsTool().getHeader(inBam) if x[0] == '@RG']
    read_groups = [dict(pair.split(':', 1) for pair in rg) for rg in read_groups]

    # Collect FASTQ pairs for each library
    lb_to_files = {}
    for rg in read_groups:
        lb_to_files.setdefault(rg.get('LB', 'none'), set())
        fname = rg['ID']
        lb_to_files[rg.get('LB', 'none')].add(os.path.join(tempDir, fname))
    log.info("found %d distinct libraries and %d read groups", len(lb_to_files), len(read_groups))

    # For each library, merge FASTQs and run rmdup for entire library
    readListAll = mkstempfname('.keep_reads_all.txt')
    per_lb_read_lists = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=util.misc.available_cpu_count()) as executor:
        futures = [executor.submit(_merge_fastqs_and_mvicuna, lb, files) for lb, files in lb_to_files.items()]
        for future in concurrent.futures.as_completed(futures):
            log.info("mvicuna finished processing library")
            try:
                readList = future.result()
                per_lb_read_lists.append(readList)
            except Exception as exc:
                log.error('mvicuna process call generated an exception: %s' % (exc))

    # merge per-library read lists together
    util.file.concat(per_lb_read_lists, readListAll)
    # remove per-library read lists
    for fl in per_lb_read_lists:
        os.unlink(fl)

    # Filter original input BAM against keep-list
    tools.picard.FilterSamReadsTool().execute(inBam, False, readListAll, outBam, JVMmemory=JVMmemory)
    return 0


def parser_rmdup_mvicuna_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('outBam', help='Output reads, BAM format.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, rmdup_mvicuna_bam, split_args=True)
    return parser


__commands__.append(('rmdup_mvicuna_bam', parser_rmdup_mvicuna_bam))



def parser_rmdup_prinseq_fastq(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq1', help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2', help='Input fastq file; 2nd end of paired-end reads.')
    parser.add_argument('outFastq1', help='Output fastq file; 1st end of paired-end reads.')
    parser.add_argument('outFastq2', help='Output fastq file; 2nd end of paired-end reads.')
    parser.add_argument(
        "--includeUnmated",
        help="Include unmated reads in the main output fastq files (default: %(default)s)",
        default=False,
        action="store_true"
    )
    parser.add_argument('--unpairedOutFastq1', default=None, help='File name of output unpaired reads from 1st end of paired-end reads (independent of --includeUnmated)')
    parser.add_argument('--unpairedOutFastq2', default=None, help='File name of output unpaired reads from 2nd end of paired-end reads (independent of --includeUnmated)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_rmdup_prinseq_fastq)
    return parser


def main_rmdup_prinseq_fastq(args):
    ''' Run prinseq-lite's duplicate removal operation on paired-end
        reads.  Also removes reads with more than one N.
    '''
    prinseq = tools.prinseq.PrinseqTool()
    prinseq.rmdup_fastq_paired(args.inFastq1, args.inFastq2, args.outFastq1, args.outFastq2, args.includeUnmated, args.unpairedOutFastq1, args.unpairedOutFastq2)
    return 0


__commands__.append(('rmdup_prinseq_fastq', parser_rmdup_prinseq_fastq))


def filter_bam_mapped_only(inBam, outBam):
    ''' Samtools to reduce a BAM file to only reads that are
        aligned (-F 4) with a non-zero mapping quality (-q 1)
        and are not marked as a PCR/optical duplicate (-F 1024).
    '''
    tools.samtools.SamtoolsTool().view(['-b', '-q', '1', '-F', '1028'], inBam, outBam)
    tools.picard.BuildBamIndexTool().execute(outBam)
    return 0


def parser_filter_bam_mapped_only(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input aligned reads, BAM format.')
    parser.add_argument('outBam', help='Output sorted indexed reads, filtered to aligned-only, BAM format.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_bam_mapped_only, split_args=True)
    return parser


__commands__.append(('filter_bam_mapped_only', parser_filter_bam_mapped_only))

# ======= Novoalign ========


def parser_novoalign(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('refFasta', help='Reference genome, FASTA format, pre-indexed by Novoindex.')
    parser.add_argument('outBam', help='Output reads, BAM format (aligned).')
    parser.add_argument('--options', default='-r Random', help='Novoalign options (default: %(default)s)')
    parser.add_argument('--min_qual', default=0, help='Filter outBam to minimum mapping quality (default: %(default)s)')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.SortSamTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--NOVOALIGN_LICENSE_PATH',
        default=None,
        dest="novoalign_license_path",
        help='A path to the novoalign.lic file. This overrides the NOVOALIGN_LICENSE_PATH environment variable. (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_novoalign)
    return parser


def main_novoalign(args):
    '''Align reads with Novoalign. Sort and index BAM output.'''
    tools.novoalign.NovoalignTool(license_path=args.novoalign_license_path).execute(
        args.inBam,
        args.refFasta,
        args.outBam,
        options=args.options.split(),
        min_qual=args.min_qual,
        JVMmemory=args.JVMmemory
    )
    return 0


__commands__.append(('novoalign', parser_novoalign))


def parser_novoindex(parser=argparse.ArgumentParser()):
    parser.add_argument('refFasta', help='Reference genome, FASTA format.')
    parser.add_argument(
        '--NOVOALIGN_LICENSE_PATH',
        default=None,
        dest="novoalign_license_path",
        help='A path to the novoalign.lic file. This overrides the NOVOALIGN_LICENSE_PATH environment variable. (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_novoindex)
    return parser

def main_novoindex(args):
    tools.novoalign.NovoalignTool(license_path=args.novoalign_license_path).index_fasta(args.refFasta)
    return 0


__commands__.append(('novoindex', parser_novoindex))

# ========= GATK ==========


def parser_gatk_ug(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('refFasta', help='Reference genome, FASTA format, pre-indexed by Picard.')
    parser.add_argument(
        'outVcf',
        help='''Output calls in VCF format. If this filename ends with .gz,
        GATK will BGZIP compress the output and produce a Tabix index file as well.'''
    )
    parser.add_argument(
        '--options',
        default='--min_base_quality_score 15 -ploidy 4',
        help='UnifiedGenotyper options (default: %(default)s)'
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.gatk.GATKTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--GATK_PATH',
        default=None,
        help='A path containing the GATK jar file. This overrides the GATK_ENV environment variable or the GATK conda package. (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_gatk_ug)
    return parser


def main_gatk_ug(args):
    '''Call genotypes using the GATK UnifiedGenotyper.'''
    tools.gatk.GATKTool(path=args.GATK_PATH).ug(
        args.inBam, args.refFasta,
        args.outVcf, options=args.options.split(),
        JVMmemory=args.JVMmemory
    )
    return 0


__commands__.append(('gatk_ug', parser_gatk_ug))


def parser_gatk_realign(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format, aligned to refFasta.')
    parser.add_argument('refFasta', help='Reference genome, FASTA format, pre-indexed by Picard.')
    parser.add_argument('outBam', help='Realigned reads.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.gatk.GATKTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--GATK_PATH',
        default=None,
        help='A path containing the GATK jar file. This overrides the GATK_ENV environment variable or the GATK conda package. (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_gatk_realign)
    parser.add_argument('--threads', type=int, default=None, help='Number of threads (default: all available cores)')
    return parser


def main_gatk_realign(args):
    '''Local realignment of BAM files with GATK IndelRealigner.'''
    tools.gatk.GATKTool(path=args.GATK_PATH).local_realign(
        args.inBam, args.refFasta,
        args.outBam, JVMmemory=args.JVMmemory,
        threads=args.threads,
    )
    return 0


__commands__.append(('gatk_realign', parser_gatk_realign))

# =========================


def align_and_fix(
    inBam, refFasta,
    outBamAll=None,
    outBamFiltered=None,
    aligner_options='',
    aligner="novoalign",
    bwa_min_score=None,
    novoalign_amplicons_bed=None,
    amplicon_window=4,
    JVMmemory=None,
    threads=None,
    skip_mark_dupes=False,
    gatk_path=None,
    novoalign_license_path=None
):
    ''' Take reads, align to reference with Novoalign, optionally mark duplicates
        with Picard, realign indels with GATK, and optionally filters
        final file to mapped/non-dupe reads.
    '''
    if not (outBamAll or outBamFiltered):
        log.warning("are you sure you meant to do nothing?")
        return

    assert aligner in ["novoalign", "bwa"]

    refFastaCopy = mkstempfname('.ref_copy.fasta')
    shutil.copyfile(refFasta, refFastaCopy)

    tools.picard.CreateSequenceDictionaryTool().execute(refFastaCopy, overwrite=True)
    tools.samtools.SamtoolsTool().faidx(refFastaCopy, overwrite=True)

    if aligner_options is None:
        if aligner=="novoalign":
            aligner_options = '-r Random'
        elif aligner=='bwa':
            aligner_options = '' # use defaults

    bam_aligned = mkstempfname('.aligned.bam')
    if aligner=="novoalign":
        if novoalign_amplicons_bed is not None:
            aligner_options += ' --amplicons {} {}'.format(novoalign_amplicons_bed, amplicon_window)

        tools.novoalign.NovoalignTool(license_path=novoalign_license_path).index_fasta(refFastaCopy)

        tools.novoalign.NovoalignTool(license_path=novoalign_license_path).execute(
            inBam, refFastaCopy, bam_aligned,
            options=aligner_options.split(),
            JVMmemory=JVMmemory
        )
    elif aligner=='bwa':
        bwa = tools.bwa.Bwa()
        bwa.index(refFastaCopy)

        opts = aligner_options.split()

        bwa.align_mem_bam(inBam, refFastaCopy, bam_aligned, min_score_to_filter=bwa_min_score, threads=threads, options=opts)

    if skip_mark_dupes:
        bam_marked = bam_aligned
    else:
        bam_marked = mkstempfname('.mkdup.bam')
        tools.picard.MarkDuplicatesTool().execute(
            [bam_aligned], bam_marked, picardOptions=['CREATE_INDEX=true'],
            JVMmemory=JVMmemory
        )
        os.unlink(bam_aligned)

    tools.samtools.SamtoolsTool().index(bam_marked)

    bam_realigned = mkstempfname('.realigned.bam')
    tools.gatk.GATKTool(path=gatk_path).local_realign(bam_marked, refFastaCopy, bam_realigned, JVMmemory=JVMmemory, threads=threads)
    os.unlink(bam_marked)

    if outBamAll:
        shutil.copyfile(bam_realigned, outBamAll)
        tools.picard.BuildBamIndexTool().execute(outBamAll)
    if outBamFiltered:
        tools.samtools.SamtoolsTool().view(['-b', '-q', '1', '-F', '1028'], bam_realigned, outBamFiltered)
        tools.picard.BuildBamIndexTool().execute(outBamFiltered)
    os.unlink(bam_realigned)


def parser_align_and_fix(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('refFasta', help='Reference genome, FASTA format; will be indexed by Picard and Novoalign.')
    parser.add_argument(
        '--outBamAll',
        default=None,
        help='''Aligned, sorted, and indexed reads.  Unmapped and duplicate reads are
                retained. By default, duplicate reads are marked. If "--skipMarkDupes"
                is specified duplicate reads are included in outout without being marked.'''
    )
    parser.add_argument(
        '--outBamFiltered',
        default=None,
        help='''Aligned, sorted, and indexed reads.  Unmapped reads are removed from this file,
                as well as any marked duplicate reads. Note that if "--skipMarkDupes" is provided,
                duplicates will be not be marked and will be included in the output.'''
    )
    parser.add_argument('--aligner_options', default=None, help='aligner options (default for novoalign: "-r Random", bwa: "-T 30"')
    parser.add_argument('--aligner', choices=['novoalign', 'bwa'], default='novoalign', help='aligner (default: %(default)s)')
    parser.add_argument('--bwa_min_score', type=int, default=None, help='BWA mem on paired reads ignores the -T parameter. Set a value here (e.g. 30) to invoke a custom post-alignment filter (default: no filtration)')
    parser.add_argument('--novoalign_amplicons_bed', default=None, help='Novoalign only: amplicon primer file (BED format) to soft clip')
    parser.add_argument('--amplicon_window', type=int, default=4, help='Novoalign only: amplicon primer window size (default: %(default)s)')
    parser.add_argument('--JVMmemory', default='4g', help='JVM virtual memory size (default: %(default)s)')
    parser.add_argument('--threads', type=int, default=None, help='Number of threads (default: all available cores)')
    parser.add_argument('--skipMarkDupes',
                        help='If specified, duplicate reads will not be marked in the resulting output file.',
                        dest="skip_mark_dupes",
                        action='store_true')
    parser.add_argument(
        '--GATK_PATH',
        default=None,
        dest="gatk_path",
        help='A path containing the GATK jar file. This overrides the GATK_ENV environment variable or the GATK conda package. (default: %(default)s)'
    )
    parser.add_argument(
        '--NOVOALIGN_LICENSE_PATH',
        default=None,
        dest="novoalign_license_path",
        help='A path to the novoalign.lic file. This overrides the NOVOALIGN_LICENSE_PATH environment variable. (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, align_and_fix, split_args=True)
    return parser


__commands__.append(('align_and_fix', parser_align_and_fix))

# =========================


def bwamem_idxstats(inBam, refFasta, outBam=None, outStats=None,
        min_score_to_filter=None, aligner_options=None):
    ''' Take reads, align to reference with BWA-MEM and perform samtools idxstats.
    '''

    assert outBam or outStats, "Either outBam or outStats must be specified"

    if outBam is None:
        bam_aligned = mkstempfname('.aligned.bam')
    else:
        bam_aligned = outBam

    samtools = tools.samtools.SamtoolsTool()
    bwa = tools.bwa.Bwa()

    ref_indexed = util.file.mkstempfname('.reference.fasta')
    shutil.copyfile(refFasta, ref_indexed)
    bwa.index(ref_indexed)

    bwa_opts = [] if aligner_options is None else aligner_options.split()
    bwa.mem(inBam, ref_indexed, bam_aligned, options=bwa_opts,
            min_score_to_filter=min_score_to_filter)

    if outStats is not None:
        samtools.idxstats(bam_aligned, outStats)

    if outBam is None:
        os.unlink(bam_aligned)


def parser_bwamem_idxstats(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('refFasta', help='Reference genome, FASTA format, pre-indexed by Picard and Novoalign.')
    parser.add_argument('--outBam', help='Output aligned, indexed BAM file', default=None)
    parser.add_argument('--outStats', help='Output idxstats file', default=None)
    parser.add_argument(
        '--minScoreToFilter',
        dest="min_score_to_filter",
        type=int,
        help=("Filter bwa alignments using this value as the minimum allowed "
              "alignment score. Specifically, sum the alignment scores across "
              "all alignments for each query (including reads in a pair, "
              "supplementary and secondary alignments) and then only include, "
              "in the output, queries whose summed alignment score is at least "
              "this value. This is only applied when the aligner is 'bwa'. "
              "The filtering on a summed alignment score is sensible for reads "
              "in a pair and supplementary alignments, but may not be "
              "reasonable if bwa outputs secondary alignments (i.e., if '-a' "
              "is in the aligner options). (default: not set - i.e., do not "
              "filter bwa's output)")
    )
    parser.add_argument(
        '--alignerOptions',
        dest="aligner_options",
        help="bwa options (default: bwa defaults)")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, bwamem_idxstats, split_args=True)
    return parser

__commands__.append(('bwamem_idxstats', parser_bwamem_idxstats))


def parser_extract_tarball(parser=argparse.ArgumentParser()):
    parser.add_argument('tarfile', help='Input tar file. May be "-" for stdin.')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('--compression',
        help='Compression type (default: %(default)s). Auto-detect is incompatible with stdin input unless pipe_hint is specified.',
        choices=('gz', 'bz2', 'lz4', 'zip', 'none', 'auto'),
        default='auto')
    parser.add_argument('--pipe_hint',
        help='If tarfile is stdin, you can provide a file-like URI string for pipe_hint which ends with a common compression file extension if you want to use compression=auto.',
        default=None)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_extract_tarball, split_args=True)
    return parser
def main_extract_tarball(*args, **kwargs):
    ''' Extract an input .tar, .tgz, .tar.gz, .tar.bz2, .tar.lz4, or .zip file
        to a given directory (or we will choose one on our own). Emit the
        resulting directory path to stdout.
    '''
    print(util.file.extract_tarball(*args, **kwargs))
__commands__.append(('extract_tarball', parser_extract_tarball))

# =========================

def fasta_read_names(in_fasta, out_read_names):
    """Save the read names of reads in a .fasta file to a text file"""
    with util.file.open_or_gzopen(in_fasta) as in_fasta_f, open(out_read_names, 'wt') as out_read_names_f:
        last_read_name = None
        for line in in_fasta_f:
            if line.startswith('>'):
                read_name = line[1:].strip()
                if read_name.endswith('/1') or read_name.endswith('/2'):
                    read_name = read_name[:-2]
                if read_name != last_read_name:
                    out_read_names_f.write(read_name+'\n')
                last_read_name = read_name


def read_names(in_reads, out_read_names, threads=None):
    """Extract read names from a sequence file"""
    _in_reads = in_reads
    with util.file.tmp_dir(suffix='_read_names.txt') as t_dir:
        if in_reads.endswith('.bam'):
            _in_reads = os.path.join(t_dir, 'reads.fasta')
            tools.samtools.SamtoolsTool().bam2fa(in_reads, _in_reads)
        fasta_read_names(_in_reads, out_read_names)

def parser_read_names(parser=argparse.ArgumentParser()):
    parser.add_argument('in_reads', help='the input reads ([compressed] fasta or bam)')
    parser.add_argument('out_read_names', help='the read names')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, read_names, split_args=True)
    return parser

__commands__.append(('read_names', parser_read_names))


# =========================

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
