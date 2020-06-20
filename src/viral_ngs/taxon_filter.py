#!/usr/bin/env python3
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''

from __future__ import print_function
__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
    + "hlevitin@broadinstitute.org"
__commands__ = []

import argparse
import glob
import logging
import subprocess
import os
import math
import tempfile
import shutil
import concurrent.futures
import contextlib

from Bio import SeqIO
import pysam

import util.cmd
import util.file
import util.misc
import tools
import tools.prinseq
import tools.picard
import tools.samtools
from util.file import mkstempfname
from errors import QCError

import classify.blast
import classify.last
import classify.bmtagger
import read_utils

log = logging.getLogger(__name__)


# =======================
# ***  deplete  ***
# =======================

def parser_deplete(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument('revertBam', nargs='?', help='Output BAM: read markup reverted with Picard.')
    parser.add_argument('bwaBam', help='Output BAM: depleted of reads with BWA.')
    parser.add_argument('bmtaggerBam', help='Output BAM: depleted of reads with BMTagger.')
    parser.add_argument(
        'blastnBam', help='Output BAM: bmtaggerBam run through another depletion of reads with BLASTN.'
    )
    parser.add_argument(
        '--bwaDbs',
        nargs='*',
        default=(),
        help='Reference databases for blast to deplete from input.'
    )
    parser.add_argument(
        '--bmtaggerDbs',
        nargs='*',
        default=(),
        help='''Reference databases to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.'''
    )
    parser.add_argument(
        '--blastDbs',
        nargs='*',
        default=(),
        help='Reference databases for blast to deplete from input.'
    )
    parser.add_argument('--srprismMemory', dest="srprism_memory", type=int, default=7168, help='Memory for srprism.')
    parser.add_argument("--chunkSize", type=int, default=1000000, help='blastn chunk size (default: %(default)s)')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size for Picard FilterSamReads (default: %(default)s)'
    )
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete)

    return parser

def main_deplete(args):
    ''' Run the entire depletion pipeline: bwa, bmtagger, blastn.
    '''

    assert len(args.bmtaggerDbs) + len(args.blastDbs) + len(args.bwaDbs) > 0

    # only RevertSam if inBam is already aligned
    # Most of the time the input will be unaligned
    # so we can save save time if we can skip RevertSam in the unaligned case
    #
    # via the SAM/BAM spec, if the file is aligned, an SQ line should be present
    # in the header. Using pysam, we can check this if header['SQ'])>0
    #   https://samtools.github.io/hts-specs/SAMv1.pdf

    # if the user has requested a revertBam

    with read_utils.revert_bam_if_aligned(              args.inBam,
                                        revert_bam    = args.revertBam,
                                        clear_tags    = args.clear_tags,
                                        tags_to_clear = args.tags_to_clear,
                                        picardOptions = ['MAX_DISCARD_FRACTION=0.5'],
                                        JVMmemory     = args.JVMmemory,
                                        sanitize      = not args.do_not_sanitize) as bamToDeplete:
        multi_db_deplete_bam(
            bamToDeplete,
            args.bwaDbs,
            deplete_bwa_bam,
            args.bwaBam,
            threads=args.threads
        )

    def bmtagger_wrapper(inBam, db, outBam, JVMmemory=None):
        return deplete_bmtagger_bam(inBam, db, outBam, srprism_memory=args.srprism_memory, JVMmemory=JVMmemory)

    multi_db_deplete_bam(
        args.bwaBam,
        args.bmtaggerDbs,
        bmtagger_wrapper,
        args.bmtaggerBam,
        JVMmemory=args.JVMmemory
    )

    # if the user has not specified saving a revertBam, we used a temp file and can remove it
    if not args.revertBam:
        os.unlink(revertBamOut)
    else:
        if os.path.getsize(args.revertBam) == 0:
            with util.file.tempfname('.empty.sam') as empty_sam:
                samtools.dumpHeader(args.inBam, empty_sam)
                samtools.view(['-b'], empty_sam, args.revertBam)

    multi_db_deplete_bam(
        args.bmtaggerBam,
        args.blastDbs,
        deplete_blastn_bam,
        args.blastnBam,
        chunkSize=args.chunkSize,
        threads=args.threads,
        JVMmemory=args.JVMmemory
    )
    return 0

__commands__.append(('deplete', parser_deplete))


# =======================
# ***  filter_lastal  ***
# =======================


def filter_lastal_bam(
    inBam,
    db,
    outBam,
    max_gapless_alignments_per_position=1,
    min_length_for_initial_matches=5,
    max_length_for_initial_matches=50,
    max_initial_matches_per_position=100,
    error_on_reads_in_neg_control=False,
    neg_control_prefixes=None, #set below: "neg","water","NTC"
    negative_control_reads_threshold=0,
    JVMmemory=None, threads=None
):
    ''' Restrict input reads to those that align to the given
        reference database using LASTAL.
    '''
    neg_control_prefixes = neg_control_prefixes or ["neg","water","NTC"]

    with util.file.tmp_dir('-lastdb') as tmp_db_dir:
        # index db if necessary
        lastdb = classify.last.Lastdb()
        if not lastdb.is_indexed(db):
            db = lastdb.build_database(db, os.path.join(tmp_db_dir, 'lastdb'))

        with util.file.tempfname('.read_ids.txt') as hitList:
            number_of_hits=0

            # look for lastal hits in BAM and write to temp file
            with open(hitList, 'wt') as outf:
                for read_id in classify.last.Lastal().get_hits(
                        inBam, db,
                        max_gapless_alignments_per_position,
                        min_length_for_initial_matches,
                        max_length_for_initial_matches,
                        max_initial_matches_per_position,
                        threads=threads
                    ):
                    number_of_hits+=1
                    outf.write(read_id + '\n')

            if error_on_reads_in_neg_control:
                sample_name=os.path.basename(inBam)
                if any(sample_name.lower().startswith(prefix.lower()) for prefix in neg_control_prefixes):
                    if number_of_hits > max(0,negative_control_reads_threshold):
                        log.warning("Error raised due to reads in negative control; re-run this without '--errorOnReadsInNegControl' if this execution should succeed.")
                        raise QCError("The sample '{}' appears to be a negative control, but it contains {} reads after filtering to desired taxa.".format(sample_name,number_of_hits))

            # filter original BAM file against keep list
            tools.picard.FilterSamReadsTool().execute(inBam, False, hitList, outBam, JVMmemory=JVMmemory)


def parser_filter_lastal_bam(parser=argparse.ArgumentParser()):
    parser.add_argument("inBam", help="Input reads")
    parser.add_argument("db", help="Database of taxa we keep")
    parser.add_argument("outBam", help="Output reads, filtered to refDb")
    parser.add_argument(
        '-n',
        dest="max_gapless_alignments_per_position",
        help='maximum gapless alignments per query position (default: %(default)s)',
        type=int,
        default=1
    )
    parser.add_argument(
        '-l',
        dest="min_length_for_initial_matches",
        help='minimum length for initial matches (default: %(default)s)',
        type=int,
        default=5
    )
    parser.add_argument(
        '-L',
        dest="max_length_for_initial_matches",
        help='maximum length for initial matches (default: %(default)s)',
        type=int,
        default=50
    )
    parser.add_argument(
        '-m',
        dest="max_initial_matches_per_position",
        help='maximum initial matches per query position (default: %(default)s)',
        type=int,
        default=100
    )
    parser.add_argument(
        '--errorOnReadsInNegControl',
        dest="error_on_reads_in_neg_control",
        help='If specified, the function will return an error if there are reads after filtering for samples with names containing: (water,neg,ntc) (default: %(default)s)',
        action="store_true",
    )
    parser.add_argument(
        '--negativeControlReadsThreshold',
        dest="negative_control_reads_threshold",
        help='maximum number of reads (single-end) or read pairs (paired-end) to tolerate in samples identified as negative controls (default: %(default)s)',
        type=int,
        default=0
    )
    parser.add_argument(
        '--negControlPrefixes',
        dest="neg_control_prefixes",
        default=["neg","water","NTC"],
        nargs='*',
        help='Bam file name prefixes to interpret as negative controls, space-separated (default: %(default)s)'
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_lastal_bam, split_args=True)
    return parser


__commands__.append(('filter_lastal_bam', parser_filter_lastal_bam))


# ==============================
# ***  deplete_bmtagger_bam  ***
# ==============================


def deplete_bmtagger_bam(inBam, db, outBam, srprism_memory=7168, JVMmemory=None):
    """
    Use bmtagger to partition the input reads into ones that match at least one
        of the databases and ones that don't match any of the databases.
    inBam: paired-end input reads in BAM format.
    db: bmtagger expects files
        db.bitmask created by bmtool, and
        db.srprism.idx, db.srprism.map, etc. created by srprism mkindex
    outBam: the output BAM files to hold the unmatched reads.
    srprism_memory: srprism memory in megabytes.
    """
    bmtaggerPath = classify.bmtagger.BmtaggerShTool().install_and_get_path()

    # bmtagger calls several executables in the same directory, and blastn;
    # make sure they are accessible through $PATH
    blastnPath = classify.blast.BlastnTool().install_and_get_path()
    path = os.environ['PATH'].split(os.pathsep)
    for t in (bmtaggerPath, blastnPath):
        d = os.path.dirname(t)
        if d not in path:
            path = [d] + path
    path = os.pathsep.join(path)
    os.environ['PATH'] = path

    with util.file.tempfname('.1.fastq') as inReads1:
        tools.samtools.SamtoolsTool().bam2fq(inBam, inReads1)

        with util.file.tempfname('.bmtagger.conf') as bmtaggerConf:
            with open(bmtaggerConf, 'w') as f:
                # Default srprismopts: "-b 100000000 -n 5 -R 0 -r 1 -M 7168"
                print('srprismopts="-b 100000000 -n 5 -R 0 -r 1 -M {srprism_memory} --paired false"'.format(srprism_memory=srprism_memory), file=f)

            with extract_build_or_use_database(db, bmtagger_build_db, 'bitmask', tmp_suffix="-bmtagger", db_prefix="bmtagger") as (db_prefix,tempDir):
                matchesFile = mkstempfname('.txt')
                cmdline = [
                    bmtaggerPath, '-b', db_prefix + '.bitmask', '-C', bmtaggerConf, '-x', db_prefix + '.srprism', '-T', tempDir, '-q1',
                    '-1', inReads1, '-o', matchesFile
                ]
                log.debug(' '.join(cmdline))
                util.misc.run_and_print(cmdline, check=True)

    tools.picard.FilterSamReadsTool().execute(inBam, True, matchesFile, outBam, JVMmemory=JVMmemory)

def parser_deplete_bam_bmtagger(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument(
        'refDbs',
        nargs='+',
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.'''
    )
    parser.add_argument('outBam', help='Output BAM file.')
    parser.add_argument('--srprismMemory', dest="srprism_memory", type=int, default=7168, help='Memory for srprism.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_bam_bmtagger)
    return parser

def main_deplete_bam_bmtagger(args):
    '''Use bmtagger to deplete input reads against several databases.'''

    def bmtagger_wrapper(inBam, db, outBam, JVMmemory=None):
        return deplete_bmtagger_bam(inBam, db, outBam, srprism_memory=args.srprism_memory, JVMmemory=JVMmemory)

    with read_utils.revert_bam_if_aligned(              args.inBam,
                                        clear_tags    = args.clear_tags,
                                        tags_to_clear = args.tags_to_clear,
                                        picardOptions = ['MAX_DISCARD_FRACTION=0.5'],
                                        JVMmemory     = args.JVMmemory,
                                        sanitize      = not args.do_not_sanitize) as bamToDeplete:
        multi_db_deplete_bam(
            args.inBam,
            args.refDbs,
            bmtagger_wrapper,
            args.outBam,
            JVMmemory=args.JVMmemory
        )

__commands__.append(('deplete_bam_bmtagger', parser_deplete_bam_bmtagger))


def multi_db_deplete_bam(inBam, refDbs, deplete_method, outBam, **kwargs):

    tmpDb = None
    if len(refDbs)>1 and not any(
            not os.path.exists(db)  # indexed db prefix
            or os.path.isdir(db)       # indexed db in directory
            or (os.path.isfile(db) and ('.tar' in db or '.tgz' in db or '.zip' in db)) # packaged indexed db
            for db in refDbs):
        # this is a scenario where all refDbs are unbuilt fasta
        # files. we can simplify and speed up execution by
        # concatenating them all and running deplete_method
        # just once
        tmpDb = mkstempfname('.fasta')
        merge_compressed_files(refDbs, tmpDb, sep='\n')
        refDbs = [tmpDb]

    samtools = tools.samtools.SamtoolsTool()
    tmpBamIn = inBam
    for db in refDbs:
        if not samtools.isEmpty(tmpBamIn):
            tmpBamOut = mkstempfname('.bam')
            deplete_method(tmpBamIn, db, tmpBamOut, **kwargs)
            if tmpBamIn != inBam:
                os.unlink(tmpBamIn)
            tmpBamIn = tmpBamOut
    shutil.copyfile(tmpBamIn, outBam)

    if tmpDb:
        os.unlink(tmpDb)


# ========================
# ***  deplete_blastn  ***
# ========================


def _run_blastn_chunk(db, input_fasta, out_hits, blast_threads):
    """ run blastn on the input fasta file. this is intended to be run in parallel
        by blastn_chunked_fasta
    """
    with util.file.open_or_gzopen(out_hits, 'wt') as outf:
        for read_id in classify.blast.BlastnTool().get_hits_fasta(input_fasta, db, threads=blast_threads):
            outf.write(read_id + '\n')

def blastn_chunked_fasta(fasta, db, out_hits, chunkSize=1000000, threads=None):
    """
    Helper function: blastn a fasta file, overcoming apparent memory leaks on
    an input with many query sequences, by splitting it into multiple chunks
    and running a new blastn process on each chunk. Return a list of output
    filenames containing hits
    """
    # the lower bound of how small a fasta chunk can be.
    # too small and the overhead of spawning a new blast process
    # will be detrimental relative to actual computation time
    MIN_CHUNK_SIZE = 20000

    # just in case blast is not installed, install it once, not many times in parallel!
    classify.blast.BlastnTool().install()

    # clamp threadcount to number of CPU cores
    threads = util.misc.sanitize_thread_count(threads)

    # determine size of input data; records in fasta file
    number_of_reads = util.file.fasta_length(fasta)
    log.debug("number of reads in fasta file %s" % number_of_reads)
    if number_of_reads == 0:
        util.file.make_empty(out_hits)

    # divide (max, single-thread) chunksize by thread count
    # to find the  absolute max chunk size per thread
    chunk_max_size_per_thread = chunkSize // threads

    # find the chunk size if evenly divided among blast threads
    reads_per_thread = number_of_reads // threads

    # use the smaller of the two chunk sizes so we can run more copies of blast in parallel
    chunkSize = min(reads_per_thread, chunk_max_size_per_thread)

    # if the chunk size is too small, impose a sensible size
    chunkSize = max(chunkSize, MIN_CHUNK_SIZE)

    log.debug("chunk_max_size_per_thread %s" % chunk_max_size_per_thread)

    # adjust chunk size so we don't have a small fraction
    # of a chunk running in its own blast process
    # if the size of the last chunk is <80% the size of the others,
    # decrease the chunk size until the last chunk is 80%
    # this is bounded by the MIN_CHUNK_SIZE
    while (number_of_reads / chunkSize) % 1 < 0.8 and chunkSize > MIN_CHUNK_SIZE:
        chunkSize = chunkSize - 1

    log.debug("blastn chunk size %s" % chunkSize)
    log.debug("number of chunks to create %s" % (number_of_reads / chunkSize))
    log.debug("blastn parallel instances %s" % threads)

    # chunk the input file. This is a sequential operation
    input_fastas = []
    with open(fasta, "rt") as fastaFile:
        record_iter = SeqIO.parse(fastaFile, "fasta")
        for batch in util.misc.batch_iterator(record_iter, chunkSize):
            chunk_fasta = mkstempfname('.fasta')

            with open(chunk_fasta, "wt") as handle:
                SeqIO.write(batch, handle, "fasta")
            batch = None
            input_fastas.append(chunk_fasta)

    num_chunks = len(input_fastas)
    log.debug("number of chunk files to be processed by blastn %d" % num_chunks)

    # run blastn on each of the fasta input chunks
    hits_files = list(mkstempfname('.hits.txt') for f in input_fastas)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # If we have so few chunks that there are cpus left over,
        # divide extra cpus evenly among chunks where possible
        # rounding to 1 if there are more chunks than extra threads.
        # Then double up this number to better maximize CPU usage.
        cpus_leftover = threads - num_chunks
        blast_threads = 2*max(1, int(cpus_leftover / num_chunks))
        for i in range(num_chunks):
            executor.submit(
                _run_blastn_chunk, db, input_fastas[i], hits_files[i], blast_threads)

    # merge results and clean up
    util.file.cat(out_hits, hits_files)
    for i in range(num_chunks):
        os.unlink(input_fastas[i])
        os.unlink(hits_files[i])


def deplete_blastn_bam(inBam, db, outBam, threads=None, chunkSize=1000000, JVMmemory=None):
#def deplete_blastn_bam(inBam, db, outBam, threads, chunkSize=0, JVMmemory=None):
    'Use blastn to remove reads that match at least one of the databases.'

    blast_hits = mkstempfname('.blast_hits.txt')

    with extract_build_or_use_database(db, blastn_build_db, 'nin', tmp_suffix="-blastn_db_unpack", db_prefix="blastn") as (db_prefix,tempDir):
        if chunkSize:
            ## chunk up input and perform blastn in several parallel threads
            with util.file.tempfname('.fasta') as reads_fasta:
                tools.samtools.SamtoolsTool().bam2fa(inBam, reads_fasta)
                log.info("running blastn on %s against %s", inBam, db)
                blastn_chunked_fasta(reads_fasta, db_prefix, blast_hits, chunkSize, threads)

        else:
            ## pipe tools together and run blastn multithreaded
            with open(blast_hits, 'wt') as outf:
                for read_id in classify.blast.BlastnTool().get_hits_bam(inBam, db_prefix, threads=threads):
                    outf.write(read_id + '\n')

    # Deplete BAM of hits
    tools.picard.FilterSamReadsTool().execute(inBam, True, blast_hits, outBam, JVMmemory=JVMmemory)
    os.unlink(blast_hits)


def parser_deplete_blastn_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument('refDbs', nargs='+', help='One or more reference databases for blast. '
                         'An ephemeral database will be created if a fasta file is provided.')
    parser.add_argument('outBam', help='Output BAM file with matching reads removed.')
    parser.add_argument("--chunkSize", type=int, default=1000000, help='FASTA chunk size (default: %(default)s)')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_blastn_bam)
    return parser


def main_deplete_blastn_bam(args):
    '''Use blastn to remove reads that match at least one of the specified databases.'''

    def wrapper(inBam, db, outBam, threads, JVMmemory=None):
        return deplete_blastn_bam(inBam, db, outBam, threads=threads, chunkSize=args.chunkSize, JVMmemory=JVMmemory)

    with read_utils.revert_bam_if_aligned(              args.inBam,
                                        clear_tags    = args.clear_tags,
                                        tags_to_clear = args.tags_to_clear,
                                        picardOptions = ['MAX_DISCARD_FRACTION=0.5'],
                                        JVMmemory     = args.JVMmemory,
                                        sanitize      = not args.do_not_sanitize) as bamToDeplete:
        multi_db_deplete_bam(bamToDeplete, args.refDbs, wrapper, args.outBam, threads=args.threads, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('deplete_blastn_bam', parser_deplete_blastn_bam))


@contextlib.contextmanager
def extract_build_or_use_database(db, db_build_command, db_extension_to_expect, tmp_suffix='db_unpack', db_prefix="db"):
    '''
    db_extension_to_expect = file extension, sans dot prefix
    '''
    with util.file.tmp_dir(tmp_suffix) as tempDbDir:
        db_dir = ""
        if os.path.exists(db):
            if os.path.isfile(db):
                # this is a single file
                if db.endswith('.fasta') or db.endswith('.fa') or ('.fasta.' in db) or ('.fa.' in db):
                    # this is an unindexed fasta file, we will need to index it
                    # function should conform to the signature:
                    # db_build_command(inputFasta, outputDirectory, outputFilePrefix)
                    # the function will need to be able to handle lz4, etc.
                    db_build_command(db, tempDbDir, db_prefix)
                    db_dir = tempDbDir
                else:
                    # this is a tarball with prebuilt indexes
                    db_dir = util.file.extract_tarball(db, tempDbDir)
            else:
                # this is a directory
                db_dir = db
            # this directory should have a .{ext} file, where {ext} is specific to the type of db
            hits = list(glob.glob(os.path.join(db_dir, '*.{ext}'.format(ext=db_extension_to_expect))))
            if len(hits) == 0:
                raise Exception("The blast database does not appear to a *.{ext} file.".format(ext=db_extension_to_expect))
            elif len(hits) == 1:
                db_prefix = hits[0][:-(len('.{ext}'.format(ext=db_extension_to_expect)))]  # remove the '.extension'
            elif len(hits) >1:
                db_prefix = os.path.commonprefix(hits).rsplit('.', 1)[0] # remove extension and split-db prefix
        else:
            # this is simply a prefix to a bunch of files, not an actual file
            db_prefix = db.rsplit('.', 1)[0] if db.endswith('.') else db

        yield (db_prefix,tempDbDir)

# ========================
# ***  deplete_bwa  ***
# ========================

def deplete_bwa_bam(inBam, db, outBam, threads=None, clear_tags=True, tags_to_clear=None, JVMmemory=None):
    'Use bwa to remove reads from an unaligned bam that match at least one of the databases.'
    tags_to_clear = tags_to_clear or []

    threads = util.misc.sanitize_thread_count(threads)

    with extract_build_or_use_database(db, bwa_build_db, 'bwt', tmp_suffix="-bwa_db_unpack", db_prefix="bwa") as (db_prefix,tempDbDir):
        with util.file.tempfname('.aligned.sam') as aligned_sam:
            tools.bwa.Bwa().align_mem_bam(inBam, db_prefix, aligned_sam, threads=threads, should_index=False, JVMmemory=JVMmemory)
        #with util.file.fifo(name='filtered.sam') as filtered_sam:
            with util.file.tempfname('.filtered.sam') as filtered_sam:
                # filter proper pairs
                tools.samtools.SamtoolsTool().view(['-h','-F0x2'], aligned_sam, filtered_sam)

                picardOptions = []
                if clear_tags:
                    for tag in tags_to_clear:
                        picardOptions.append("ATTRIBUTE_TO_CLEAR={}".format(tag))
                tools.picard.RevertSamTool().execute(
                   filtered_sam,
                   outBam,
                   picardOptions=['SORT_ORDER=queryname'] + picardOptions,
                    JVMmemory=JVMmemory
                )
            # TODO: consider using Bwa().mem() so the input bam is not broken out by read group
            # TODO: pipe bwa input directly to samtools process (need to use Bwa().mem() directly, )
            #       with Popen to background bwa process

def parser_deplete_bwa_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument('refDbs', nargs='+', help='One or more reference databases for bwa. '
                         'An ephemeral database will be created if a fasta file is provided.')
    parser.add_argument('outBam', help='Ouput BAM file with matching reads removed.')
    parser = read_utils.parser_revert_sam_common(parser)
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_bwa_bam)
    return parser

def main_deplete_bwa_bam(args):
    '''Use BWA to remove reads that match at least one of the specified databases.'''
    with read_utils.revert_bam_if_aligned(              args.inBam,
                                        clear_tags    = args.clear_tags,
                                        tags_to_clear = args.tags_to_clear,
                                        picardOptions = ['MAX_DISCARD_FRACTION=0.5'],
                                        JVMmemory     = args.JVMmemory,
                                        sanitize      = not args.do_not_sanitize) as bamToDeplete:

        #def wrapper(inBam, db, outBam, threads, JVMmemory=None):
        #    return deplete_bwa_bam(inBam, db, outBam, threads=threads, )
        multi_db_deplete_bam(bamToDeplete, args.refDbs, deplete_bwa_bam, args.outBam, threads=args.threads, clear_tags=args.clear_tags, tags_to_clear=args.tags_to_clear, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('deplete_bwa_bam', parser_deplete_bwa_bam))


# ========================
# ***  lastal_build_db  ***
# ========================


def lastal_build_db(inputFasta, outputDirectory, outputFilePrefix):
    ''' build a database for use with last based on an input fasta file '''
    if outputFilePrefix:
        outPrefix = outputFilePrefix
    else:
        baseName = os.path.basename(inputFasta)
        fileNameSansExtension = os.path.splitext(baseName)[0]
        outPrefix = fileNameSansExtension

    classify.last.Lastdb().build_database(inputFasta, os.path.join(outputDirectory, outPrefix))


def parser_lastal_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('inputFasta', help='Location of the input FASTA file')
    parser.add_argument('outputDirectory', help='Location for the output files (default is cwd: %(default)s)')
    parser.add_argument(
        '--outputFilePrefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, lastal_build_db, split_args=True)
    return parser


__commands__.append(('lastal_build_db', parser_lastal_build_db))

# ================================
# ***  merge_compressed_files  ***
# ================================

def merge_compressed_files(inFiles, outFile, sep=''):
    ''' Take a collection of input text files, possibly compressed,
        and concatenate into a single output text file.
    '''
    with util.file.open_or_gzopen(outFile, 'wt') as outf:
        first = True
        for infname in inFiles:
            if not first:
                if sep:
                    outf.write(sep)
            else:
                first = False
            with util.file.open_or_gzopen(infname, 'rt', newline=None) as inf:
                shutil.copyfileobj(inf, outf)

# ========================
# ***  bwa_build_db  ***
# ========================


def bwa_build_db(inputFasta, outputDirectory, outputFilePrefix):
    """ Create a database for use with bwa from an input reference FASTA file
    """

    new_fasta = None
    if not (inputFasta.endswith('.fasta') or inputFasta.endswith('.fa')):
        new_fasta = util.file.mkstempfname('.fasta')
        with util.file.open_or_gzopen(inputFasta, 'rt', newline=None) as inf, open(new_fasta, 'wt') as outf:
            shutil.copyfileobj(inf, outf)
        inputFasta = new_fasta

    # make the output path if it does not exist
    util.file.mkdir_p(outputDirectory)

    if outputFilePrefix:
        outPrefix = outputFilePrefix
    else:
        baseName = os.path.basename(inputFasta)
        fileNameSansExtension = os.path.splitext(baseName)[0]
        outPrefix = fileNameSansExtension

    tools.bwa.Bwa().index(inputFasta, output=os.path.join(outputDirectory, outPrefix))

    if new_fasta is not None:
        os.unlink(new_fasta)


def parser_bwa_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('inputFasta', help='Location of the input FASTA file')
    parser.add_argument('outputDirectory', help='Location for the output files')
    parser.add_argument(
        '--outputFilePrefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, bwa_build_db, split_args=True)
    return parser


__commands__.append(('bwa_build_db', parser_bwa_build_db))


# ========================
# ***  blastn_build_db  ***
# ========================


def blastn_build_db(inputFasta, outputDirectory, outputFilePrefix):
    """ Create a database for use with blastn from an input reference FASTA file
    """

    new_fasta = None
    if not (inputFasta.endswith('.fasta') or inputFasta.endswith('.fa')):
        new_fasta = util.file.mkstempfname('.fasta')
        with util.file.open_or_gzopen(inputFasta, 'rt', newline=None) as inf, open(new_fasta, 'wt') as outf:
            shutil.copyfileobj(inf, outf)
        inputFasta = new_fasta

    if outputFilePrefix:
        outPrefix = outputFilePrefix
    else:
        baseName = os.path.basename(inputFasta)
        fileNameSansExtension = os.path.splitext(baseName)[0]
        outPrefix = fileNameSansExtension

    blastdb_path = classify.blast.MakeblastdbTool().build_database(inputFasta, os.path.join(outputDirectory, outPrefix))

    if new_fasta is not None:
        os.unlink(new_fasta)


def parser_blastn_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('inputFasta', help='Location of the input FASTA file')
    parser.add_argument('outputDirectory', help='Location for the output files')
    parser.add_argument(
        '--outputFilePrefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, blastn_build_db, split_args=True)
    return parser


__commands__.append(('blastn_build_db', parser_blastn_build_db))

# ========================
# ***  bmtagger_build_db  ***
# ========================


def bmtagger_build_db(inputFasta, outputDirectory, outputFilePrefix, word_size=18):
    """ Create a database for use with Bmtagger from an input FASTA file.
    """

    new_fasta = None
    if not (inputFasta.endswith('.fasta') or inputFasta.endswith('.fa')):
        new_fasta = util.file.mkstempfname('.fasta')
        with util.file.open_or_gzopen(inputFasta, 'rt', newline=None) as inf, open(new_fasta, 'wt') as outf:
            shutil.copyfileobj(inf, outf)
        inputFasta = new_fasta

    if outputFilePrefix:
        outPrefix = outputFilePrefix
    else:
        baseName = os.path.basename(inputFasta)
        fileNameSansExtension = os.path.splitext(baseName)[0]
        outPrefix = fileNameSansExtension

    log.debug("building bmtagger and srprism databases on {}".format(os.path.join(outputDirectory, outPrefix)))
    bmtooldb_path = classify.bmtagger.BmtoolTool().build_database(
        inputFasta, os.path.join(outputDirectory, outPrefix + ".bitmask"), word_size=word_size
    )
    srprismdb_path = classify.bmtagger.SrprismTool().build_database(
        inputFasta, os.path.join(outputDirectory, outPrefix + ".srprism")
    )

    if new_fasta is not None:
        os.unlink(new_fasta)


def parser_bmtagger_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('inputFasta', help='Location of the input FASTA file')
    parser.add_argument(
        'outputDirectory',
        help='Location for the output files (Where *.bitmask and *.srprism files will be stored)'
    )
    parser.add_argument(
        '--outputFilePrefix',
        help='Prefix for the output file name (default: inputFasta name, sans ".fasta" extension)'
    )
    parser.add_argument(
        '--word_size',
        type=int,
        default=18,
        help='Database word size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, bmtagger_build_db, split_args=True)
    return parser


__commands__.append(('bmtagger_build_db', parser_bmtagger_build_db))

# ========================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
