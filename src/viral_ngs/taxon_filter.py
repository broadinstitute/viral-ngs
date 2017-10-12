#!/usr/bin/env python
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''

from __future__ import print_function
__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
    + "hlevitin@broadinstitute.org"
__commands__ = []

import argparse
import logging
import subprocess
import os
import math
import tempfile
import shutil
import functools
import concurrent.futures

from Bio import SeqIO
import pysam

import util.cmd
import util.file
import util.misc
import tools
import tools.blast
import tools.last
import tools.prinseq
import tools.bmtagger
import tools.picard
import tools.samtools
from util.file import mkstempfname
import read_utils

log = logging.getLogger(__name__)

# =======================
# ***  deplete_human  ***
# =======================


def parser_deplete_human(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument('revertBam', nargs='?', help='Output BAM: read markup reverted with Picard.')
    parser.add_argument('bmtaggerBam', help='Output BAM: depleted of human reads with BMTagger.')
    parser.add_argument('rmdupBam', help='Output BAM: bmtaggerBam run through M-Vicuna duplicate removal.')
    parser.add_argument(
        'blastnBam', help='Output BAM: rmdupBam run through another depletion of human reads with BLASTN.'
    )
    parser.add_argument(
        '--taxfiltBam',
        help='Output BAM: blastnBam run through taxonomic selection via LASTAL.',
        default=None
    )
    parser.add_argument(
        '--bmtaggerDbs',
        nargs='+',
        required=True,
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.'''
    )
    parser.add_argument(
        '--blastDbs',
        nargs='+',
        required=True,
        help='One or more reference databases for blast to deplete from input.'
    )
    parser.add_argument(
        '--lastDb',
        help='One reference database for last (required if --taxfiltBam is specified).',
        default=None
    )
    parser.add_argument('--threads', type=int, default=4, help='The number of threads to use in running blastn.')
    parser.add_argument('--srprismMemory', dest="srprism_memory", type=int, default=7168, help='Memory for srprism.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size for Picard FilterSamReads (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_human)
    return parser


def main_deplete_human(args):
    ''' Run the entire depletion pipeline: bmtagger, mvicuna, blastn.
        Optionally, use lastal to select a specific taxon of interest.'''

    # only RevertSam if inBam is already aligned
    # Most of the time the input will be unaligned
    # so we can save save time if we can skip RevertSam in the unaligned case
    #
    # via the SAM/BAM spec, if the file is aligned, an SQ line should be present
    # in the header. Using pysam, we can check this if header['SQ'])>0
    #   https://samtools.github.io/hts-specs/SAMv1.pdf

    # if the user has requested a revertBam
    revertBamOut = args.revertBam if args.revertBam else mkstempfname('.bam')

    bamToDeplete = args.inBam
    with pysam.AlignmentFile(args.inBam, 'rb', check_sq=False) as bam:
        # if it looks like the bam is aligned, revert it
        if 'SQ' in bam.header and len(bam.header['SQ'])>0:
            tools.picard.RevertSamTool().execute(
                args.inBam, revertBamOut, picardOptions=['SORT_ORDER=queryname', 'SANITIZE=true']
            )
            bamToDeplete = revertBamOut
        else:
            # if we don't need to produce a revertBam file
            # but the user has specified one anyway
            # simply touch the output
            if args.revertBam:
                log.warning("An output was specified for 'revertBam', but the input is unaligned, so RevertSam was not needed. Touching the output.")
                util.file.touch(revertBamOut)
                # TODO: error out? run RevertSam anyway?

    def bmtagger_wrapper(inBam, db, outBam, threads, JVMmemory=None):
        return deplete_bmtagger_bam(inBam, db, outBam, threads=threads, srprism_memory=args.srprism_memory, JVMmemory=JVMmemory)

    multi_db_deplete_bam(
        bamToDeplete,
        args.bmtaggerDbs,
        bmtagger_wrapper,
        args.bmtaggerBam,
        threads=args.threads,
        JVMmemory=args.JVMmemory
    )

    # if the user has not specified saving a revertBam, we used a temp file and can remove it
    if not args.revertBam:
        os.unlink(revertBamOut)
    read_utils.rmdup_mvicuna_bam(args.bmtaggerBam, args.rmdupBam, JVMmemory=args.JVMmemory)

    multi_db_deplete_bam(
        args.rmdupBam,
        args.blastDbs,
        deplete_blastn_bam,
        args.blastnBam,
        threads=args.threads,
        JVMmemory=args.JVMmemory
    )
    if args.taxfiltBam and args.lastDb:
        filter_lastal_bam(args.blastnBam, args.lastDb, args.taxfiltBam, JVMmemory=args.JVMmemory)
    return 0


__commands__.append(('deplete_human', parser_deplete_human))



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
    JVMmemory=None
):
    ''' Restrict input reads to those that align to the given
        reference database using LASTAL.
    '''

    with util.file.tmp_dir('-lastdb') as tmp_db_dir:
        # index db if necessary
        lastdb = tools.last.Lastdb()
        if not lastdb.is_indexed(db):
            db = lastdb.build_database(db, os.path.join(tmp_db_dir, 'lastdb'))

        with util.file.tempfname('.read_ids.txt') as hitList:
            # look for lastal hits in BAM and write to temp file
            with open(hitList, 'wt') as outf:
                for read_id in tools.last.Lastal().get_hits(
                        inBam, db,
                        max_gapless_alignments_per_position,
                        min_length_for_initial_matches,
                        max_length_for_initial_matches,
                        max_initial_matches_per_position
                    ):
                    outf.write(read_id + '\n')

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
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_lastal_bam, split_args=True)
    return parser


__commands__.append(('filter_lastal_bam', parser_filter_lastal_bam))


# ==============================
# ***  deplete_bmtagger_bam  ***
# ==============================


def deplete_bmtagger_bam(inBam, db, outBam, threads=None, srprism_memory=7168, JVMmemory=None):
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
    bmtaggerPath = tools.bmtagger.BmtaggerShTool().install_and_get_path()

    # bmtagger calls several executables in the same directory, and blastn;
    # make sure they are accessible through $PATH
    blastnPath = tools.blast.BlastnTool().install_and_get_path()
    path = os.environ['PATH'].split(os.pathsep)
    for t in (bmtaggerPath, blastnPath):
        d = os.path.dirname(t)
        if d not in path:
            path = [d] + path
    path = os.pathsep.join(path)
    os.environ['PATH'] = path

    inReads1 = mkstempfname('.1.fastq')
    tools.samtools.SamtoolsTool().bam2fq(inBam, inReads1)

    bmtaggerConf = mkstempfname('.bmtagger.conf')
    with open(bmtaggerConf, 'w') as f:
        # Default srprismopts: "-b 100000000 -n 5 -R 0 -r 1 -M 7168"
        print('srprismopts="-b 100000000 -n 5 -R 0 -r 1 -M {srprism_memory} --paired false"'.format(srprism_memory=srprism_memory), file=f)
    tempDir = tempfile.mkdtemp()
    matchesFile = mkstempfname('.txt')
    cmdline = [
        bmtaggerPath, '-b', db + '.bitmask', '-C', bmtaggerConf, '-x', db + '.srprism', '-T', tempDir, '-q1',
        '-1', inReads1, '-o', matchesFile
    ]
    log.debug(' '.join(cmdline))
    util.misc.run_and_print(cmdline, check=True)
    os.unlink(inReads1)
    os.unlink(bmtaggerConf)

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
    parser.add_argument('--threads', type=int, default=4, help='The number of threads to use in running blastn.')
    parser.add_argument('--srprismMemory', dest="srprism_memory", type=int, default=7168, help='Memory for srprism.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_bam_bmtagger)
    return parser

def main_deplete_bam_bmtagger(args):
    '''Use bmtagger to deplete input reads against several databases.'''

    def bmtagger_wrapper(inBam, db, outBam, threads, JVMmemory=None):
        return deplete_bmtagger_bam(inBam, db, outBam, threads=threads, srprism_memory=args.srprism_memory, JVMmemory=JVMmemory)

    multi_db_deplete_bam(
        args.inBam,
        args.refDbs,
        bmtagger_wrapper,
        args.outBam,
        threads=args.threads,
        JVMmemory=args.JVMmemory
    )

__commands__.append(('deplete_bam_bmtagger', parser_deplete_bam_bmtagger))


def multi_db_deplete_bam(inBam, refDbs, deplete_method, outBam, threads=1, JVMmemory=None, **kwargs):
    samtools = tools.samtools.SamtoolsTool()
    tmpBamIn = inBam
    for db in refDbs:
        if not samtools.isEmpty(tmpBamIn):
            tmpBamOut = mkstempfname('.bam')
            deplete_method(tmpBamIn, db, tmpBamOut, threads=threads, JVMmemory=JVMmemory, **kwargs)
            if tmpBamIn != inBam:
                os.unlink(tmpBamIn)
            tmpBamIn = tmpBamOut
    shutil.copyfile(tmpBamIn, outBam)

# ========================
# ***  deplete_blastn  ***
# ========================


def run_blastn(blastn_path, db, input_fasta, blast_threads=1):
    """ run blastn on the input fasta file. this is intended to be run in parallel """
    chunk_hits = mkstempfname('.hits.txt.gz')

    blastnCmd = [
        blastn_path, '-db', db, '-word_size', '16', '-num_threads', str(blast_threads), '-evalue', '1e-6', '-outfmt',
        '6', '-max_target_seqs', '1', '-query', input_fasta,
    ]
    log.debug(' '.join(blastnCmd))
    blast_pipe = subprocess.Popen(blastnCmd, stdout=subprocess.PIPE)

    with util.file.open_or_gzopen(chunk_hits, 'wt') as outf:
        # strip tab output to just query read ID names and emit
        last_read_id = None
        for line in blast_pipe.stdout:
            line = line.decode('UTF-8').rstrip('\n\r')
            read_id = line.split('\t')[0]
            # only emit if it is not a duplicate of the previous read ID
            if read_id != last_read_id:
                last_read_id = read_id
                outf.write(read_id + '\n')

    if blast_pipe.poll():
        raise CalledProcessError()
    os.unlink(input_fasta)

    return chunk_hits


def blastn_chunked_fasta(fasta, db, chunkSize=1000000, threads=1):
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

    # get the blastn path here so we don't run conda install checks multiple times
    blastnPath = tools.blast.BlastnTool().install_and_get_path()

    # clamp threadcount to number of CPUs minus one
    threads = max(min(util.misc.available_cpu_count() - 1, threads), 1)

    # determine size of input data; records in fasta file
    number_of_reads = util.file.fasta_length(fasta)
    log.debug("number of reads in fasta file %s" % number_of_reads)
    if number_of_reads == 0:
        return [mkstempfname('.hits.txt')]

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

    log.debug("number of chunk files to be processed by blastn %s" % len(input_fastas))

    hits_files = []
    # run blastn on each of the fasta input chunks
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # if we have so few chunks that there are cpus left over,
        # divide extra cpus evenly among chunks where possible
        # rounding to 1 if there are more chunks than extra threads
        cpus_leftover = (threads - len(input_fastas))
        blast_threads = max(1, int(cpus_leftover / len(input_fastas)))

        futs = [
            executor.submit(functools.partial(run_blastn, blastnPath, db, input_fasta, blast_threads))
            for input_fasta in input_fastas
        ]
        hits_files = [fut.result() for fut in concurrent.futures.as_completed(futs)]

    return hits_files


def deplete_blastn_bam(inBam, db, outBam, threads, chunkSize=1000000, JVMmemory=None):
#def deplete_blastn_bam(inBam, db, outBam, threads, chunkSize=0, JVMmemory=None):
    'Use blastn to remove reads that match at least one of the databases.'

    blast_hits = mkstempfname('.blast_hits.txt')

    if chunkSize:
        ## chunk up input and perform blastn in several parallel threads

        # Initial BAM -> FASTA sequences
        fasta = mkstempfname('.fasta')
        tools.samtools.SamtoolsTool().bam2fa(inBam, fasta)

        # Find BLAST hits
        log.info("running blastn on %s against %s", inBam, db)
        blastOutFiles = blastn_chunked_fasta(fasta, db, chunkSize, threads)
        util.file.cat(blast_hits, blastOutFiles)
        os.unlink(fasta)

    else:
        ## pipe tools together and run blastn multithreaded
        with open(blast_hits, 'wt') as outf:
            for read_id in tools.blast.BlastnTool().get_hits(inBam, db, threads=threads):
                outf.write(read_id + '\n')

    # Deplete BAM of hits
    tools.picard.FilterSamReadsTool().execute(inBam, True, blast_hits, outBam, JVMmemory=JVMmemory)
    os.unlink(blast_hits)


def parser_deplete_blastn_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input BAM file.')
    parser.add_argument('refDbs', nargs='+', help='One or more reference databases for blast.')
    parser.add_argument('outBam', help='Output BAM file with matching reads removed.')
    parser.add_argument('--threads', type=int, default=4, help='The number of threads to use in running blastn.')
    parser.add_argument("--chunkSize", type=int, default=1000000, help='FASTA chunk size (default: %(default)s)')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_deplete_blastn_bam)
    return parser


def main_deplete_blastn_bam(args):
    '''Use blastn to remove reads that match at least one of the specified databases.'''

    def wrapper(inBam, db, outBam, threads, JVMmemory=None):
        return deplete_blastn_bam(inBam, db, outBam, threads=threads, chunkSize=args.chunkSize, JVMmemory=JVMmemory)

    multi_db_deplete_bam(args.inBam, args.refDbs, wrapper, args.outBam, threads=args.threads, JVMmemory=args.JVMmemory)
    return 0


__commands__.append(('deplete_blastn_bam', parser_deplete_blastn_bam))

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

    tools.last.Lastdb().build_database(inputFasta, os.path.join(outputDirectory, outPrefix))


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

# ========================
# ***  blastn_build_db  ***
# ========================


def blastn_build_db(inputFasta, outputDirectory, outputFilePrefix):
    """ Create a database for use with blastn from an input reference FASTA file
    """

    if outputFilePrefix:
        outPrefix = outputFilePrefix
    else:
        baseName = os.path.basename(inputFasta)
        fileNameSansExtension = os.path.splitext(baseName)[0]
        outPrefix = fileNameSansExtension

    blastdb_path = tools.blast.MakeblastdbTool().build_database(inputFasta, os.path.join(outputDirectory, outPrefix))


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


def bmtagger_build_db(inputFasta, outputDirectory, outputFilePrefix):
    """ Create a database for use with Bmtagger from an input FASTA file.
    """
    if outputFilePrefix:
        outPrefix = outputFilePrefix
    else:
        baseName = os.path.basename(inputFasta)
        fileNameSansExtension = os.path.splitext(baseName)[0]
        outPrefix = fileNameSansExtension

    bmtooldb_path = tools.bmtagger.BmtoolTool().build_database(
        inputFasta, os.path.join(outputDirectory, outPrefix + ".bitmask")
    )
    srprismdb_path = tools.bmtagger.SrprismTool().build_database(
        inputFasta, os.path.join(outputDirectory, outPrefix + ".srprism")
    )


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
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, bmtagger_build_db, split_args=True)
    return parser


__commands__.append(('bmtagger_build_db', parser_bmtagger_build_db))

# ========================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
