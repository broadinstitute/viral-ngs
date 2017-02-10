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
import tools.trimmomatic
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

    multi_db_deplete_bam(
        bamToDeplete,
        args.bmtaggerDbs,
        deplete_bmtagger_bam,
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

# ==========================
# ***  trim_trimmomatic  ***
# ==========================


def trimmomatic(
    inFastq1,
    inFastq2,
    pairedOutFastq1,
    pairedOutFastq2,
    clipFasta,
    unpairedOutFastq1=None,
    unpairedOutFastq2=None,
    leading_q_cutoff=15,
    trailing_q_cutoff=15,
    minlength_to_keep=30,
    sliding_window_size=4,
    sliding_window_q_cutoff=25
):
    '''Trim read sequences with Trimmomatic.'''
    trimmomaticPath = tools.trimmomatic.TrimmomaticTool().install_and_get_path()
    unpairedFastq1 = unpairedOutFastq1 or mkstempfname()
    unpairedFastq2 = unpairedOutFastq2 or mkstempfname()

    javaCmd = []

    # the conda version wraps the jar file with a shell script
    if trimmomaticPath.endswith(".jar"):
        #  This java program has a lot of argments...
        javaCmd.extend(
            [
                'java', '-Xmx2g', '-Djava.io.tmpdir=' + tempfile.tempdir, '-classpath', trimmomaticPath,
                'org.usadellab.trimmomatic.TrimmomaticPE'
            ]
        )
    else:
        javaCmd.extend([trimmomaticPath, "PE"])

    # Explicitly use Phred-33 quality scores
    javaCmd.extend(['-phred33'])

    javaCmd.extend(
        [
            inFastq1, inFastq2, pairedOutFastq1, unpairedFastq1, pairedOutFastq2, unpairedFastq2,
            'LEADING:{leading_q_cutoff}'.format(leading_q_cutoff=leading_q_cutoff),
            'TRAILING:{trailing_q_cutoff}'.format(trailing_q_cutoff=trailing_q_cutoff),
            'SLIDINGWINDOW:{sliding_window_size}:{sliding_window_q_cutoff}'.format(
                sliding_window_size=sliding_window_size,
                sliding_window_q_cutoff=sliding_window_q_cutoff,
            ), 
            'MINLEN:{minlength_to_keep}'.format(minlength_to_keep=minlength_to_keep),
            'ILLUMINACLIP:{clipFasta}:2:30:12'.format(clipFasta=clipFasta)
        ]
    )

    log.debug(' '.join(javaCmd))
    util.misc.run_and_print(javaCmd, check=True)

    if not unpairedOutFastq1:
        os.unlink(unpairedFastq1)
    if not unpairedOutFastq2:
        os.unlink(unpairedFastq2)


def parser_trim_trimmomatic(parser=argparse.ArgumentParser()):
    parser.add_argument("inFastq1", help="Input reads 1")
    parser.add_argument("inFastq2", help="Input reads 2")
    parser.add_argument("pairedOutFastq1", help="Paired output 1")
    parser.add_argument("pairedOutFastq2", help="Paired output 2")
    parser.add_argument("--unpairedOutFastq1", help="Unpaired output 1 (default: write to temp and discard)")
    parser.add_argument("--unpairedOutFastq2", help="Unpaired output 2 (default: write to temp and discard)")
    parser.add_argument(
        '--leadingQCutoff',
        dest="leading_q_cutoff",
        help='minimum quality required to keep a base on the leading end (default: %(default)s)',
        type=int,
        default=15
    )
    parser.add_argument(
        '--trailingQCutoff',
        dest="trailing_q_cutoff",
        help='minimum quality required to keep a base on the trailing end (default: %(default)s)',
        type=int,
        default=15
    )
    parser.add_argument(
        '--minlengthToKeep',
        dest="minlength_to_keep",
        help='minimum length of reads to be kept (default: %(default)s)',
        type=int,
        default=30
    )
    parser.add_argument(
        '--slidingWindowSize',
        dest="sliding_window_size",
        help='the number of bases to average across (default: %(default)s)',
        type=int,
        default=4
    )
    parser.add_argument(
        '--slidingWindowQCutoff',
        dest="sliding_window_q_cutoff",
        help='specifies the average quality required in the sliding window (default: %(default)s)',
        type=int,
        default=25
    )
    parser.add_argument("clipFasta", help="Fasta file with adapters, PCR sequences, etc. to clip off")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, trimmomatic, split_args=True)
    return parser


__commands__.append(('trim_trimmomatic', parser_trim_trimmomatic))

# =======================
# ***  filter_lastal  ***
# =======================


def lastal_chunked_fastq(
    inFastq,
    db,
    outFastq,
    max_gapless_alignments_per_position=1,
    min_length_for_initial_matches=5,
    max_length_for_initial_matches=50,
    max_initial_matches_per_position=100,
    chunk_size=100000
):

    lastal_path = tools.last.Lastal().install_and_get_path()
    mafsort_path = tools.last.MafSort().install_and_get_path()
    mafconvert_path = tools.last.MafConvert().install_and_get_path()
    no_blast_like_hits_path = os.path.join(util.file.get_scripts_path(), 'noBlastLikeHits.py')

    filtered_fastq_files = []
    with open(inFastq, "rt") as fastqFile:
        record_iter = SeqIO.parse(fastqFile, "fastq")
        for batch in util.misc.batch_iterator(record_iter, chunk_size):

            chunk_fastq = mkstempfname('.fastq')
            with open(chunk_fastq, "wt") as handle:
                SeqIO.write(batch, handle, "fastq")
            batch = None

            lastal_out = mkstempfname('.lastal')
            with open(lastal_out, 'wt') as outf:
                cmd = [lastal_path, '-Q1', '-P0']
                cmd.extend(
                    [
                        '-n', max_gapless_alignments_per_position, '-l', min_length_for_initial_matches, '-L',
                        max_length_for_initial_matches, '-m', max_initial_matches_per_position
                    ]
                )
                cmd = [str(x) for x in cmd]
                cmd.extend([db, chunk_fastq])
                log.debug(' '.join(cmd) + ' > ' + lastal_out)
                util.misc.run_and_save(cmd, outf=outf)
            # everything below this point in this method should be replaced with
            # our own code that just reads lastal output and makes a list of read names

            mafsort_out = mkstempfname('.mafsort')
            with open(mafsort_out, 'wt') as outf:
                with open(lastal_out, 'rt') as inf:
                    cmd = [mafsort_path, '-n2']
                    log.debug('cat ' + lastal_out + ' | ' + ' '.join(cmd) + ' > ' + mafsort_out)
                    subprocess.check_call(cmd, stdin=inf, stdout=outf)
            os.unlink(lastal_out)

            mafconvert_out = mkstempfname('.mafconvert')
            with open(mafconvert_out, 'wt') as outf:
                cmd = ["python", mafconvert_path, 'tab', mafsort_out]
                log.debug(' '.join(cmd) + ' > ' + mafconvert_out)
                subprocess.check_call(cmd, stdout=outf)
            os.unlink(mafsort_out)

            filtered_fastq_chunk = mkstempfname('.filtered.fastq')
            with open(filtered_fastq_chunk, 'wt') as outf:
                cmd = [no_blast_like_hits_path, '-b', mafconvert_out, '-r', chunk_fastq, '-m', 'hit']
                log.debug(' '.join(cmd) + ' > ' + filtered_fastq_chunk)
                subprocess.check_call(cmd, stdout=outf)
                filtered_fastq_files.append(filtered_fastq_chunk)
            os.unlink(mafconvert_out)

    # concatenate filtered fastq files to outFastq
    util.file.concat(filtered_fastq_files, outFastq)

    # remove temp fastq files
    for tempfastq in filtered_fastq_files:
        os.unlink(tempfastq)


def lastal_get_hits(
    inFastq,
    db,
    outList,
    max_gapless_alignments_per_position=1,
    min_length_for_initial_matches=5,
    max_length_for_initial_matches=50,
    max_initial_matches_per_position=100
):
    filteredFastq = mkstempfname('.filtered.fastq')
    lastal_chunked_fastq(
        inFastq,
        db,
        filteredFastq,
        max_gapless_alignments_per_position=max_gapless_alignments_per_position,
        min_length_for_initial_matches=min_length_for_initial_matches,
        max_length_for_initial_matches=max_length_for_initial_matches,
        max_initial_matches_per_position=max_initial_matches_per_position
    )

    with open(outList, 'wt') as outf:
        with open(filteredFastq, 'rt') as inf:
            line_num = 0
            for line in inf:
                if (line_num % 4) == 0:
                    seq_id = line.rstrip('\n\r')[1:]
                    if seq_id.endswith('/1') or seq_id.endswith('/2'):
                        seq_id = seq_id[:-2]
                    outf.write(seq_id + '\n')
                line_num += 1

    os.unlink(filteredFastq)


def parser_lastal_generic(parser=argparse.ArgumentParser()):
    # max_gapless_alignments_per_position, min_length_for_initial_matches, max_length_for_initial_matches, max_initial_matches_per_position
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
    return parser


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
    # convert BAM to paired FASTQ
    inReads1 = mkstempfname('.1.fastq')
    inReads2 = mkstempfname('.2.fastq')
    tools.samtools.SamtoolsTool().bam2fq(inBam, inReads1, inReads2)

    # look for hits in inReads1 and inReads2
    hitList1 = mkstempfname('.1.hits')
    hitList2 = mkstempfname('.2.hits')
    lastal_get_hits(
        inReads1, db, hitList1, max_gapless_alignments_per_position, min_length_for_initial_matches,
        max_length_for_initial_matches, max_initial_matches_per_position
    )
    os.unlink(inReads1)
    lastal_get_hits(
        inReads2, db, hitList2, max_gapless_alignments_per_position, min_length_for_initial_matches,
        max_length_for_initial_matches, max_initial_matches_per_position
    )
    os.unlink(inReads2)

    # merge hits
    hitList = mkstempfname('.hits')
    with open(hitList, 'wt') as outf:
        subprocess.check_call(['sort', '-u', hitList1, hitList2], stdout=outf)
    os.unlink(hitList1)
    os.unlink(hitList2)

    # filter original BAM file against keep list
    tools.picard.FilterSamReadsTool().execute(inBam, False, hitList, outBam, JVMmemory=JVMmemory)
    os.unlink(hitList)


def parser_filter_lastal_bam(parser=argparse.ArgumentParser()):
    parser = parser_lastal_generic(parser)
    parser.add_argument("inBam", help="Input reads")
    parser.add_argument("db", help="Database of taxa we keep")
    parser.add_argument("outBam", help="Output reads, filtered to refDb")
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_lastal_bam, split_args=True)
    return parser


__commands__.append(('filter_lastal_bam', parser_filter_lastal_bam))


def filter_lastal(
    inFastq,
    refDb,
    outFastq,
    max_gapless_alignments_per_position=1,
    min_length_for_initial_matches=5,
    max_length_for_initial_matches=50,
    max_initial_matches_per_position=100
):
    ''' Restrict input reads to those that align to the given
        reference database using LASTAL. Also, remove duplicates with prinseq.
    '''
    assert outFastq.endswith('.fastq')

    filtered_fastq = mkstempfname('.filtered.fastq')
    lastal_chunked_fastq(
        inFastq,
        refDb,
        filtered_fastq,
        max_gapless_alignments_per_position=max_gapless_alignments_per_position,
        min_length_for_initial_matches=min_length_for_initial_matches,
        max_length_for_initial_matches=max_length_for_initial_matches,
        max_initial_matches_per_position=max_initial_matches_per_position
    )

    # remove duplicate reads and reads with multiple Ns
    if os.path.getsize(filtered_fastq) == 0:
        # prinseq-lite fails on empty file input (which can happen in real life
        # if no reads match the refDb) so handle this scenario specially
        log.info("output is empty: no reads in input match refDb")
        shutil.copyfile(filtered_fastq, outFastq)
    else:
        prinseqPath = tools.prinseq.PrinseqTool().install_and_get_path()
        prinseqCmd = [
            'perl', prinseqPath, '-ns_max_n', '1', '-derep', '1', '-fastq', filtered_fastq, '-out_bad', 'null',
            '-line_width', '0', '-out_good', outFastq[:-6]
        ]
        log.debug(' '.join(prinseqCmd))
        util.misc.run_and_print(prinseqCmd, check=True)
    os.unlink(filtered_fastq)


def parser_filter_lastal(parser=argparse.ArgumentParser()):
    parser = parser_lastal_generic(parser)
    parser.add_argument("inFastq", help="Input fastq file")
    parser.add_argument("refDb", help="Reference database to retain from input")
    parser.add_argument("outFastq", help="Output fastq file")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_lastal, split_args=True)
    return parser


__commands__.append(('filter_lastal', parser_filter_lastal))

# ============================
# ***  partition_bmtagger  ***
# ============================


def deplete_bmtagger_bam(inBam, db, outBam, threads=None, JVMmemory=None):
    """
    Use bmtagger to partition the input reads into ones that match at least one
        of the databases and ones that don't match any of the databases.
    inBam: paired-end input reads in BAM format.
    db: bmtagger expects files
        db.bitmask created by bmtool, and
        db.srprism.idx, db.srprism.map, etc. created by srprism mkindex
    outBam: the output BAM files to hold the unmatched reads.
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
    inReads2 = mkstempfname('.2.fastq')
    tools.samtools.SamtoolsTool().bam2fq(inBam, inReads1, inReads2)

    bmtaggerConf = mkstempfname('.bmtagger.conf')
    with open(bmtaggerConf, 'w') as f:
        # Default srprismopts: "-b 100000000 -n 5 -R 0 -r 1 -M 7168"
        print('srprismopts="-b 100000000 -n 5 -R 0 -r 1 -M 7168 --paired false"', file=f)
    tempDir = tempfile.mkdtemp()
    matchesFile = mkstempfname('.txt')
    cmdline = [
        bmtaggerPath, '-b', db + '.bitmask', '-C', bmtaggerConf, '-x', db + '.srprism', '-T', tempDir, '-q1', '-1',
        inReads1, '-2', inReads2, '-o', matchesFile
    ]
    log.debug(' '.join(cmdline))
    util.misc.run_and_print(cmdline, check=True)

    tools.picard.FilterSamReadsTool().execute(inBam, True, matchesFile, outBam, JVMmemory=JVMmemory)
    os.unlink(matchesFile)


def select_reads(inFastq, outFastq, selectorFcn):
    """
        selectorFcn: Bio.SeqRecord.SeqRecord -> bool
        Output in outFastq all reads from inFastq for which
            selectorFcn returns True.
        TO DO: change this to use Picard FilterSamReads (and operate
            on BAM files) which is likely much faster. This is the
            slowest step of partition_bmtagger currently.
    """
    with util.file.open_or_gzopen(inFastq, 'rt') as inFile:
        with util.file.open_or_gzopen(outFastq, 'wt') as outFile:
            for rec in SeqIO.parse(inFile, 'fastq'):
                if selectorFcn(rec):
                    SeqIO.write([rec], outFile, 'fastq')


def partition_bmtagger(inFastq1, inFastq2, databases, outMatch=None, outNoMatch=None):
    """
    Use bmtagger to partition the input reads into ones that match at least one
        of the databases and ones that don't match any of the databases.
    inFastq1, inFastq2: paired-end input reads in fastq format
        The names of the reads must be in one-to-one correspondence.
    databases: for each db in databases bmtagger expects files
        db.bitmask created by bmtool, and
        db.srprism.idx, db.srprism.map, etc. created by srprism mkindex
    outMatch, outNoMatch: either may be None, otherwise a pair of files to
        hold the matching or unmatched reads.
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

    # bmtagger's list of matches strips /1 and /2 from ends of reads
    strip12 = lambda id : id[:-2] if id.endswith('/1') or id.endswith('/2') \
        else id

    tempDir = tempfile.mkdtemp()
    matchesFiles = [mkstempfname() for db in databases]
    curReads1, curReads2 = inFastq1, inFastq2
    for count, (db, matchesFile) in \
            enumerate(zip(databases, matchesFiles)):

        # Loop invariants:
        #     At the end of the kth loop, curReadsN has the original reads
        #     depleted by all matches to the first k databases, and
        #     matchesFiles[:k] contain the list of matching read names.

        cmdline = [
            bmtaggerPath, '-b', db + '.bitmask', '-x', db + '.srprism', '-T', tempDir, '-q1', '-1', curReads1, '-2',
            curReads2, '-o', matchesFile
        ]
        log.debug(' '.join(cmdline))
        util.misc.run_and_print(cmdline, check=True)
        prevReads1, prevReads2 = curReads1, curReads2
        if count < len(databases) - 1:
            curReads1, curReads2 = mkstempfname(), mkstempfname()
        elif outNoMatch is not None:
            # Final time through loop, output depleted to requested files
            curReads1, curReads2 = outNoMatch[0], outNoMatch[1]
        else:
            # No need to calculate final depleted file. No one asked for it.
            # Technically, this violates the loop invariant ;-)
            break
        log.debug("starting select_reads")
        with open(matchesFile) as inf:
            matches = set(line.strip() for line in inf)
        noMatchFcn = lambda rec: strip12(rec.id) not in matches
        select_reads(prevReads1, curReads1, noMatchFcn)
        select_reads(prevReads2, curReads2, noMatchFcn)
    if outMatch is not None:
        log.debug("preparing outMatch files")
        allMatches = set()
        for matchesFile in matchesFiles:
            with open(matchesFile) as inf:
                newMatches = set(line.strip() for line in inf)
            allMatches = allMatches.union(newMatches)
        matchFcn = lambda rec: strip12(rec.id) in allMatches
        select_reads(inFastq1, outMatch[0], matchFcn)
        select_reads(inFastq2, outMatch[1], matchFcn)
    log.debug("partition_bmtagger complete")


def deplete_bmtagger(inFastq1, inFastq2, databases, outFastq1, outFastq2):
    """
    Use bmtagger to partition the input reads into ones that match at least one
        of the databases and ones that don't match any of the databases.
    inFastq1, inFastq2: paired-end input reads in fastq format
        The names of the reads must be in one-to-one correspondence.
    databases: for each db in databases bmtagger expects files
        db.bitmask created by bmtool, and
        db.srprism.idx, db.srprism.map, etc. created by srprism mkindex
    outFastq1, outFastq2: pair of output fastq files depleted of reads present
        in the databases
    This version is optimized for the case of only requiring depletion, which
    allows us to avoid time-intensive lookups.
    """
    bmtaggerPath = tools.bmtagger.BmtaggerShTool().install_and_get_path()
    blastnPath = tools.blast.BlastnTool().install_and_get_path()

    # bmtagger calls several executables in the same directory, and blastn;
    # make sure they are accessible through $PATH
    path = os.environ['PATH'].split(os.pathsep)
    for t in (bmtaggerPath, blastnPath):
        d = os.path.dirname(t)
        if d not in path:
            path = [d] + path
    path = os.pathsep.join(path)
    os.environ['PATH'] = path

    tempDir = tempfile.mkdtemp()
    curReads1, curReads2 = inFastq1, inFastq2
    tempfiles = []
    for db in databases:
        outprefix = mkstempfname()
        cmdline = [
            bmtaggerPath, '-X', '-b', db + '.bitmask', '-x', db + '.srprism', '-T', tempDir, '-q1', '-1', curReads1,
            '-2', curReads2, '-o', outprefix
        ]
        log.debug(' '.join(cmdline))
        util.misc.run_and_print(cmdline, check=True)
        curReads1, curReads2 = [outprefix + suffix for suffix in ('_1.fastq', '_2.fastq')]
        tempfiles += [curReads1, curReads2]
    shutil.copyfile(curReads1, outFastq1)
    shutil.copyfile(curReads2, outFastq2)
    for fn in tempfiles:
        os.unlink(fn)
    log.debug("deplete_bmtagger complete")


def parser_partition_bmtagger(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq1', help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument(
        'inFastq2',
        help='Input fastq file; 2nd end of paired-end reads.  '
        'Must have same names as inFastq1'
    )
    parser.add_argument(
        'refDbs',
        nargs='+',
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.
             '''
    )
    parser.add_argument('--outMatch', nargs=2, help='Filenames for fastq output of matching reads.')
    parser.add_argument('--outNoMatch', nargs=2, help='Filenames for fastq output of unmatched reads.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_partition_bmtagger)
    return parser


def main_partition_bmtagger(args):
    ''' Use bmtagger to partition input reads into ones that
        match at least one of several databases and ones that don't match
        any of the databases.
    '''
    inFastq1 = args.inFastq1
    inFastq2 = args.inFastq2
    databases = args.refDbs
    outMatch = args.outMatch
    outNoMatch = args.outNoMatch
    assert outMatch or outNoMatch
    # comment this out until we can figure out why bmtagger -X fails only on Travis
    # if outMatch==None:
    #    deplete_bmtagger(inFastq1, inFastq2, databases, outNoMatch[0], outNoMatch[1])
    # else:
    #    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)
    # return 0
    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)


__commands__.append(('partition_bmtagger', parser_partition_bmtagger))


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
    multi_db_deplete_bam(
        args.inBam,
        args.refDbs,
        deplete_bmtagger_bam,
        args.outBam,
        threads=args.threads,
        JVMmemory=args.JVMmemory
    )


__commands__.append(('deplete_bam_bmtagger', parser_deplete_bam_bmtagger))


def multi_db_deplete_bam(inBam, refDbs, deplete_method, outBam, threads=1, JVMmemory=None):
    samtools = tools.samtools.SamtoolsTool()
    tmpBamIn = inBam
    for db in refDbs:
        if not samtools.isEmpty(tmpBamIn):
            tmpBamOut = mkstempfname('.bam')
            deplete_method(tmpBamIn, db, tmpBamOut, threads=threads, JVMmemory=JVMmemory)
            if tmpBamIn != inBam:
                os.unlink(tmpBamIn)
            tmpBamIn = tmpBamOut
    shutil.copyfile(tmpBamIn, outBam)

# ========================
# ***  deplete_blastn  ***
# ========================


def run_blastn(blastn_path, db, input_fasta, blast_threads=1):
    """ run blastn on the input fasta file. this is intended to be run in parallel """
    chunk_hits = mkstempfname('.hits.txt')
    blastnCmd = [
        blastn_path, '-db', db, '-word_size', '16', '-num_threads', str(blast_threads), '-evalue', '1e-6', '-outfmt',
        '6', '-max_target_seqs', '2', '-query', input_fasta, '-out', chunk_hits
    ]
    log.debug(' '.join(blastnCmd))
    util.misc.run_and_print(blastnCmd, check=True)

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


def no_blast_hits(blastOutCombined, inFastq, outFastq):
    '''
        outputs to outFastq: reads that have no blast hits
    '''

    blastReads = {}
    with open(blastOutCombined, 'r') as blastFile:
        for line in blastFile:
            blastReads[(line[0:line.find('\t')])] = True

    with util.file.open_or_gzopen(outFastq, 'wt') as outf:
        with open(inFastq, 'r') as readsFile:
            nohit = True
            isFastq = inFastq.endswith('.fastq')
            while True:
                line1 = readsFile.readline()
                line2 = readsFile.readline()
                if not line2:
                    break
                line3 = ''
                line4 = ''
                if isFastq:
                    line3 = readsFile.readline()
                    if not line3:
                        break
                    line4 = readsFile.readline()
                    if not line4:
                        break
                if nohit != (line1[1:line1.find('\n')] in blastReads):
                    outf.write(line1 + line2 + line3 + line4)


def deplete_blastn(inFastq, outFastq, refDbs, threads=1, chunkSize=1000000):
    'Use blastn to remove reads that match at least one of the databases.'

    # Convert to fasta
    inFasta = mkstempfname('.fasta')
    read_utils.fastq_to_fasta(inFastq, inFasta)

    # Run blastn using each of the databases in turn
    blastOutFiles = []
    for db in refDbs:
        log.info("running blastn on %s against %s", inFastq, db)
        blastOutFiles += blastn_chunked_fasta(inFasta, db, chunkSize, threads)

    # Combine results from different databases
    blastOutCombined = mkstempfname('.txt')
    catCmd = ['cat'] + blastOutFiles
    log.debug(' '.join(catCmd) + '> ' + blastOutCombined)
    with open(blastOutCombined, 'wt') as outf:
        subprocess.check_call(catCmd, stdout=outf)

    # extract reads with no blast hits
    no_blast_hits(blastOutCombined, inFastq, outFastq)


def parser_deplete_blastn(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq', help='Input fastq file.')
    parser.add_argument('outFastq', help='Output fastq file with matching reads removed.')
    parser.add_argument('refDbs', nargs='+', help='One or more reference databases for blast.')
    parser.add_argument('--threads', type=int, default=4, help='The number of threads to use in running blastn.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, deplete_blastn, split_args=True)
    return parser


__commands__.append(('deplete_blastn', parser_deplete_blastn))


def deplete_blastn_paired(infq1, infq2, outfq1, outfq2, refDbs, threads):
    'Use blastn to remove reads that match at least one of the databases.'
    tmpfq1_a = mkstempfname('.fastq')
    tmpfq1_b = mkstempfname('.fastq')
    tmpfq2_b = mkstempfname('.fastq')
    tmpfq2_c = mkstempfname('.fastq')
    # deplete fq1
    deplete_blastn(infq1, tmpfq1_a, refDbs)
    # purge fq2 of read pairs lost in fq1
    # (this should significantly speed up the second run of deplete_blastn)
    read_utils.purge_unmated(tmpfq1_a, infq2, tmpfq1_b, tmpfq2_b)
    # deplete fq2
    deplete_blastn(tmpfq2_b, tmpfq2_c, refDbs, threads)
    # purge fq1 of read pairs lost in fq2
    read_utils.purge_unmated(tmpfq1_b, tmpfq2_c, outfq1, outfq2)


def parser_deplete_blastn_paired(parser=argparse.ArgumentParser()):
    parser.add_argument('infq1', help='Input fastq file.')
    parser.add_argument('infq2', help='Input fastq file.')
    parser.add_argument('outfq1', help='Output fastq file with matching reads removed.')
    parser.add_argument('outfq2', help='Output fastq file with matching reads removed.')
    parser.add_argument('refDbs', nargs='+', help='One or more reference databases for blast.')
    parser.add_argument('--threads', type=int, default=4, help='The number of threads to use in running blastn.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, deplete_blastn_paired, split_args=True)
    return parser


__commands__.append(('deplete_blastn_paired', parser_deplete_blastn_paired))


def deplete_blastn_bam(inBam, db, outBam, threads, chunkSize=1000000, JVMmemory=None):
    'Use blastn to remove reads that match at least one of the databases.'

    #blastnPath = tools.blast.BlastnTool().install_and_get_path()
    fastq1 = mkstempfname('.1.fastq')
    fastq2 = mkstempfname('.2.fastq')
    fasta = mkstempfname('.1.fasta')
    blast_hits = mkstempfname('.blast_hits.txt')
    halfBam = mkstempfname('.half.bam')
    blastOutFile = mkstempfname('.hits.txt')

    # Initial BAM -> FASTQ pair
    tools.samtools.SamtoolsTool().bam2fq(inBam, fastq1, fastq2)

    # Find BLAST hits against FASTQ1
    read_utils.fastq_to_fasta(fastq1, fasta)
    os.unlink(fastq1)
    os.unlink(fastq2)
    log.info("running blastn on %s pair 1 against %s", inBam, db)
    blastOutFiles = blastn_chunked_fasta(fasta, db, chunkSize, threads)
    with open(blast_hits, 'wt') as outf:
        for blastOutFile in blastOutFiles:
            with open(blastOutFile, 'rt') as inf:
                for line in inf:
                    idVal = line.split('\t')[0].strip()
                    if idVal.endswith('/1') or idVal.endswith('/2'):
                        idVal = idVal[:-2]
                    outf.write(idVal + '\n')
            os.unlink(blastOutFile)

    # Deplete BAM of hits in FASTQ1
    tools.picard.FilterSamReadsTool().execute(inBam, True, blast_hits, halfBam, JVMmemory=JVMmemory)

    # Depleted BAM -> FASTQ pair
    tools.samtools.SamtoolsTool().bam2fq(halfBam, fastq1, fastq2)

    # Find BLAST hits against FASTQ2 (which is already smaller than before)
    read_utils.fastq_to_fasta(fastq2, fasta)
    os.unlink(fastq1)
    os.unlink(fastq2)
    log.info("running blastn on %s pair 2 against %s", inBam, db)
    blastOutFiles = blastn_chunked_fasta(fasta, db, chunkSize, threads)
    with open(blast_hits, 'wt') as outf:
        for blastOutFile in blastOutFiles:
            with open(blastOutFile, 'rt') as inf:
                for line in inf:
                    idVal = line.split('\t')[0].strip()
                    if idVal.endswith('/1') or idVal.endswith('/2'):
                        idVal = idVal[:-2]
                    outf.write(idVal + '\n')
            os.unlink(blastOutFile)

    # Deplete BAM of hits against FASTQ2
    tools.picard.FilterSamReadsTool().execute(halfBam, True, blast_hits, outBam, JVMmemory=JVMmemory)

    # Clean up
    for fn in (fasta, blast_hits, halfBam):
        os.unlink(fn)


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
