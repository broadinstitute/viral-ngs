#!/usr/bin/env python
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, os, tempfile, errno, shutil
from Bio import SeqIO
import util.cmd, util.file
import tools.last, tools.prinseq, tools.trimmomatic, tools.bmtagger, \
       tools.blast, tools.mvicuna
from util.file import mkstempfname
from read_utils import fastq_to_fasta

log = logging.getLogger(__name__)

# ==========================
# ***  trim_trimmomatic  ***
# ==========================

def trimmomatic(inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2,
        clipFasta):
    """
    TODO: docstring here
    """
    trimmomaticPath = tools.trimmomatic.TrimmomaticTool() \
        .install_and_get_path()
    tmpUnpaired1 = mkstempfname()
    tmpUnpaired2 = mkstempfname()

    #  This java program has a lot of argments...
    javaCmd = ' '.join( [ 'java -Xmx2g -classpath',
        trimmomaticPath,
        'org.usadellab.trimmomatic.TrimmomaticPE',
        inFastq1,
        inFastq2,
        pairedOutFastq1,
        tmpUnpaired1,
        pairedOutFastq2,
        tmpUnpaired2,
        'LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:30',
        'ILLUMINACLIP:{}:2:30:12'.format(clipFasta)
        ])

    fullCmd = "{javaCmd} && rm {tmpUnpaired1} {tmpUnpaired2}".format(**locals())
    log.info("Running command: {}".format(fullCmd))
    assert not os.system(fullCmd)

def parser_trim_trimmomatic():
    """
    TODO: docstring here
    original clipFasta =
        /idi/sabeti-scratch/kandersen/references/contaminants/contaminants.fasta
    """
    parser = argparse.ArgumentParser(
        description='''Trim read sequences with Trimmomatic.''')
    parser.add_argument("inFastq1", help = "Input reads 1")
    parser.add_argument("inFastq2", help = "Input reads 2")
    parser.add_argument("pairedOutFastq1", help = "Paired output 1")
    parser.add_argument("pairedOutFastq2", help = "Paired output 2")
    parser.add_argument("clipFasta", help = "Fasta file with adapters, PCR "
                        + "sequences, etc. to clip off")

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser

    # Future: handle BAM input and output; handle multiple databases.
    # Will need to implement bam->fastq->bam wrappers that maintain the read
    # metadata
    #parser.add_argument("inBam", help="Input BAM file")
    #parser.add_argument("outBam", help="Output BAM file")

def main_trim_trimmomatic(args):
    '''
        Perhaps move this to a separate script of general bam/alignment utility
        functions?...
    '''
    inFastq1 = args.inFastq1
    inFastq2 = args.inFastq2
    pairedOutFastq1 = args.pairedOutFastq1
    pairedOutFastq2 = args.pairedOutFastq2
    clipFasta = args.clipFasta
    trimmomatic(inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2, clipFasta)
    return 0

__commands__.append(('trim_trimmomatic', main_trim_trimmomatic,
                     parser_trim_trimmomatic))

# =======================
# ***  filter_lastal  ***
# =======================

def filter_lastal(inFastq, refDbs, outFastq):
    """
    TODO: docstring here
    """
    tempFilePath = mkstempfname()
    lastalPath = tools.last.Lastal().install_and_get_path()
    mafSortPath = tools.last.MafSort().install_and_get_path()
    mafConvertPath = tools.last.MafConvert().install_and_get_path()
    prinseqPath = tools.prinseq.PrinseqTool().install_and_get_path()
    noBlastLikeHitsPath = os.path.join( util.file.get_scripts_path(),
                                        'noBlastLikeHits.py')

    # each pipe separated cmd gets own line
    # unfortunately, it doesn't seem to work to do .format(**locals()) on the
    # final string as opposed to the individual parts.
    lastalCmd = ' '.join([
        '{lastalPath} -Q1 {refDbs} {inFastq}'.format(**locals()),
        '| {mafSortPath} -n2'.format(**locals()),
        '| {mafConvertPath} tab /dev/stdin > {tempFilePath}'.format(**locals()),
        ])

    # each option/flag on own line
    noBlastLikeHitsCmd = ' '.join([
        'python', noBlastLikeHitsPath,
            '-b', tempFilePath,
            '-r', inFastq,
            '-m hit' ])

    prinseqCmd = ' '.join([
        'perl', prinseqPath,
            '-ns_max_n 1',
            '-derep 1',
            '-fastq stdin',
            '-out_bad null',
            '-line_width 0',
            '-out_good', outFastq,
        '&& rm', tempFilePath
        ])

    fullCmd = "{lastalCmd} && {noBlastLikeHitsCmd} | {prinseqCmd}" \
        .format(**locals())

    log.debug(fullCmd)
    assert not os.system(fullCmd)

def parser_filter_lastal():
    parser = argparse.ArgumentParser(
        description = '''Restrict input reads to those that align to the given
        reference database using LASTAL.''')
    parser.add_argument("inFastq", help="Input fastq file")
    parser.add_argument("refDbs",
                        help="Reference database to retain from input")
    parser.add_argument("outFastq", help = "Output fastq file")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser

def main_filter_lastal(args):
    inFastq = args.inFastq
    refDbs = args.refDbs
    outFastq = args.outFastq
    filter_lastal(inFastq, refDbs, outFastq)
    return 0
__commands__.append(('filter_lastal', main_filter_lastal, parser_filter_lastal))


# ============================
# ***  partition_bmtagger  ***
# ============================

def select_reads(inFastq, outFastq, selectorFcn) :
    """
    selectorFcn: Bio.SeqRecord.SeqRecord -> bool
    Output in outFastq all reads from inFastq for which
        selectorFcn returns True.
    """
    inFile  = util.file.open_or_gzopen(inFastq)
    outFile = util.file.open_or_gzopen(outFastq, 'w')
    for rec in SeqIO.parse(inFile, 'fastq') :
        if selectorFcn(rec) :
            SeqIO.write([rec], outFile, 'fastq')
    outFile.close()

def partition_bmtagger(inFastq1, inFastq2, databases,
                       outMatch = None, outNoMatch = None) :
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
    blastnDir = os.path.dirname(blastnPath)
    bmtaggerDir = os.path.dirname(bmtaggerPath)
    envStr = 'PATH={bmtaggerDir}:{blastnDir}:$PATH'.format(**locals())
    
    # bmtagger's list of matches strips /1 and /2 from ends of reads
    strip12 = lambda id : id[:-2] if id.endswith('/1') or id.endswith('/2') \
                          else id
    
    tempDir = tempfile.mkdtemp()
    matchesFiles = [mkstempfname() for db in databases]
    curReads1, curReads2 = inFastq1, inFastq2
    for count, (db, matchesFile) in \
            enumerate(zip(databases, matchesFiles)) :
        """
        Loop invariants:
            At the end of the kth loop, curReadsN has the original reads
            depleted by all matches to the first k databases, and
            matchesFiles[:k] contain the list of matching read names.
        """
        cmdline = ('{envStr} {bmtaggerPath} '
                   '-b {db}.bitmask -x {db}.srprism -T {tempDir} '
                   '-q1 -1 {curReads1} -2 {curReads2} '
                   '-o {matchesFile}').format(**locals())
        log.debug(cmdline)
        assert not os.system(cmdline)
        prevReads1, prevReads2 = curReads1, curReads2
        if count < len(databases) - 1 :
            curReads1, curReads2 = mkstempfname(), mkstempfname()
        elif outNoMatch != None :
            # Final time through loop, output depleted to requested files
            curReads1, curReads2 = outNoMatch[0], outNoMatch[1]
        else :
            # No need to calculate final depleted file. No one asked for it.
            # Technically, this violates the loop invariant ;-)
            break
        log.debug("starting select_reads")
        matches = set(line.strip() for line in open(matchesFile))
        noMatchFcn = lambda rec : strip12(rec.id) not in matches
        select_reads(prevReads1, curReads1, noMatchFcn)
        select_reads(prevReads2, curReads2, noMatchFcn)
    if outMatch != None :
        log.debug("preparing outMatch files")
        allMatches = set(line.strip()
                         for matchesFile in matchesFiles
                         for line in open(matchesFile))
        matchFcn = lambda rec : strip12(rec.id) in allMatches
        select_reads(inFastq1, outMatch[0], matchFcn)
        select_reads(inFastq2, outMatch[1], matchFcn)
    log.debug("partition_bmtagger complete")

def parser_partition_bmtagger() :
    parser = argparse.ArgumentParser(
        description='''Use bmtagger to partition input reads into ones that 
            match at least one of several databases and ones that don't match 
            any of the databases.''')
    parser.add_argument('inFastq1',
        help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2',
        help='Input fastq file; 2nd end of paired-end reads.  '\
             'Must have same names as inFastq1')
    parser.add_argument('refDbs', nargs='+',
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.
             ''')
    parser.add_argument('--outMatch', nargs = 2,
        help='Filenames for fastq output of matching reads.')
    parser.add_argument('--outNoMatch', nargs = 2,
        help='Filenames for fastq output of unmatched reads.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser
def main_partition_bmtagger(args) :
    inFastq1 = args.inFastq1
    inFastq2 = args.inFastq2
    databases = args.refDbs
    outMatch = args.outMatch
    outNoMatch = args.outNoMatch
    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)
    return 0
__commands__.append(('partition_bmtagger', main_partition_bmtagger,
                     parser_partition_bmtagger))


# ============================
# ***  dup_remove_mvicuna  ***
# ============================

def dup_remove_mvicuna(inPair, outPair, outUnpaired=None):
    """
    Run mvicuna's duplicate removal operation on paired-end input reads in
        fastq format, producing various outputs in fastq format.
    inPair, pairedOutPair are pairs of file names,
        while unpairedOut is a single file name.
    Notes on weird behaviors of M-Vicuna DupRm:
        For some reason, it requires you to specify both -opfq and -drm_op.
        The -drm_op pair (here, tmp1OutPair) is where the output is initially written.
        Then M-Vicuna renames the -drm_op files to the -opfq filenames (here, tmp2OutPair).
        So obviously, we use throwaway (tempfile) names for -drm_op and use the real
        desired output file names in -opfq.
        The problem is that M-Vicuna uses a rename/move operating system function, which
        means your tempfiles cannot be on a different file system than the final output
        files. This is our typical use case (local disks for tempfile.tempdir and
        network file systems for final output). Hence, our wrapper sets -opfq to yet
        another set of temp file names, and we move the final output files ourselves
        using shutil.move (which is capable of crossing file system boundaries).
    """
    if not outUnpaired:
        outUnpaired = mkstempfname(suffix='.unpaired.fastq')
    tmp1OutPair = (mkstempfname(suffix='.tmp1out.1.fastq'), mkstempfname(suffix='.tmp1out.2.fastq'))
    tmp2OutPair = (mkstempfname(suffix='.tmp2out.1.fastq'), mkstempfname(suffix='.tmp2out.2.fastq'))
    mvicunaPath = tools.mvicuna.MvicunaTool().install_and_get_path()
    input       = ','.join(inPair)
    pairedOut    = ','.join(tmp2OutPair)
    postDupRmOut = ','.join(tmp1OutPair)
    cmdline = ('{mvicunaPath} '
               '-ipfq {input} '
               '-opfq {pairedOut} '
               '-osfq {outUnpaired} '
               '-drm_op {postDupRmOut} '
               '-tasks DupRm').format(**locals())
    log.debug(cmdline)
    assert not os.system(cmdline)
    for tmpfname, outfname in zip(tmp2OutPair, outPair):
        shutil.move(tmpfname, outfname)

def parser_dup_remove_mvicuna() :
    parser = argparse.ArgumentParser(
        description='''Run mvicuna's duplicate removal operation on paired-end 
                       reads.''')
    parser.add_argument('inFastq1',
        help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2',
        help='Input fastq file; 2nd end of paired-end reads.')
    parser.add_argument('pairedOutFastq1',
        help='Output fastq file; 1st end of paired-end reads.')
    parser.add_argument('pairedOutFastq2',
        help='Output fastq file; 2nd end of paired-end reads.')
    parser.add_argument('--unpairedOutFastq',
        default=None,
        help='File name of output unpaired reads')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser
def main_dup_remove_mvicuna(args) :
    inPair = (args.inFastq1, args.inFastq2)
    pairedOutPair = (args.pairedOutFastq1, args.pairedOutFastq2)
    unpairedOutFastq = args.unpairedOutFastq
    dup_remove_mvicuna(inPair, pairedOutPair, unpairedOutFastq)
    return 0

__commands__.append(('dup_remove_mvicuna', main_dup_remove_mvicuna,
                     parser_dup_remove_mvicuna))


# ========================
# ***  deplete_blastn  ***
# ========================

def deplete_blastn(inFastq, outFastq, refDbs) :
    'Use blastn to remove reads that match at least one of the databases.'
    
    ## Get tools
    blastnPath = tools.blast.BlastnTool().install_and_get_path()
    noBlastHits_v3Path = os.path.join(util.file.get_scripts_path(),
                                      'noBlastHits_v3.py')
    
    ## Convert to fasta
    inFasta = mkstempfname() + '.fasta'
    fastq_to_fasta(inFastq, inFasta)

    ## Run blastn using each of the databases in turn
    blastOutFiles = [mkstempfname() for db in refDbs]
    for db, blastOutFile in zip(refDbs, blastOutFiles) :
        blastnCmd = '{blastnPath} -db {db} '                 \
                    '-word_size 16 -evalue 1e-6 -outfmt 6 '  \
                    '-num_descriptions 2 -num_alignments 2 ' \
                    '-query {inFasta} -out {blastOutFile}'.format(**locals())
        log.debug(blastnCmd)
        assert not os.system(blastnCmd)

    ## Combine results from different databases
    blastOutFilesStr = ' '.join(blastOutFiles)
    blastOutCombined = mkstempfname()
    catCmd = 'cat {blastOutFilesStr} > {blastOutCombined}'.format(**locals())
    log.debug(catCmd)
    assert not os.system(catCmd)

    ## run noBlastHits_v3.py to extract reads with no blast hits
    noBlastHitsCmd = 'python {noBlastHits_v3Path} -b {blastOutCombined} ' \
                     '-r {inFastq} -m nohit > {outFastq}'.format(**locals())
    log.debug(noBlastHitsCmd)
    assert not os.system(noBlastHitsCmd)

def parser_deplete_blastn() :
    parser = argparse.ArgumentParser(
        description='''Use blastn to remove reads that match at least
                       one of the specified databases.''')
    parser.add_argument('inFastq',
        help='Input fastq file.')
    parser.add_argument('outFastq',
        help='Output fastq file with matching reads removed.')
    parser.add_argument('refDbs', nargs='+',
        help='One or more reference databases for blast.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser
def main_deplete_blastn(args) :
    inFastq = args.inFastq
    outFastq = args.outFastq
    refDbs = args.refDbs
    deplete_blastn(inFastq, outFastq, refDbs)
    return 0
__commands__.append(('deplete_blastn', main_deplete_blastn,
                     parser_deplete_blastn))

# ========================

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
