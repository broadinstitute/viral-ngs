#!/usr/bin/env python
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, subprocess, os, tempfile, errno, shutil
from Bio import SeqIO
import util.cmd, util.file
import tools, tools.blast
import tools.last, tools.prinseq, tools.trimmomatic, tools.bmtagger
from util.file import mkstempfname
import read_utils

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
    tempdir = tempfile.gettempdir()

    #  This java program has a lot of argments...
    javaCmd = ' '.join( ['java -Xmx2g',
        '-Djava.io.tmpdir={}'.format(tempdir),
        '-classpath',
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
    assert outFastq.endswith('.fastq')
    outFastq = outFastq[:-6]
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
    log.debug(lastalCmd)
    assert not os.system(lastalCmd)

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
            '-out_good', outFastq
        ])

    fullCmd = "{noBlastLikeHitsCmd} | {prinseqCmd}".format(**locals())
    log.debug(fullCmd)
    assert not os.system(fullCmd)
    
    log.debug("done")

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
    with util.file.open_or_gzopen(inFastq, 'rt') as inFile :
        with util.file.open_or_gzopen(outFastq, 'wt') as outFile :
            for rec in SeqIO.parse(inFile, 'fastq') :
                if selectorFcn(rec) :
                    SeqIO.write([rec], outFile, 'fastq')

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
            enumerate(zip(databases, matchesFiles)) :
        """
        Loop invariants:
            At the end of the kth loop, curReadsN has the original reads
            depleted by all matches to the first k databases, and
            matchesFiles[:k] contain the list of matching read names.
        """
        cmdline = [bmtaggerPath,
                   '-b', db+'.bitmask', '-x', db+'.srprism' '-T', tempDir,
                   '-q1', '-1', curReads1, '-2', curReads2,
                   '-o', matchesFile]
        log.debug(' '.join(cmdline))
        subprocess.check_call(cmdline)
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
        cmdline = [bmtaggerPath, '-X',
                   '-b', db+'.bitmask', '-x', db+'.srprism' '-T', tempDir,
                   '-q1', '-1', curReads1, '-2', curReads2,
                   '-o', outprefix]
        log.debug(' '.join(cmdline))
        subprocess.check_call(cmdline)
        curReads1, curReads2 = [outprefix+suffix for suffix in ('_1.fastq', '_2.fastq')]
        tempfiles += [curReads1, curReads2]
    shutil.copyfile(curReads1, outFastq1)
    shutil.copyfile(curReads2, outFastq2)
    map(os.unlink, tempfiles)
    log.debug("deplete_bmtagger complete")

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
    assert outMatch or outNoMatch
    # comment this out until we can figure out why bmtagger -X fails only on Travis
    #if outMatch==None:
    #    deplete_bmtagger(inFastq1, inFastq2, databases, outNoMatch[0], outNoMatch[1])
    #else:
    #    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)
    #return 0
    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)
__commands__.append(('partition_bmtagger', main_partition_bmtagger,
                     parser_partition_bmtagger))



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
    read_utils.fastq_to_fasta(inFastq, inFasta)
    
    ## Run blastn using each of the databases in turn
    blastOutFiles = [mkstempfname() for db in refDbs]
    for db, blastOutFile in zip(refDbs, blastOutFiles) :
        log.info("running blastn on {} against {}".format(inFastq, db))
        blastnCmd = [blastnPath, '-db', db,
                    '-word_size', '16', '-evalue', '1e-6', '-outfmt', '6',
                    '-num_descriptions', '2', '-num_alignments', '2',
                    '-query', inFasta, '-out', blastOutFile]
        log.debug(' '.join(blastnCmd))
        subprocess.check_call(blastnCmd)

    ## Combine results from different databases
    blastOutCombined = mkstempfname('.txt')
    catCmd = ['cat'] + blastOutFiles
    log.debug(' '.join(catCmd) + '> ' + blastOutCombined)
    with open(blastOutCombined, 'wt') as outf:
        subprocess.check_call(catCmd, stdout = outf)

    ## run noBlastHits_v3.py to extract reads with no blast hits
    # TODO: slurp the small amount of code in this script into here
    noBlastHitsCmd = ['python', noBlastHits_v3Path, '-b', blastOutCombined,
                     '-r', inFastq, '-m', 'nohit']
    log.debug(' '.join(noBlastHitsCmd) + '> ' + outFastq)
    with util.file.open_or_gzopen(outFastq, 'wt') as outf :
        subprocess.check_call(noBlastHitsCmd, stdout = outf)

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

def deplete_blastn_paired(infq1, infq2, outfq1, outfq2, refDbs):
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
    deplete_blastn(tmpfq2_b, tmpfq2_c, refDbs)
    # purge fq1 of read pairs lost in fq2
    read_utils.purge_unmated(tmpfq1_b, tmpfq2_c, outfq1, outfq2)

def parser_deplete_blastn_paired() :
    parser = argparse.ArgumentParser(
        description='''Use blastn to remove reads that match at least
                       one of the specified databases.''')
    parser.add_argument('inFastq1',
        help='Input fastq file.')
    parser.add_argument('inFastq2',
        help='Input fastq file.')
    parser.add_argument('outFastq1',
        help='Output fastq file with matching reads removed.')
    parser.add_argument('outFastq2',
        help='Output fastq file with matching reads removed.')
    parser.add_argument('refDbs', nargs='+',
        help='One or more reference databases for blast.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser
def main_deplete_blastn_paired(args) :
    deplete_blastn_paired(args.inFastq1, args.inFastq2,
        args.outFastq1, args.outFastq2, args.refDbs)
    return 0
__commands__.append(('deplete_blastn_paired', main_deplete_blastn_paired,
                     parser_deplete_blastn_paired))

# ========================

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
