#!/usr/bin/env python
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"
__commands__ = []

import argparse, logging, subprocess, os, tempfile, errno, shutil
from Bio import SeqIO
import util.cmd, util.file
import tools, tools.blast
import tools.last, tools.prinseq, tools.trimmomatic, tools.bmtagger, tools.picard
from util.file import mkstempfname
import read_utils

log = logging.getLogger(__name__)

# =======================
# ***  deplete_human  ***
# =======================

def parser_deplete_human(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam',
        help='Input BAM file.')
    parser.add_argument('revertBam',
        help='Output BAM file.')
    parser.add_argument('bmtaggerBam',
        help='Output BAM file.')
    parser.add_argument('rmdupBam',
        help='Output BAM file.')
    parser.add_argument('blastnBam',
        help='Output BAM file.')
    parser.add_argument('--taxfiltBam',
        help='Output BAM file.',
        default=None)
    parser.add_argument('--bmtaggerDbs', nargs='+', required=True,
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.''')
    parser.add_argument('--blastDbs', nargs='+', required=True,
        help='One or more reference databases for blast.')
    parser.add_argument('--lastDb',
        help='One reference database for last.',
        default=None)
    parser.add_argument('--JVMmemory', default = tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_deplete_human)
    return parser
def main_deplete_human(args):
    '''Run the entire depletion pipeline: bmtagger, mvicuna, blastn, and maybe lastal'''
    tools.picard.RevertSamTool().execute(args.inBam, args.revertBam,
        picardOptions=['SORT_ORDER=queryname', 'SANITIZE=true'])
    multi_db_deplete_bam(args.revertBam, args.bmtaggerDbs, deplete_bmtagger_bam, args.bmtaggerBam, JVMmemory=args.JVMmemory)
    read_utils.rmdup_mvicuna_bam(args.bmtaggerBam, args.rmdupBam, JVMmemory=args.JVMmemory)
    multi_db_deplete_bam(args.rmdupBam, args.blastDbs, deplete_blastn_bam, args.blastnBam, JVMmemory=args.JVMmemory)
    if args.taxfiltBam and args.lastDb:
        filter_lastal_bam(args.blastnBam, args.lastDb, args.taxfiltBam, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('deplete_human', parser_deplete_human))


# ==========================
# ***  trim_trimmomatic  ***
# ==========================

def trimmomatic(inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2,
        clipFasta):
    '''Trim read sequences with Trimmomatic.'''
    trimmomaticPath = tools.trimmomatic.TrimmomaticTool().install_and_get_path()
    tmpUnpaired1 = mkstempfname()
    tmpUnpaired2 = mkstempfname()

    #  This java program has a lot of argments...
    javaCmd = ['java', '-Xmx2g',
        '-Djava.io.tmpdir='+tempfile.tempdir,
        '-classpath',
        trimmomaticPath,
        'org.usadellab.trimmomatic.TrimmomaticPE',
        inFastq1,
        inFastq2,
        pairedOutFastq1,
        tmpUnpaired1,
        pairedOutFastq2,
        tmpUnpaired2,
        'LEADING:20', 'TRAILING:20', 'SLIDINGWINDOW:4:25', 'MINLEN:30',
        'ILLUMINACLIP:{}:2:30:12'.format(clipFasta)
        ]

    log.debug(' '.join(javaCmd))
    subprocess.check_call(javaCmd)
    os.unlink(tmpUnpaired1)
    os.unlink(tmpUnpaired2)

def parser_trim_trimmomatic(parser=argparse.ArgumentParser()):
    parser.add_argument("inFastq1", help = "Input reads 1")
    parser.add_argument("inFastq2", help = "Input reads 2")
    parser.add_argument("pairedOutFastq1", help = "Paired output 1")
    parser.add_argument("pairedOutFastq2", help = "Paired output 2")
    parser.add_argument("clipFasta",
        help = "Fasta file with adapters, PCR sequences, etc. to clip off")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, trimmomatic, split_args=True)
    return parser
__commands__.append(('trim_trimmomatic', parser_trim_trimmomatic))

# =======================
# ***  filter_lastal  ***
# =======================


def lastal_get_hits(inFastq, db, outList):
    lastalPath = tools.last.Lastal().install_and_get_path()
    mafSortPath = tools.last.MafSort().install_and_get_path()
    mafConvertPath = tools.last.MafConvert().install_and_get_path()
    prinseqPath = tools.prinseq.PrinseqTool().install_and_get_path()
    noBlastLikeHitsPath = os.path.join( util.file.get_scripts_path(),
                                        'noBlastLikeHits.py')
    
    lastalOut = mkstempfname('.lastal')
    with open(lastalOut, 'wt') as outf:
        cmd = [lastalPath, '-Q1', db, inFastq]
        log.debug(' '.join(cmd) + ' > ' + lastalOut)
        subprocess.check_call(cmd, stdout=outf)
    # everything below this point in this method should be replaced with
    # our own code that just reads lastal output and makes a list of read names
    
    mafSortOut = mkstempfname('.mafsort')
    with open(mafSortOut, 'wt') as outf:
        with open(lastalOut, 'rt') as inf:
            cmd = [mafSortPath, '-n2']
            log.debug('cat ' + lastalOut + ' | ' + ' '.join(cmd) + ' > ' + mafSortOut)
            subprocess.check_call(cmd, stdin=inf, stdout=outf)
    os.unlink(lastalOut)
    
    mafConvertOut = mkstempfname('.mafconvert')
    with open(mafConvertOut, 'wt') as outf:
        cmd = [mafConvertPath, 'tab', mafSortOut]
        log.debug(' '.join(cmd) + ' > ' + mafConvertOut)
        subprocess.check_call(cmd, stdout=outf)
    os.unlink(mafSortOut)
    
    filteredFastq = mkstempfname('.filtered.fastq')
    with open(filteredFastq, 'wt') as outf:
        cmd = [noBlastLikeHitsPath, '-b', mafConvertOut, '-r', inFastq, '-m', 'hit']
        log.debug(' '.join(cmd) + ' > ' + filteredFastq)
        subprocess.check_call(cmd, stdout=outf)
    os.unlink(mafConvertOut)
    
    with open(outList, 'wt') as outf:
        with open(filteredFastq, 'rt') as inf:
            line_num = 0
            for line in inf:
                if (line_num % 4) == 0:
                    id = line.rstrip('\n\r')[1:]
                    if id.endswith('/1') or id.endswith('/2'):
                        id = id[:-2]
                    outf.write(id+'\n')
                line_num += 1

def filter_lastal_bam(inBam, db, outBam, JVMmemory=None):
    ''' Restrict input reads to those that align to the given
        reference database using LASTAL.
    '''
    # convert BAM to paired FASTQ
    inReads1 = mkstempfname('.1.fastq')
    inReads2 = mkstempfname('.2.fastq')
    tools.picard.SamToFastqTool().execute(inBam, inReads1, inReads2)
    
    # look for hits in inReads1 and inReads2
    hitList1 = mkstempfname('.1.hits')
    hitList2 = mkstempfname('.2.hits')
    lastal_get_hits(inReads1, db, hitList1)
    os.unlink(inReads1)
    lastal_get_hits(inReads2, db, hitList2)
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
    parser.add_argument("inBam",  help="Input reads")
    parser.add_argument("db",  help="Database of taxa we keep")
    parser.add_argument("outBam", help="Output reads, filtered to refDb")
    parser.add_argument('--JVMmemory', default = tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, filter_lastal_bam, split_args=True)
    return parser
__commands__.append(('filter_lastal_bam', parser_filter_lastal_bam))


def filter_lastal(inFastq, refDb, outFastq):
    ''' Restrict input reads to those that align to the given
        reference database using LASTAL.  Also, remove duplicates with prinseq.
    '''
    assert outFastq.endswith('.fastq')
    tempFilePath = mkstempfname('.hits')
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
        '{lastalPath} -Q1 {refDb} {inFastq}'.format(**locals()),
        '| {mafSortPath} -n2'.format(**locals()),
        '| {mafConvertPath} tab /dev/stdin > {tempFilePath}'.format(**locals()),
        ])
    log.debug(lastalCmd)
    assert not os.system(lastalCmd)

    # filter inFastq against lastal hits
    filteredFastq = mkstempfname('.filtered.fastq')
    with open(filteredFastq, 'wt') as outf:
        noBlastLikeHitsCmd = [
            noBlastLikeHitsPath, '-b', tempFilePath, '-r', inFastq, '-m', 'hit']
        log.debug(' '.join(noBlastLikeHitsCmd) + ' > ' + filteredFastq)
        subprocess.check_call(noBlastLikeHitsCmd, stdout=outf)
    
    # remove duplicate reads and reads with multiple Ns
    if os.path.getsize(filteredFastq) == 0:
        # prinseq-lite fails on empty file input (which can happen in real life
        # if no reads match the refDb) so handle this scenario specially
        log.info("output is empty: no reads in input match refDb")
        shutil.copyfile(filteredFastq, outFastq)
    else:
        prinseqCmd = [
            'perl', prinseqPath,
                '-ns_max_n', '1',
                '-derep', '1',
                '-fastq', filteredFastq,
                '-out_bad', 'null',
                '-line_width', '0',
                '-out_good', outFastq[:-6]
            ]
        log.debug(' '.join(prinseqCmd))
        subprocess.check_call(prinseqCmd)
    os.unlink(filteredFastq)

def parser_filter_lastal(parser=argparse.ArgumentParser()):
    parser.add_argument("inFastq", help="Input fastq file")
    parser.add_argument("refDb",
                        help="Reference database to retain from input")
    parser.add_argument("outFastq", help = "Output fastq file")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, filter_lastal, split_args=True)
    return parser
__commands__.append(('filter_lastal', parser_filter_lastal))


# ============================
# ***  partition_bmtagger  ***
# ============================

def deplete_bmtagger_bam(inBam, db, outBam, JVMmemory=None) :
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
    tools.picard.SamToFastqTool().execute(inBam, inReads1, inReads2)
    
    tempDir = tempfile.mkdtemp()
    matchesFile = mkstempfname('.txt')
    cmdline = [bmtaggerPath,
               '-b', db+'.bitmask', '-x', db+'.srprism', '-T', tempDir,
               '-q1', '-1', inReads1, '-2', inReads2,
               '-o', matchesFile]
    log.debug(' '.join(cmdline))
    subprocess.check_call(cmdline)
    
    tools.picard.FilterSamReadsTool().execute(inBam, True, matchesFile, outBam, JVMmemory=JVMmemory)
    os.unlink(matchesFile)

def select_reads(inFastq, outFastq, selectorFcn) :
    """
        selectorFcn: Bio.SeqRecord.SeqRecord -> bool
        Output in outFastq all reads from inFastq for which
            selectorFcn returns True.
        TO DO: change this to use Picard FilterSamReads (and operate
            on BAM files) which is likely much faster. This is the
            slowest step of partition_bmtagger currently.
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
                   '-b', db+'.bitmask', '-x', db+'.srprism', '-T', tempDir,
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
                   '-b', db+'.bitmask', '-x', db+'.srprism', '-T', tempDir,
                   '-q1', '-1', curReads1, '-2', curReads2,
                   '-o', outprefix]
        log.debug(' '.join(cmdline))
        subprocess.check_call(cmdline)
        curReads1, curReads2 = [outprefix+suffix for suffix in ('_1.fastq', '_2.fastq')]
        tempfiles += [curReads1, curReads2]
    shutil.copyfile(curReads1, outFastq1)
    shutil.copyfile(curReads2, outFastq2)
    for fn in tempfile:
        os.unlink(fn)
    log.debug("deplete_bmtagger complete")

def parser_partition_bmtagger(parser=argparse.ArgumentParser()):
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
    util.cmd.attach_main(parser, main_partition_bmtagger)
    return parser
def main_partition_bmtagger(args) :
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
    #if outMatch==None:
    #    deplete_bmtagger(inFastq1, inFastq2, databases, outNoMatch[0], outNoMatch[1])
    #else:
    #    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)
    #return 0
    partition_bmtagger(inFastq1, inFastq2, databases, outMatch, outNoMatch)
__commands__.append(('partition_bmtagger', parser_partition_bmtagger))

def parser_deplete_bam_bmtagger(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam',
        help='Input BAM file.')
    parser.add_argument('refDbs', nargs='+',
        help='''Reference databases (one or more) to deplete from input.
                For each db, requires prior creation of db.bitmask by bmtool,
                and db.srprism.idx, db.srprism.map, etc. by srprism mkindex.''')
    parser.add_argument('outBam',
        help='Output BAM file.')
    parser.add_argument('--JVMmemory', default = tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_deplete_bam_bmtagger)
    return parser
def main_deplete_bam_bmtagger(args) :
    '''Use bmtagger to deplete input reads against several databases.'''
    multi_db_deplete_bam(args.inBam, args.refDbs, deplete_bmtagger_bam, args.outBam, JVMmemory=args.JVMmemory)
__commands__.append(('deplete_bam_bmtagger', parser_deplete_bam_bmtagger))


def multi_db_deplete_bam(inBam, refDbs, deplete_method, outBam, JVMmemory=None):
    tmpBamIn = inBam
    for db in refDbs:
        tmpBamOut = mkstempfname('.bam')
        deplete_method(tmpBamIn, db, tmpBamOut, JVMmemory=JVMmemory)
        if tmpBamIn != inBam:
            os.unlink(tmpBamIn)
        tmpBamIn = tmpBamOut
    shutil.copyfile(tmpBamIn, outBam)


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
    inFasta = mkstempfname('.fasta')
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

def parser_deplete_blastn(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq',
        help='Input fastq file.')
    parser.add_argument('outFastq',
        help='Output fastq file with matching reads removed.')
    parser.add_argument('refDbs', nargs='+',
        help='One or more reference databases for blast.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, deplete_blastn, split_args=True)
    return parser
__commands__.append(('deplete_blastn', parser_deplete_blastn))

def deplete_blastn_paired(infq1, infq2, outfq1, outfq2, refDbs):
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
    deplete_blastn(tmpfq2_b, tmpfq2_c, refDbs)
    # purge fq1 of read pairs lost in fq2
    read_utils.purge_unmated(tmpfq1_b, tmpfq2_c, outfq1, outfq2)

def parser_deplete_blastn_paired(parser=argparse.ArgumentParser()):
    parser.add_argument('infq1',
        help='Input fastq file.')
    parser.add_argument('infq2',
        help='Input fastq file.')
    parser.add_argument('outfq1',
        help='Output fastq file with matching reads removed.')
    parser.add_argument('outfq2',
        help='Output fastq file with matching reads removed.')
    parser.add_argument('refDbs', nargs='+',
        help='One or more reference databases for blast.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, deplete_blastn_paired, split_args=True)
    return parser
__commands__.append(('deplete_blastn_paired', parser_deplete_blastn_paired))


def deplete_blastn_bam(inBam, db, outBam, JVMmemory=None):
    'Use blastn to remove reads that match at least one of the databases.'
    
    blastnPath = tools.blast.BlastnTool().install_and_get_path()
    fastq1 = mkstempfname('.1.fastq')
    fastq2 = mkstempfname('.2.fastq')
    fasta = mkstempfname('.1.fasta')
    blast_hits = mkstempfname('.blast_hits.txt')
    halfBam = mkstempfname('.half.bam')
    blastOutFile = mkstempfname('.hits.txt')
    
    # Initial BAM -> FASTQ pair
    tools.picard.SamToFastqTool().execute(inBam, fastq1, fastq2)
    
    # Find BLAST hits against FASTQ1
    read_utils.fastq_to_fasta(fastq1, fasta)
    os.unlink(fastq1)
    os.unlink(fastq2)
    log.info("running blastn on {} pair 1 against {}".format(inBam, db))
    blastnCmd = [blastnPath, '-db', db,
                '-word_size', '16', '-evalue', '1e-6', '-outfmt', '6',
                '-num_descriptions', '2', '-num_alignments', '2',
                '-query', fasta, '-out', blastOutFile]
    log.debug(' '.join(blastnCmd))
    subprocess.check_call(blastnCmd)
    with open(blast_hits, 'wt') as outf:
        with open(blastOutFile, 'rt') as inf:
            for line in inf:
                id = line.split('\t')[0].strip()
                if id.endswith('/1') or id.endswith('/2'):
                    id = id[:-2]
                outf.write(id+'\n')
        os.unlink(blastOutFile)
    
    # Deplete BAM of hits in FASTQ1
    tools.picard.FilterSamReadsTool().execute(inBam, True, blast_hits, halfBam, JVMmemory=JVMmemory)
    
    # Depleted BAM -> FASTQ pair
    tools.picard.SamToFastqTool().execute(halfBam, fastq1, fastq2)
    
    # Find BLAST hits against FASTQ2 (which is already smaller than before)
    read_utils.fastq_to_fasta(fastq2, fasta)
    os.unlink(fastq1)
    os.unlink(fastq2)
    log.info("running blastn on {} pair 2 against {}".format(inBam, db))
    blastnCmd = [blastnPath, '-db', db,
                '-word_size', '16', '-evalue', '1e-6', '-outfmt', '6',
                '-num_descriptions', '2', '-num_alignments', '2',
                '-query', fasta, '-out', blastOutFile]
    log.debug(' '.join(blastnCmd))
    subprocess.check_call(blastnCmd)
    with open(blast_hits, 'wt') as outf:
        with open(blastOutFile, 'rt') as inf:
            for line in inf:
                id = line.split('\t')[0].strip()
                if id.endswith('/1') or id.endswith('/2'):
                    id = id[:-2]
                outf.write(id+'\n')
        os.unlink(blastOutFile)
    
    # Deplete BAM of hits against FASTQ2
    tools.picard.FilterSamReadsTool().execute(halfBam, True, blast_hits, outBam, JVMmemory=JVMmemory)
    
    # Clean up
    for fn in (fasta, blast_hits, halfBam):
        os.unlink(fn)
    
def parser_deplete_blastn_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam',
        help='Input BAM file.')
    parser.add_argument('refDbs', nargs='+',
        help='One or more reference databases for blast.')
    parser.add_argument('outBam',
        help='Output BAM file with matching reads removed.')
    parser.add_argument('--JVMmemory', default = tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_deplete_blastn_bam)
    return parser
def main_deplete_blastn_bam(args):
    '''Use blastn to remove reads that match at least one of the specified databases.'''
    multi_db_deplete_bam(args.inBam, args.refDbs, deplete_blastn_bam, args.outBam, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('deplete_blastn_bam', parser_deplete_blastn_bam))
                     
# ========================

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
