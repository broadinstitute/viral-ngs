#!/usr/bin/env python
'''This script contains a number of utilities for filtering NGS reads based
on membership or non-membership in a species / genus / taxonomic grouping.
'''

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
                + "hlevitin@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, os, tempfile, errno
from Bio import SeqIO
import util.cmd, util.file, util.vcf, util.misc
import tools.last, tools.prinseq, tools.trimmomatic, tools.bmtagger, \
       tools.samtools, tools.picard, tools.blast, tools.mvicuna
from util.file import mkstempfname

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

    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
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
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))

    # Future: handle BAM input and output; handle multiple databases.
    # Will need to implement bam->fastq->bam wrappers that maintain the read
    # metadata
    #parser.add_argument("inBam", help="Input BAM file")
    #parser.add_argument("outBam", help="Output BAM file")

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
        matches = set(line.strip() for line in open(matchesFile))
        noMatchFcn = lambda rec : rec.id not in matches
        select_reads(prevReads1, curReads1, noMatchFcn)
        select_reads(prevReads2, curReads2, noMatchFcn)
    if outMatch != None :
        allMatches = set(line.strip()
                         for matchesFile in matchesFiles
                         for line in open(matchesFile))
        matchFcn = lambda rec : rec.id in allMatches
        select_reads(inFastq1, outMatch[0], matchFcn)
        select_reads(inFastq2, outMatch[1], matchFcn)

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
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
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

def dup_remove_mvicuna(inPair, pairedOutPair, unpairedOut, postDupRmPair) :
    """
    Run mvicuna's duplicate removal operation on paired-end input reads in
        fastq format, producing various outputs in fastq format.
    inPair, pairedOutPair, and postDupRmOutPair are pairs of file names,
        while unpairedOut is a single file name.
    """
    mvicunaPath = tools.mvicuna.MvicunaTool().install_and_get_path()
    input       = ','.join(inPair)
    pairedOut    = ','.join(pairedOutPair)
    postDupRmOut = ','.join(postDupRmPair)
    cmdline = ('{mvicunaPath} '
               '-ipfq {input} '
               '-opfq {pairedOut} '
               '-osfq {unpairedOut} '
               '-drm_op {postDupRmOut} '
               '-tasks DupRm').format(**locals())
    log.debug(cmdline)
    assert not os.system(cmdline)

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
        help='File name of output unpaired reads')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    return parser
def main_dup_remove_mvicuna(args) :
    inPair = (args.inFastq1, args.inFastq2)
    pairedOutPair = (args.pairedOutFastq1, args.pairedOutFastq2)
    unpairedOutFastq = args.unpairedOutFastq
    if unpairedOutFastq == None :
        unpairedOutFastq = mkstempfname()
    postDupRmPair = (mkstempfname(), mkstempfname())
    dup_remove_mvicuna(inPair, pairedOutPair, unpairedOutFastq, postDupRmPair)
    return 0

__commands__.append(('dup_remove_mvicuna', main_dup_remove_mvicuna,
                     parser_dup_remove_mvicuna))


# ========================
# ***  deplete_blastn  ***
# ========================

def deplete_blastn(inFastq, outFastq, refDbs) :
    'Use blastn to remove reads that match at least one of the databases.'
    
    ## Get tools
    prinseqLitePath = tools.prinseq.PrinseqTool().install_and_get_path()
    blastnPath = tools.blast.BlastnTool().install_and_get_path()
    noBlastHits_v3Path = os.path.join(util.file.get_scripts_path(),
                                      'noBlastHits_v3.py')
    
    ## Convert to fasta
    inFastaBase = mkstempfname()
    inFasta = inFastaBase + '.fasta' # prinseqLite adds the .fasta
    prinseqCmd = '{prinseqLitePath} -out_format 1 -line_width 0 '          \
                 '-fastq {inFastq} -out_good {inFastaBase} -out_bad null'. \
                 format(**locals())
    log.debug(prinseqCmd)
    assert not os.system(prinseqCmd)

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
    return parser
def main_deplete_blastn(args) :
    inFastq = args.inFastq
    outFastq = args.outFastq
    refDbs = args.refDbs
    deplete_blastn(inFastq, outFastq, refDbs)
    return 0
__commands__.append(('deplete_blastn', main_deplete_blastn,
                     parser_deplete_blastn))


# =======================
# ***  purge_unmated  ***
# =======================

def purge_unmated(inFastq1, inFastq2, outFastq1, outFastq2) :
    """Use mergeShuffledFastqSeqs to purge unmated reads, and put corresponding
       reads in the same order."""
    tempOutput = mkstempfname()
    mergeShuffledFastqSeqsPath = os.path.join(util.file.get_scripts_path(),
                                              'mergeShuffledFastqSeqs.pl')
    # The regular expression that follow says that the sequence identifiers
    # of corresponding sequences must be of the form SEQID/1 and SEQID/2
    cmdline = "{mergeShuffledFastqSeqsPath} -t -r '^@(\S+)/[1|2]$' " \
              "-f1 {inFastq1} -f2 {inFastq2} -o {tempOutput}".format(**locals())
    log.debug(cmdline)
    assert not os.system(cmdline)
    os.rename(tempOutput + '.1.fastq', outFastq1)
    os.rename(tempOutput + '.2.fastq', outFastq2)

def parser_purge_unmated() :
    parser = argparse.ArgumentParser(
        description='''Use mergeShuffledFastqSeqs to purge unmated reads, and
                       put corresponding reads in the same order.
                       Corresponding sequences must have sequence identifiers
                       of the form SEQID/1 and SEQID/2.
                    ''')
    parser.add_argument('inFastq1',
        help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2',
        help='Input fastq file; 2nd end of paired-end reads.')
    parser.add_argument('outFastq1',
        help='Output fastq file; 1st end of paired-end reads.')
    parser.add_argument('outFastq2',
        help='Output fastq file; 2nd end of paired-end reads.')
    return parser

def main_purge_unmated(args) :
    inFastq1 = args.inFastq1
    inFastq2 = args.inFastq2
    outFastq1 = args.outFastq1
    outFastq2 = args.outFastq2
    purge_unmated(inFastq1, inFastq2, outFastq1, outFastq2)
    return 0

__commands__.append(('purge_unmated', main_purge_unmated,
                     parser_purge_unmated))


''' KGA's "recipe" for human read depletion (with some notes by Irwin)
###-------- CLEANING OF READS FOR SRA SUBMISSION --------#
# MAKE REQUIRED SUB-DIRECTORIES - DON'T DO THIS IF YOU ALREADY CREATED THESE WITH ANOTHER PIPELINE
for directory in
do
bsub -o ~/log.txt -P sabeti_meta "mkdir $directory/_logs $directory/_temp $directory/_bams $directory/_reports $directory/_pileup $directory/_meta $directory/_reads"
done

## BMTAGGER REMOVAL OF HUMAN READS AND CONTAMINANTS
# Rodent sequences can be removed using mm9_mn, nt_rodent.1 and nt_rodent.2
for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
for db1 in GRCh37.68_ncRNA-GRCh37.68_transcripts-HS_rRNA_mitRNA
do
for db2 in hg19
do
for db3 in metagenomics_contaminants_v3 # If you want to remove Lassa, use metagenomics_contaminants_v3_w_lassa, ZEBOV use metagenomics_contaminants_v3_w_zebov
do

# Tools and files needed
bmtaggerPath = /idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh
bitmaskDb1 = /idi/sabeti-scratch/kandersen/references/bmtagger/$db1.bitmask
srprismDb1 = /idi/sabeti-scratch/kandersen/references/bmtagger/$db1.srprism
bitmaskDb2 = /idi/sabeti-scratch/kandersen/references/bmtagger/$db2.bitmask
srprismDb2 = /idi/sabeti-scratch/kandersen/references/bmtagger/$db2.srprism
bitmaskDb3 = /idi/sabeti-scratch/kandersen/references/bmtagger/$db3.bitmask
srprismDb3 = /idi/sabeti-scratch/kandersen/references/bmtagger/$db3.srprism
inreads1 = $directory/_reads/$sample.reads1.fastq
inreads2 = $directory/_reads/$sample.reads2.fastq
tempMrna = $temp/$sample.bmtagger.mrna
tempMrna1 = $temp/$sample.bmtagger.mrna.1.fastq # Presumably produced by previous call to bmtagger
tempMrna2 = $temp/$sample.bmtagger.mrna.2.fastq
tempHg19 = $temp/$sample.bmtagger.hg19
tempHg19_1 = $temp/$sample.bmtagger.hg19.1.fastq
tempHg10_2 = $temp/$sample.bmtagger.hg19.2.fastq
tempContaminants = $temp/$sample.bmtagger.contaminants

# Commands to execute
bmtaggerPath -X -b bitmaskDb1 -x srprismDb1 -T $temp -q1 -1 inreads1 -2 inreads2 -o tempMrna && 
bmtaggerPath -X -b bitmaskDb2 -x srprismDb2 -T $temp -q1 -1 tempMrna1 -2 tempMrna2 -o tempHg19 && 
bmtaggerPath -X -b bitmaskDb3 -x srprismDb3 -T $temp -q1 -1 tempHg19_1 -2 tempHg19_2 -o tempContaminants


bsub -R "rusage[mem=8]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.bt "/idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh -X -b /idi/sabeti-scratch/kandersen/references/bmtagger/$db1.bitmask -x /idi/sabeti-scratch/kandersen/references/bmtagger/$db1.srprism -T $temp -q1 -1 $directory/_reads/$sample.reads1.fastq -2 $directory/_reads/$sample.reads2.fastq -o $temp/$sample.bmtagger.mrna && /idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh -X -b /idi/sabeti-scratch/kandersen/references/bmtagger/$db2.bitmask -x /idi/sabeti-scratch/kandersen/references/bmtagger/$db2.srprism -T $temp -q1 -1 $temp/$sample.bmtagger.mrna.1.fastq -2 $temp/$sample.bmtagger.mrna.2.fastq -o $temp/$sample.bmtagger.hg19 && /idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh -X -b /idi/sabeti-scratch/kandersen/references/bmtagger/$db3.bitmask -x /idi/sabeti-scratch/kandersen/references/bmtagger/$db3.srprism -T $temp -q1 -1 $temp/$sample.bmtagger.hg19.1.fastq -2 $temp/$sample.bmtagger.hg19.2.fastq -o $temp/$sample.bmtagger.contaminants"
done
done
done
done
done
done

## REMOVE DUPLICATES
for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do


# Tools and files needed
mvicunaPath = /gsap/garage-viral/viral/analysis/xyang/programs/M-Vicuna/bin/mvicuna
contaminants1fastq = $temp/$sample.bmtagger.contaminants.1.fastq
contaminants2fastq = $temp/$sample.bmtagger.contaminants.2.fastq
cleaned1fastq = $temp/$sample.cleaned_reads.prinseq.1.fastq
cleaned2fastq = $temp/$sample.cleaned_reads.prinseq.2.fastq
cleanedUnpairedfastq = $temp/$sample.cleaned_reads.unpaired.fastq
hg19temp1fastq = $temp/$sample.bmtagger.hg19.temp1.fastq
hg19temp2fastq = $temp/$sample.bmtagger.hg19.temp2.fastq
# old path: /seq/viral/analysis/xyang/programs/M-Vicuna/bin/mvicuna

# Commands to execute...
mvicunaPath -ipfq contaminants1fastq,contaminants2fastq -opfq cleaned1fastq,cleaned2fastq -osfq cleanedUnpairedfastq -drm_op hg19temp1fastq,hg19temp2fastq -tasks DupRm



bsub -R "rusage[mem=$memory]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.p1 "/seq/viral/analysis/xyang/programs/M-Vicuna/bin/mvicuna -ipfq $temp/$sample.bmtagger.contaminants.1.fastq,$temp/$sample.bmtagger.contaminants.2.fastq -opfq $temp/$sample.cleaned_reads.prinseq.1.fastq,$temp/$sample.cleaned_reads.prinseq.2.fastq -osfq $temp/$sample.cleaned_reads.unpaired.fastq -drm_op $temp/$sample.bmtagger.hg19.temp1.fastq,$temp/$sample.bmtagger.hg19.temp2.fastq -tasks DupRm"
done
done
done

## REMOVE HUMAN READS AND CONTAMINANTS USING NOVOALIGN
for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do

# NOTE: 10/23/14 Danny said this step could be eliminated

# Tools and files needed
novoalignPath = /idi/sabeti-scratch/kandersen/bin/novocraft/novoalign
sortSamPath = /seq/software/picard/current/bin/SortSam.jar
samToFastqPath = /seq/software/picard/current/bin/SamToFastq.jar
prinseqLitePath = /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl
prinseq1fastq = $temp/$sample.cleaned_reads.prinseq.1.fastq
prinseq2fastq = $temp/$sample.cleaned_reads.prinseq.2.fastq
novalignlog = $directory/_logs/$sample.log.novoalign.txt
sampleSortedBam = $directory/_temp/$sample.sorted.bam
sampleDepletedBam = $directory/_temp/$sample.depleted.bam
novoDepleted1Fastq = $directory/_temp/$sample.novo.depleted.reads1.fastq
novoDepleted2Fastq = $directory/_temp/$sample.novo.depleted.reads2.fastq
prinseq1 = $directory/_temp/$sample.prinseq.1
prinseq2 = $directory/_temp/$sample.prinseq.2
refnix = /idi/sabeti-scratch/kandersen/references/novo_clean/metag_v3.ncRNA.mRNA.mitRNA.consensus.nix

# Commands to execute
novoalignPath -c 4 -f prinseq1fastq prinseq2fastq -r Random -l 30 -g 20 -x 6 -t 502 -F STDFQ -d refnix -o SAM $'@RG\tID:140813.$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> novalignlog |
java -Xmx2g -jar sortSamPath SO=coordinate I=/dev/stdin O=sampleSortedBam CREATE_INDEX=true &&
samtools view -b -f 4 -u sampleSortedBam > sampleDepletedBam &&
java -Xmx2g -jar samToFastqPath INPUT=sampleDepletedBam FASTQ=novoDepleted1Fastq SECOND_END_FASTQ=novoDepleted2Fastq VALIDATION_STRINGENCY=SILENT &&
prinseqLitePath -out_format 1 -line_width 0 -fastq novoDepleted1Fastq -out_good prinseq1 -out_bad null &&
prinseqLitePath -out_format 1 -line_width 0 -fastq novoDepleted2Fastq -out_good prinseq2 -out_bad null


bsub -W 4:00 -q hour -R "rusage[mem=4]" -n 4 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.al "/idi/sabeti-scratch/kandersen/bin/novocraft/novoalign -c 4 -f $temp/$sample.cleaned_reads.prinseq.1.fastq $temp/$sample.cleaned_reads.prinseq.2.fastq -r Random -l 30 -g 20 -x 6 -t 502 -F STDFQ -d /idi/sabeti-scratch/kandersen/references/novo_clean/metag_v3.ncRNA.mRNA.mitRNA.consensus.nix -o SAM $'@RG\tID:140813.$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> $directory/_logs/$sample.log.novoalign.txt | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.sorted.bam CREATE_INDEX=true && samtools view -b -f 4 -u $directory/_temp/$sample.sorted.bam > $directory/_temp/$sample.depleted.bam && java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$directory/_temp/$sample.depleted.bam FASTQ=$directory/_temp/$sample.novo.depleted.reads1.fastq SECOND_END_FASTQ=$directory/_temp/$sample.novo.depleted.reads2.fastq VALIDATION_STRINGENCY=SILENT && /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -out_format 1 -line_width 0 -fastq $directory/_temp/$sample.novo.depleted.reads1.fastq -out_good $directory/_temp/$sample.prinseq.1 -out_bad null && /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -out_format 1 -line_width 0 -fastq $directory/_temp/$sample.novo.depleted.reads2.fastq -out_good $directory/_temp/$sample.prinseq.2 -out_bad null"
done
done
done

## SPLIT FILES FOR BLASTN ANALYSIS


# Tools and files needed
prinseq1Fasta = $directory/_temp/$sample.prinseq.1.fasta
prinseq2Fasta = $directory/_temp/$sample.prinseq.2.fasta
prinseq1Split = $temp/$sample.prinseq.1.split.
prinseq1Split = $temp/$sample.prinseq.2.split.

# Commands to execute
split -a 3 -l 20000 prinseq1Fasta prinseq1Split
split -a 3 -l 20000 prinseq2Fasta prinseq2Split

for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
bsub -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.sf "split -a 3 -l 20000 $directory/_temp/$sample.prinseq.1.fasta $temp/$sample.prinseq.1.split."
bsub -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.sf "split -a 3 -l 20000 $directory/_temp/$sample.prinseq.2.fasta $temp/$sample.prinseq.2.split."
done
done
done

## RUN BLASTN ANALYSIS

# Tools and files needed
referenceDb = /idi/sabeti-scratch/kandersen/references/blast/metag_v3.ncRNA.mRNA.mitRNA.consensus
prinseq1Split = $temp/$sample.prinseq.1.split.XXX
prinseq2Split = $temp/$sample.prinseq.2.split.XXX

sample1XXX = $temp/$sample.1.$db.$((i++)).txt
sample2XXX = $temp/$sample.2.$db.$((i++)).txt

# Commands to execute
blastn -db referenceDb -word_size 16 -evalue 1e-6 -outfmt 6 -num_descriptions 2 -num_alignments 2 -query prinseq1Split -out sample1XXX
blastn -db referenceDb -word_size 16 -evalue 1e-6 -outfmt 6 -num_descriptions 2 -num_alignments 2 -query prinseq2Split -out sample2XXX


for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
for db in metag_v3.ncRNA.mRNA.mitRNA.consensus
do
i=1
j=1
for a in $temp/$sample.prinseq.1.split.*
do
bsub -R "rusage[mem=2]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.$((j++)).bn "blastn -db /idi/sabeti-scratch/kandersen/references/blast/$db -word_size 16 -evalue 1e-6 -outfmt 6 -num_descriptions 2 -num_alignments 2 -query $a -out $temp/$sample.1.$db.$((i++)).txt"
done
done
done
done
done
for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
for db in metag_v3.ncRNA.mRNA.mitRNA.consensus
do
i=1
j=1
for b in $temp/$sample.prinseq.2.split.*
do
bsub -R "rusage[mem=$memory]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.$((j++)).bn "blastn -db /idi/sabeti-scratch/kandersen/references/blast/$db -word_size 16 -evalue 1e-6 -outfmt 6 -num_descriptions 2 -num_alignments 2 -query $b -out $temp/$sample.2.$db.$((i++)).txt"
done
done
done
done
done

## CONCATENATE BLASTN RESULTS
for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
bsub -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.c1 "cat $temp/$sample.1.*.*.txt > $directory/_temp/$sample.blast.1.txt"
bsub -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.c2 "cat $temp/$sample.2.*.*.txt > $directory/_temp/$sample.blast.2.txt"
done
done
done

## EXTRACT READS WITH NO BLAST HITS

# Tools and files needed
noBlastHits_v3Path = /idi/sabeti-scratch/kandersen/bin/scripts/noBlastHits_v3.py
blastOutput1Txt = $directory/_temp/$sample.blast.1.txt
in1Fastq = $directory/_temp/$sample.novo.depleted.reads1.fastq
out1Fastq = $temp/$sample.nohits.1.fastq


# Commands to execute
python noBlastHits_v3Path -b blastOutput1Txt -r in1Fastq -m nohit > out1Fastq
Same with 2.

for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
bsub -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.nb "python /idi/sabeti-scratch/kandersen/bin/scripts/noBlastHits_v3.py -b $directory/_temp/$sample.blast.1.txt -r $directory/_temp/$sample.novo.depleted.reads1.fastq -m nohit > $temp/$sample.nohits.1.fastq"
bsub -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.nb "python /idi/sabeti-scratch/kandersen/bin/scripts/noBlastHits_v3.py -b $directory/_temp/$sample.blast.2.txt -r $directory/_temp/$sample.novo.depleted.reads2.fastq -m nohit > $temp/$sample.nohits.2.fastq"
done
done
done

## FIX MATE-PAIR INFORMATION

# Tools and files needed
mergeShuffledFastqSeqsPath = /idi/sabeti-scratch/kandersen/bin/scripts/mergeShuffledFastqSeqs.pl
in1Fastq = $temp/$sample.nohits.1.fastq
in2Fastq = $temp/$sample.nohits.2.fastq
outFastq = $directory/_temp/$sample.cleaned (makes ${outFastq}.[12].fastq)

# Commands to execute
mergeShuffledFastqSeqsPath -t -r '^@(\S+)/[1|2]$' -f1 in1Fastq -f2 in2Fastq -o outFastq

for sample in
do
for directory in
do
for temp in /broad/hptmp/andersen
do
bsub -R "rusage[mem=4]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.fm "/idi/sabeti-scratch/kandersen/bin/scripts/mergeShuffledFastqSeqs.pl -t -r '^@(\S+)/[1|2]$' -f1 $temp/$sample.nohits.1.fastq -f2 $temp/$sample.nohits.2.fastq -o $directory/_temp/$sample.cleaned"
done
done
done

## GET NON-VIRAL AND VIRAL READS
for sample in
do
for directory in
do
for reference in $directory/_refs/zaire_guinea.nix
do

# Tools and files needed
novoalign3Path = /idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign

# Inputs 1st call:
$directory/_reads/$sample.reads[12].fastq
$reference

# Outputs 1st call:
$directory/_temp/$sample.viral.reads[12].fastq

# Inputs 2nd call:
$directory/_temp/$sample.cleaned.[12].fastq
$reference

# Outputs 2nd call:
$directory/_temp/$sample.viral.depleted.reads[12].fastq


# Commands to execute
novoalign3Path -c 1 -f $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq -r Random -l 40 -g 20 -x 6 -t 502 -F STDFQ -d $reference -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> $directory/_logs/$sample.log.viral.txt |
java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.sorted3.bam CREATE_INDEX=true &&
samtools view -b -q 1 -u $directory/_temp/$sample.sorted3.bam |
java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.mapped3.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT &&
java -Xmx2g -jar /seq/software/picard/current/bin/MarkDuplicates.jar I=$directory/_temp/$sample.mapped3.bam O=$directory/_temp/$sample.mappedNoDub3.bam METRICS_FILE=$directory/_temp/$sample.log.markdups3.txt CREATE_INDEX=true REMOVE_DUPLICATES=true &&
java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$directory/_temp/$sample.mappedNoDub3.bam FASTQ=$directory/_temp/$sample.viral.reads1.fastq SECOND_END_FASTQ=$directory/_temp/$sample.viral.reads2.fastq VALIDATION_STRINGENCY=SILENT

novoalign3Path -c 1 -f $directory/_temp/$sample.cleaned.1.fastq $directory/_temp/$sample.cleaned.2.fastq -r Random -l 40 -g 20 -x 6 -t 502 -F STDFQ -d $reference -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> $directory/_logs/$sample.log.viral-deplete.txt |
java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.sorted2.bam CREATE_INDEX=true &&
samtools view -b -f 4 -u $directory/_temp/$sample.sorted2.bam > $directory/_temp/$sample.depleted2.bam &&
java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$directory/_temp/$sample.depleted2.bam FASTQ=$directory/_temp/$sample.viral.depleted.reads1.fastq SECOND_END_FASTQ=$directory/_temp/$sample.viral.depleted.reads2.fastq VALIDATION_STRINGENCY=SILENT





bsub -q week -W 24:00 -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.a2 "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -c 1 -f $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq -r Random -l 40 -g 20 -x 6 -t 502 -F STDFQ -d $reference -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> $directory/_logs/$sample.log.viral.txt | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.sorted3.bam CREATE_INDEX=true && samtools view -b -q 1 -u $directory/_temp/$sample.sorted3.bam | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.mapped3.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && java -Xmx2g -jar /seq/software/picard/current/bin/MarkDuplicates.jar I=$directory/_temp/$sample.mapped3.bam O=$directory/_temp/$sample.mappedNoDub3.bam METRICS_FILE=$directory/_temp/$sample.log.markdups3.txt CREATE_INDEX=true REMOVE_DUPLICATES=true && java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$directory/_temp/$sample.mappedNoDub3.bam FASTQ=$directory/_temp/$sample.viral.reads1.fastq SECOND_END_FASTQ=$directory/_temp/$sample.viral.reads2.fastq VALIDATION_STRINGENCY=SILENT"
bsub -q hour -W 4:00 -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.a1 "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -c 1 -f $directory/_temp/$sample.cleaned.1.fastq $directory/_temp/$sample.cleaned.2.fastq -r Random -l 40 -g 20 -x 6 -t 502 -F STDFQ -d $reference -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> $directory/_logs/$sample.log.viral-deplete.txt | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.sorted2.bam CREATE_INDEX=true && samtools view -b -f 4 -u $directory/_temp/$sample.sorted2.bam > $directory/_temp/$sample.depleted2.bam && java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$directory/_temp/$sample.depleted2.bam FASTQ=$directory/_temp/$sample.viral.depleted.reads1.fastq SECOND_END_FASTQ=$directory/_temp/$sample.viral.depleted.reads2.fastq VALIDATION_STRINGENCY=SILENT"
done
done
done

## COMBINE READS
for sample in
do
for directory in
do
bsub -R "rusage[mem=4]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.fm "cat $directory/_temp/$sample.viral.depleted.reads1.fastq $directory/_temp/$sample.viral.reads1.fastq > $directory/_reads/$sample.cleaned.1.fastq && cat $directory/_temp/$sample.viral.depleted.reads2.fastq $directory/_temp/$sample.viral.reads2.fastq > $directory/_reads/$sample.cleaned.2.fastq"
done
done

## CONVERT TO BAM FILE

# Tools and files needed
FastqToSamTool = /seq/software/picard/current/bin/FastqToSam.jar
in1Fastq, in2Fastq = $directory/_reads/$sample.cleaned.[12].fastq
outBam = $directory/_bams/$sample.bam
sampleName = $sample

# Commands to execute
java -Xmx2g -jar FastqToSamTool FASTQ=in1Fastq FASTQ2=in2Fastq OUTPUT=outBam SAMPLE_NAME=sampleName LIBRARY_NAME=sampleName PLATFORM=illumina SEQUENCING_CENTER=broad RUN_DATE=$date CREATE_MD5_FILE=True



for sample in
do
for directory in
do
for date in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.BF "java -Xmx2g -jar /seq/software/picard/current/bin/FastqToSam.jar FASTQ=$directory/_reads/$sample.cleaned.1.fastq FASTQ2=$directory/_reads/$sample.cleaned.2.fastq OUTPUT=$directory/_bams/$sample.bam SAMPLE_NAME=$sample LIBRARY_NAME=$sample PLATFORM=illumina SEQUENCING_CENTER=broad RUN_DATE=$date CREATE_MD5_FILE=True"
done
done
done
    '''
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
