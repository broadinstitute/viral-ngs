#!/usr/bin/env python
"""
Utilities for working with sequence reads, such as converting formats and
fixing mate pairs.
"""

__author__ = "irwin@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, os
from Bio import SeqIO
import util.cmd, util.file
from util.file import mkstempfname
import tools.picard, tools.samtools

log = logging.getLogger(__name__)


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

# =========================
# ***  fastq_to_fasta   ***
# =========================
def fastq_to_fasta(inFastq, outFasta) :
    'Convert from fastq format to fasta format.'
    'Warning: output reads might be split onto multiple lines.'
    
    # Do this with biopython rather than prinseq, because if the latter fails
    #    it doesn't return an error. (On the other hand, prinseq
    #    can guarantee that output lines are not split...)
    inFile  = util.file.open_or_gzopen(inFastq)
    outFile = util.file.open_or_gzopen(outFasta, 'w')
    for rec in SeqIO.parse(inFile, 'fastq') :
        SeqIO.write([rec], outFile, 'fasta')
    outFile.close()

def parser_fastq_to_fasta() :
    parser = argparse.ArgumentParser(
        description='Convert from fastq format to fasta format.')
    parser.add_argument('inFastq', help='Input fastq file.')
    parser.add_argument('outFasta', help='Output fasta file.')
    return parser

def main_fastq_to_fasta(args) :
    inFastq = args.inFastq
    outFasta = args.outFasta
    fastq_to_fasta(inFastq, outFasta)
    return 0
__commands__.append(('fastq_to_fasta', main_fastq_to_fasta,
                     parser_fastq_to_fasta))

# =======================
# ***  bam_to_fastq   ***
# =======================
# TBD...


# =======================
# ***  fastq_to_bam   ***
# =======================

def fastq_to_bam(sampleName, inFastq1, inFastq2, outBam,
                 header = None, JVMmemory = None, libraryName = None,
                 runDate = None, platform = None, sequencingCenter = None) :
    'Convert a pair of fastq paired-end read files and optional text header ' \
    'to a single bam file.'
    
    FastqToSamTool = tools.picard.FastqToSamTool().install_and_get_path()
    if header :
        fastqToSamOut = mkstempfname()
    else :
        fastqToSamOut = outBam
    memStr = JVMmemory if JVMmemory else '2g'
    options = []
    if libraryName :
        options.append('LIBRARY_NAME=' + libraryName)
    if runDate :
        options.append('RUN_DATE=' + runDate)
    if platform :
        options.append('PLATFORM=' + platform)
    if sequencingCenter :
        options.append('SEQUENCING_CENTER=' + sequencingCenter)
    fastqToSamCmd = ('java -Xmx{memStr} -jar {FastqToSamTool} '
                     'FASTQ={inFastq1} FASTQ2={inFastq2} '
                     'OUTPUT={fastqToSamOut} '
                     'SAMPLE_NAME={sampleName} '
                     'CREATE_MD5_FILE=True ' +
                     ' '.join(options)).format(**locals())
    log.debug(fastqToSamCmd)
    assert not os.system(fastqToSamCmd)
    if header :
        samtoolsPath = tools.samtools.SamtoolsTool().install_and_get_path()
        reheadCmd = '{samtoolsPath} reheader {header} {fastqToSamOut} >' \
                    '{outBam}'.format(**locals())
        log.debug(reheadCmd)
        assert not os.system(reheadCmd)
        os.system('md5 {outBam} > {outBam}.md5'.format(**locals()))

def parser_fastq_to_bam() :
    defJvmMem = '2g'
    parser = argparse.ArgumentParser(
        description='Convert a pair of fastq paired-end read files and '
                    'optional text header to a single bam file.')
    parser.add_argument('sampleName',
        help='Sample name to insert into the read group header.')
    parser.add_argument('inFastq1',
        help='Input fastq file; 1st end of paired-end reads.')
    parser.add_argument('inFastq2',
        help='Input fastq file; 2nd end of paired-end reads.')
    parser.add_argument('outBam', help='Output bam file.')
    parser.add_argument('--header',
        help='Optional text file containing header.')
    parser.add_argument('--JVMmemory', default = defJvmMem,
        help='JVM virtual memory size (default: {})'.format(defJvmMem))
    parser.add_argument('--libraryName',
        help='Library name to put in LB attribute in read group header.')
    parser.add_argument('--runDate',
        help='Run date in Iso8601 format (e.g., 2014-10-29).')
    parser.add_argument('--platform',
        help='The platform type (e.g. illumina, solid).')
    parser.add_argument('--sequencingCenter',
        help='The sequencing center from which the data originated.')
    return parser

def main_fastq_to_bam(args) :
    sampleName = args.sampleName
    inFastq1 = args.inFastq1
    inFastq2 = args.inFastq2
    outBam = args.outBam
    header = args.header
    JVMmemory = args.JVMmemory
    libraryName = args.libraryName
    runDate = args.runDate
    platform = args.platform
    sequencingCenter = args.sequencingCenter
    fastq_to_bam(sampleName, inFastq1, inFastq2, outBam, header,
                 JVMmemory, libraryName, runDate, platform, sequencingCenter)
    return 0
__commands__.append(('fastq_to_bam', main_fastq_to_bam,
                     parser_fastq_to_bam))




# =======================

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
