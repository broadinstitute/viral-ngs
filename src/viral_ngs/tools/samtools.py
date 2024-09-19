'''
    The Samtools package.

    TO DO: much of this stuff can be all eliminated by using pysam instead, as
    pysam (in its current versions) is meant to be the complete python
    implementation of htslib/samtools.

    http://pysam.readthedocs.org/en/latest/usage.html#using-samtools-commands-within-python

    Current bug with pysam 0.8.1: nosetests does not work unless you use --nocapture.
    python -m unittest works. Something about redirecting stdout.
    Actually, Travis CI still has issues with pysam and stdout even with --nocapture.

    pysam 0.9.1 stops redirecting stdout which makes things much easier,
    but it's a pain to pip install.
'''

import logging
import shutil
import os
import re
import os.path
import subprocess
import tempfile
import contextlib
from collections import OrderedDict
from decimal import *

import pysam

import tools
import util.file
import util.misc

TOOL_NAME = 'samtools'

log = logging.getLogger(__name__)

class SamtoolsTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(SamtoolsTool, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '--version']).decode('UTF-8').split('\n')[0].split()[1]

    def execute(self, command, args, stdout=None, stderr=None, background=False):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        if stderr:
            stderr = open(stderr, 'w')

        env = os.environ.copy()
        env.pop('JAVA_HOME', None)
        if background:
            subprocess.Popen(tool_cmd, stdout=stdout, stderr=stderr, env=env)
        else:
            subprocess.check_call(tool_cmd, stdout=stdout, stderr=stderr, env=env)

            if stdout:
                stdout.close()
            if stderr:
                stderr.close()

    def view(self, args, inFile, outFile, regions=None, threads=None, background=False):
        regions = regions or []
        args    = args or []

        # -@ seems to result in segfaults in some cases
        # so this is commented out for now...

        self.execute('view', args + ['-o', outFile, inFile] + regions, background=background)
        #opts = args + ['-o', outFile, inFile] + regions
        #pysam.view(*opts)

    def bam2fq(self, inBam, outFq1, outFq2=None):
        if outFq2 is None:
            self.execute('bam2fq', ['-n', inBam], stdout=outFq1)
        else:
            self.execute('bam2fq', ['-1', outFq1, '-2', outFq2, inBam])

    def bam2fq_pipe(self, inBam):
        tool_cmd = [self.install_and_get_path(), 'bam2fq', '-n', inBam]
        log.debug(' '.join(tool_cmd) + ' |')
        p = subprocess.Popen(tool_cmd, stdout=subprocess.PIPE)
        return p

    def bam2fa(self, inBam, outFa1, outFa2=None, outFa0=None, append_mate_num=False):
        if outFa2 is None:
            args = ['-N' if append_mate_num else '-n']
        else:
            args = ['-1', outFa1, '-2', outFa2]
        if outFa0:
            args += ['-0', outFa0]
        args += [inBam]
        if outFa2 is None:
            self.execute('fasta', args, stdout=outFa1)
        else:
            self.execute('fasta', args)

    def bam2fa_pipe(self, inBam):
        tool_cmd = [self.install_and_get_path(), 'fasta', '-n', inBam]
        log.debug(' '.join(tool_cmd) + ' |')
        p = subprocess.Popen(tool_cmd, stdout=subprocess.PIPE)
        return p.stdout

    @contextlib.contextmanager
    def bam2fq_tmp(self, inBam):
        with util.file.tempfnames(('.1.fq', '.2.fq')) as (reads1, reads2):
            self.bam2fq(inBam, reads1, reads2)
            yield reads1, reads2

    @contextlib.contextmanager
    def bam2fa_tmp(self, inBam):
        with util.file.tempfnames(('.1.fa', '.2.fa')) as (reads1, reads2):
            self.bam2fa(inBam, reads1, reads2)
            yield reads1, reads2

    def sort(self, inFile, outFile, args=None, threads=None):
        # inFile can be .sam, .bam, .cram
        # outFile can be .sam, .bam, .cram
        args = args or []
        if '-@' not in args:
            args.extend(('-@', str(util.misc.sanitize_thread_count(threads))))
        if '-T' not in args and os.path.isdir(tempfile.tempdir):
            args.extend(('-T', tempfile.tempdir))

        args.extend(('-o', outFile, inFile))
        self.execute('sort', args)

    def merge(self, inFiles, outFile, options=None):
        "Merge a list of inFiles to create outFile."
        options = options or ['-f']

        # We are using -f for now because mkstempfname actually makes an empty
        # file, and merge fails with that as output target without the -f.
        # When mkstempfname is fixed, we should remove the -f.
        self.execute('merge', options + [outFile] + inFiles)

    def index(self, inBam):
        # inBam can be .bam or .cram
        self.execute('index', [inBam])

    def faidx(self, inFasta, overwrite=False):
        ''' Index reference genome for samtools '''
        outfname = inFasta + '.fai'
        if os.path.isfile(outfname):
            if overwrite:
                os.unlink(outfname)
            else:
                return
        #pysam.faidx(inFasta)
        with util.file.fastas_with_sanitized_ids(inFasta, use_tmp=True) as sanitized_fastas:
            sanitized_fasta = sanitized_fastas[0]
            self.execute('faidx', [sanitized_fasta])
            shutil.copyfile(sanitized_fasta + '.fai', outfname)

    def depth(self, inBam, outFile, options=None):
        """ Write a TSV file with coverage depth by position """
        options = options or ["-aa", "-m", "1000000"]

        self.execute('depth', options + [inBam], stdout=outFile)

    def idxstats(self, inBam, statsFile):
        self.execute('idxstats', [inBam], stdout=statsFile)

    def reheader(self, inBam, headerFile, outBam):
        self.execute('reheader', [headerFile, inBam], stdout=outBam)

    def dumpHeader(self, inBam, outHeader):
        opts = []
        if inBam.endswith('.bam'):
            opts = ['-H']
        elif inBam.endswith('.sam'):
            opts = ['-H', '-S']
        #header = pysam.view(*opts)
        self.view(opts, inBam, outHeader)

    def removeDoublyMappedReads(self, inBam, outBam):
        #opts = ['-b', '-f 2']
        opts = ['-b', '-F' '1028', '-f', '2', '-@', '3']
        self.view(opts, inBam, outBam)

    def filter_to_proper_primary_mapped_reads(self, inBam, outBam, require_pairs_to_be_proper=True, reject_singletons=True, reject_duplicates=True):
        '''
            This function writes a bam file filtered to include only reads that are:
              - not flagged as duplicates
              - not secondary or supplementary (split/chimeric reads)
              - For paired-end reads:
              -   marked as proper pair (if require_pairs_to_be_proper=True) OR
                  both not unmapped (if require_pairs_to_be_proper=False) OR
                  not a member of a pair with a singleton (if reject_singletons=True)
              - For single-end reads:
                  mapped

        '''

        with pysam.AlignmentFile(inBam, 'rb', check_sq=False) as inb:
            with pysam.AlignmentFile(outBam, 'wb', header=inb.header) as outf:
                # process the lines individually and write them or not, depending on the flags
                # For explanation of flags, see:
                # https://broadinstitute.github.io/picard/explain-flags.html
                # https://pysam.readthedocs.io/en/latest/api.html
                # https://samtools.github.io/hts-specs/SAMv1.pdf
                # https://github.com/pysam-developers/pysam/blob/31183d7fac52b529b304bdf61ff933818ae4a71f/samtools/stats.c#L72-L81

                for read in inb:
                    # check if a read is paired
                    is_single_end=not read.is_paired

                    # if a PCR/optical duplicate, do not write
                    if reject_duplicates and read.is_duplicate:
                        continue

                    # if a read is a secondary or supplementary mapping (split/chimeric), do not write
                    if read.is_secondary or read.is_supplementary:
                        continue

                    # do not write if
                    # paired-end
                    if (read.is_paired and 
                            # reject anything not marked as proper pair (this bit is not guaranteed)
                            (require_pairs_to_be_proper and not read.is_proper_pair) or
                            # reject pairs where both mates are unmapped
                            (read.mate_is_unmapped and read.is_unmapped) or
                            # reject reads where only one mate is mapped (singletons)
                            (reject_singletons and (read.mate_is_unmapped or read.is_unmapped))):
                        continue 

                    if is_single_end and read.is_unmapped: # or if this is single-end and unmapped, reject
                        continue

                    # otherwise write the read to the output
                    outf.write(read)

    def filterByCigarString(self, inBam, outBam, 
                            regexToMatchForRemoval='^((?:[0-9]+[ID]){1}(?:[0-9]+[MNIDSHPX=])+)|((?:[0-9]+[MNIDSHPX=])+(?:[0-9]+[ID]){1})$', 
                            invertResult=False):
        '''
            This function applies a regex to the cigar string for each read.
            If the regex matches, the read is not written to the output unless
            invertResult is True.
            The default regex matches cigar strings with trailing or leading indels:
              '^((?:[0-9]+[ID]){1}(?:[0-9]+[MNIDSHPX=])+)|((?:[0-9]+[MNIDSHPX=])+(?:[0-9]+[ID]){1})$'

        '''
        regex = re.compile(regexToMatchForRemoval)
        with pysam.AlignmentFile(inBam, 'rb', check_sq=False) as inb:
            with pysam.AlignmentFile(outBam, 'wb', header=inb.header) as outf:
                # process the lines individually and write them or not, depending on 
                # whether they match regexToMatchForRemoval
                for read in inb:
                    if read.cigarstring:
                        # perform a regex match
                        matches = regex.search(read.cigarstring)
                        # if the regex was found (or not, if inverted)
                        if (not invertResult and matches) or (invertResult and not matches):
                            # continue to the next read (don't write out this one)
                            continue
                    else:
                        if invertResult:
                            continue
                    # otherwise write out the line
                    outf.write(read)


    def downsample(self, inBam, outBam, probability):

        if not probability:
            raise Exception("Probability must be defined")
        if float(probability) <= 0 or float(probability) > 1:
            raise Exception("Probability must be in range (0,1]. This value was given: %s" % probability)

        opts = ['-s', str(1) + '.' + str(probability).split('.')[1], '-@', '3']    # -s subsamples: seed.fraction
        self.view(opts, inBam, outBam)

    def downsample_to_approx_count(self, inBam, outBam, read_count):
        total_read_count = self.count(inBam)
        probability = Decimal(int(read_count)) / Decimal(total_read_count)
        probability = 1 if probability > 1 else probability

        if probability < 1:
            self.downsample(inBam, outBam, probability)
        else:
            _log.info("Requested downsample count exceeds number of reads. Including all reads in output.")
            shutil.copyfile(inBam, outBam)

    def getHeader(self, inBam):
        ''' fetch BAM header as a list of tuples (already split on tabs) '''
        tmpf = util.file.mkstempfname('.txt')
        self.dumpHeader(inBam, tmpf)
        with open(tmpf, 'rb') as inf:
            header = list(line.decode("latin-1").rstrip('\n').split('\t') for line in inf)
        os.unlink(tmpf)
        return header

    def getReadGroups(self, inBam):
        ''' fetch all read groups from the BAM header as an OrderedDict of
            RG ID -> RG dict.  The RG dict is a mapping of read group keyword
            (like ID, DT, PU, LB, SM, CN, PL, etc) to value.  ID is included
            and not stripped out. ID is required for all read groups.
            Resulting keys are in same order as @RG lines in bam file.
        '''
        rgs = [
            dict(x.split(':', 1) for x in row[1:]) for row in self.getHeader(inBam) if len(row) > 0 and row[0] == '@RG'
        ]
        return OrderedDict((rg['ID'], rg) for rg in rgs)

    def count(self, inBam, opts=None, regions=None):
        opts = opts or []
        regions = regions or []

        cmd = [self.install_and_get_path(), 'view', '-c'] + opts + [inBam] + regions
        return int(subprocess.check_output(cmd).strip())
        #if inBam.endswith('.sam') and '-S' not in opts:
        #    opts = ['-S'] + opts
        #cmd = ['-c'] + opts + [inBam] + regions
        #return int(pysam.view(*cmd)[0].strip())

    def mpileup(self, inBam, outPileup, opts=None):
        opts = opts or []

        self.execute('mpileup', opts + [inBam], stdout=outPileup, stderr='/dev/null')    # Suppress info messages

    def isEmpty(self, inBam):
        if not os.path.isfile(inBam):
            return True
        else:
            tmpf = util.file.mkstempfname('.txt')
            self.dumpHeader(inBam, tmpf)
            header_size = os.path.getsize(tmpf)
            if (os.path.getsize(inBam) > (100 + 5 * header_size)):
                # large BAM file: assume it is not empty
                # a BAM file definitely has reads in it if its filesize is larger
                # than just the header itself
                return False
            else:
                # small BAM file: just count and see if it's non-zero
                return (0 == self.count(inBam))
