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
from collections import OrderedDict
from decimal import *

#import pysam

import tools
import util.file
import util.misc

TOOL_NAME = 'samtools'
TOOL_VERSION = '1.3.1'

log = logging.getLogger(__name__)


class SamtoolsTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, args, stdout=None, stderr=None):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        if stderr:
            stderr = open(stderr, 'w')
        subprocess.check_call(tool_cmd, stdout=stdout, stderr=stderr)
        if stdout:
            stdout.close()
        if stderr:
            stderr.close()

    def view(self, args, inFile, outFile, regions=None):
        regions = regions or []

        self.execute('view', args + ['-o', outFile, inFile] + regions)
        #opts = args + ['-o', outFile, inFile] + regions
        #pysam.view(*opts)

    def bam2fq(self, inBam, outFq1, outFq2=None):
        if outFq2 is None:
            self.execute('bam2fq', ['-n', '-0', outFq1, inBam])
        else:
            self.execute('bam2fq', ['-1', outFq1, '-2', outFq2, inBam])

    def sort(self, inFile, outFile, args=None, threads=None):
        # inFile can be .sam, .bam, .cram
        # outFile can be .sam, .bam, .cram
        args = args or []
        threads = threads or util.misc.available_cpu_count()
        if '-@' not in args:
            args.extend(('-@', str(threads)))
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
        self.execute('faidx', [inFasta])

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
        opts = ['-b', '-F' '1028', '-f', '2']
        self.view(opts, inBam, outBam)

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
        in_args = [self.install_and_get_path(), 'view', '-h', inBam]

        regex = re.compile(regexToMatchForRemoval)

        # an input samtools process to read a bam file
        # bufsize 0=unbuffered, 1=linebuffered, -1=system buffered (usually fully)
        in_process = subprocess.Popen(in_args, bufsize=1, stdout=subprocess.PIPE,
                                    universal_newlines=True)

        # an output samtools process to write filtered output via stdin
        out_args = [self.install_and_get_path(), 'view', '-h', '-o', outBam]
        out_process = subprocess.Popen(out_args, stdin=subprocess.PIPE)

        # process the lines individually and write them or not, depending on 
        # whether they match regexToMatchForRemoval
        # Use while and readline() instead of "for line in in_process.stdout"
        # to avoid Python readahead bug: https://bugs.python.org/issue3907
        while True:
            line = in_process.stdout.readline()
            if not line or line=='':
                break
            # always write header lines, so skip the cigar string check
            if not line.startswith('@'):
                # cigar strings are in column 6
                cigar_string = line.split('\t')[5]

                # perform a regex match
                matches = regex.search(cigar_string)
                # if the regex was found (or not, if inverted)
                if (not invertResult and matches) or (invertResult and not matches):
                    # continue to the next read (don't write out this one)
                    continue

            # otherwise write out the line to samtools' stdin
            out_process.stdin.write(line.encode("utf-8"))

        # close stdin stream of the output process so it can terminate
        out_process.stdin.close()


    def downsample(self, inBam, outBam, probability):

        if not probability:
            raise Exception("Probability must be defined")
        if float(probability) <= 0 or float(probability) > 1:
            raise Exception("Probability must be in range (0,1]. This value was given: %s" % probability)

        opts = ['-s', str(1) + '.' + str(probability).split('.')[1]]    # -s subsamples: seed.fraction
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
