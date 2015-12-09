'''
    The Samtools package.

    TO DO: much of this stuff can be all eliminated by using pysam instead, as
    pysam (in its current versions) is meant to be the complete python
    implementation of htslib/samtools.

    http://pysam.readthedocs.org/en/latest/usage.html#using-samtools-commands-within-python

    Current bug with pysam 0.8.1: nosetests does not work unless you use --nocapture.
    python -m unittest works. Something about redirecting stdout.
    Actually, Travis CI still has issues with pysam and stdout even with --nocapture.
'''

import logging
import tools
import util.file
import os
import os.path
import subprocess
from collections import OrderedDict
#import pysam

TOOL_NAME = 'samtools'
TOOL_VERSION = '0.1.19'
CONDA_TOOL_VERSION = '1.2'
TOOL_URL = 'http://sourceforge.net/projects/samtools/files/samtools/' \
    + '{ver}/samtools-{ver}.tar.bz2'.format(ver=TOOL_VERSION)

log = logging.getLogger(__name__)


class SamtoolsTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append( tools.CondaPackage(TOOL_NAME, version=CONDA_TOOL_VERSION) )
            install_methods.append(tools.DownloadPackage(TOOL_URL,
                                      'samtools-{}/samtools'.format(TOOL_VERSION),
                                      post_download_command='cd samtools-{}; make -s'.format(TOOL_VERSION)))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, args, stdin=None, stdout=None, stderr=None): # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdin:
            stdin = open(stdin, 'r')
        if stdout:
            stdout = open(stdout, 'w')
        if stderr:
            stderr = open(stderr, 'w')
        subprocess.check_call(tool_cmd, stdin=stdin, stdout=stdout, stderr=stderr)
        if stdin:
            stdin.close()
        if stdout:
            stdout.close()
        if stderr:
            stderr.close()

    def view(self, args, inFile, outFile, regions=None):
        regions = regions or []

        self.execute('view', args + ['-o', outFile, inFile] + regions)

    def merge(self, inFiles, outFile, options=None):
        "Merge a list of inFiles to create outFile."
        options = options or ['-f']
        
        # We are using -f for now because mkstempfname actually makes an empty
        # file, and merge fails with that as output target without the -f.
        # When mkstempfname is fixed, we should remove the -f.
        self.execute('merge', options + [outFile] + inFiles)

    def index(self, inBam):
        self.execute('index', [inBam])

    def faidx(self, inFasta, overwrite=False):
        ''' Index reference genome for samtools '''
        outfname = inFasta + '.fai'
        if os.path.isfile(outfname):
            if overwrite:
                os.unlink(outfname)
            else:
                return
        # pysam.faidx(inFasta)
        self.execute('faidx', [inFasta])

    def reheader(self, inBam, headerFile, outBam):
        self.execute('reheader', [headerFile, inBam], stdout=outBam)

    def dumpHeader(self, inBam, outHeader):
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

    def getHeader(self, inBam):
        ''' fetch BAM header as a list of tuples (already split on tabs) '''
        tmpf = util.file.mkstempfname('.txt')
        self.dumpHeader(inBam, tmpf)
        with open(tmpf, 'rt') as inf:
            header = list(line.rstrip('\n').split('\t') for line in inf)
        os.unlink(tmpf)
        return header

    def getReadGroups(self, inBam):
        ''' fetch all read groups from the BAM header as an OrderedDict of
            RG ID -> RG dict.  The RG dict is a mapping of read group keyword
            (like ID, DT, PU, LB, SM, CN, PL, etc) to value.  ID is included
            and not stripped out. ID is required for all read groups.
            Resulting keys are in same order as @RG lines in bam file.
        '''
        rgs = [dict(x.split(':', 1) for x in row[1:]) for row in self.getHeader(inBam)
               if len(row) > 0 and row[0] == '@RG']
        return OrderedDict((rg['ID'], rg) for rg in rgs)

    def count(self, inBam, opts=None, regions=None):
        opts = opts or []
        regions = regions or []

        cmd = [self.install_and_get_path(), 'view', '-c'] + opts + [inBam] + regions
        # return int(pysam.view(*cmd)[0].strip())
        return int(subprocess.check_output(cmd).strip())

    def mpileup(self, inBam, outPileup, opts=None):
        opts = opts or []

        self.execute('mpileup', opts + [inBam], stdout=outPileup, stderr='/dev/null')  # Suppress info messages
