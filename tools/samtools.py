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

import logging, tools, util.file
import os, os.path, subprocess
#import pysam

tool_version = '0.1.19'
url = 'http://sourceforge.net/projects/samtools/files/samtools/' \
    + '{ver}/samtools-{ver}.tar.bz2'.format(ver=tool_version)

log = logging.getLogger(__name__)

class SamtoolsTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = [
                tools.DownloadPackage(url, 'samtools-{}/samtools'.format(tool_version),
                    post_download_command='cd samtools-{}; make -s'.format(tool_version))]
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def execute(self, command, args, stdin=None, stdout=None):
        toolCmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(toolCmd))
        if stdin:
            stdin = open(stdin, 'r')
        if stdout:
            stdout = open(stdout, 'w')
        subprocess.check_call(toolCmd, stdin=stdin, stdout=stdout)
        if stdin:
            stdin.close()
        if stdout:
            stdout.close()

    def view(self, args, inFile, outFile, regions=[]):
        self.execute('view', args + ['-o', outFile, inFile] + regions)
    
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
    
    def reheader(self, inBam, headerFile, outBam):
        self.execute('reheader', [headerFile, inBam], stdout=outBam)
    
    def dumpHeader(self, inBam, outHeader) :
        if inBam.endswith('.bam'):
            opts = ['-H']
        elif inBam.endswith('.sam'):
            opts = ['-H', '-S']
        #header = pysam.view(*opts)
        self.view(opts, inBam, outHeader)
    
    def getHeader(self, inBam):
        tmpf = util.file.mkstempfname('.txt')
        self.dumpHeader(inBam, tmpf)
        with open(tmpf, 'rt') as inf:
            header = list(line.rstrip('\n').split('\t') for line in inf)
        os.unlink(tmpf)
        return header
    
    def count(self, inBam, opts=[], regions=[]):
        cmd = [self.install_and_get_path(), 'view', '-c'] + opts + [inBam] + regions
        #return int(pysam.view(*cmd)[0].strip())
        return int(subprocess.check_output(cmd).strip())
