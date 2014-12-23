'''
    The Samtools package.
    
    TO DO: much of this stuff can be all eliminated by using pysam instead, as
    pysam (in its current versions) is meant to be the complete python
    implementation of htslib/samtools.
    
    http://pysam.readthedocs.org/en/latest/usage.html#using-samtools-commands-within-python
    
    Current bug with pysam 0.8.1: nosetests does not work unless you use --nocapture.
    python -m unittest works. Something about redirecting stdout.
'''

import logging, tools, util.file
import os, os.path, subprocess
import pysam

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
    
    def faidx(self, inFasta, overwrite=False):
        ''' Index reference genome for samtools '''
        outfname = inFasta + '.fai'
        if os.path.isfile(outfname):
            if overwrite:
                os.unlink(outfname)
            else:
                return
        pysam.faidx(inFasta)
    
    def reheader(self, inBam, headerFile, outBam):
        self.execute('reheader', [headerFile, inBam], stdout=outBam)
    
    def dumpHeader(self, inBam, outHeader):
        with open(outHeader, 'wt') as outf:
            for row in self.getHeader(inBam):
                outf.write('\t'.join(row)+'\n')
    
    def getHeader(self, inBam):
        if inBam.endswith('.bam'):
            header = pysam.view('-H', inBam)
        elif inBam.endswith('.sam'):
            header = pysam.view('-H', '-S', inBam)
        return list(line.rstrip('\n').split('\t') for line in header)
    
    def count(self, inBam, opts=[], regions=[]):
        cmd = opts + ['-c', inBam] + regions
        if inBam.endswith('.sam') and '-S' not in opts:
            cmd = ['-S'] + cmd
        return int(pysam.view(*cmd)[0].strip())
