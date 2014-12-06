import logging, tools, util.file
import os, os.path, subprocess

'''
TO DO: much of this stuff can be all eliminated by using pysam instead:
http://pysam.readthedocs.org/en/latest/usage.html#using-samtools-commands-within-python
'''


tool_version = '0.1.19'
url = 'http://sourceforge.net/projects/samtools/files/samtools/' \
    + '{ver}/samtools-{ver}.tar.bz2'.format(ver=tool_version)

log = logging.getLogger(__name__)

class SamtoolsTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = [
                tools.DownloadPackage(url, 'samtools-{}/samtools'.format(tool_version),
                    post_download_command='cd samtools-{}; make'.format(tool_version))]
            #path = '/idi/sabeti-data/software/samtools/samtools-0.1.19/samtools',
            #install_methods.append(tools.PrexistingUnixCommand(path))
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
    
    def faidx(self, inFasta, overwrite=False) :
        ''' Index reference genome for samtools '''
        outfname = inFasta + '.fai'
        if os.path.isfile(outfname):
            if overwrite:
                os.unlink(outfname)
            else:
                return
        self.execute('faidx', [inFasta])
    
    def reheader(self, inBam, headerFile, outBam) :
        self.execute('reheader', [headerFile, inBam], stdout=outBam)
    
    def dumpHeader(self, inBam, outHeader) :
        self.execute('view', ['-H', inBam], stdout=outHeader)
    
    def count(self, inBam, opts=[]) :
        cmd = [self.install_and_get_path(), 'view', '-c'] + opts + [inBam]
        return int(subprocess.check_output(cmd).strip())
