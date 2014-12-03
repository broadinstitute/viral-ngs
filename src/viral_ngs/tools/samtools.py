import logging, tools, util.file
import os, os.path, subprocess

url = 'http://sourceforge.net/projects/samtools/files/samtools/0.1.19/' \
        + 'samtools-0.1.19.tar.bz2'

log = logging.getLogger(__name__)

class SamtoolsTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = []
            install_methods.append(
                tools.DownloadPackage(url, 'samtools-0.1.19/samtools',
                    post_download_command='cd samtools-0.1.19; make'))
            #path = '/idi/sabeti-data/software/samtools/samtools-0.1.19/samtools',
            #install_methods.append(tools.PrexistingUnixCommand(path))
        tools.Tool.__init__(self, install_methods = install_methods)
    
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
        outfname = fasta + '.fai'
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
        tmp = util.file.mkstempfname('.count')
        self.execute('view', ['-c'] + opts, stdout=tmp)
        with open(tmp, 'rt') as inf:
            return int(inf.readline().strip())
