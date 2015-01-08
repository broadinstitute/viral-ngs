'''
    The Trinity RNA-SEQ assembler
    
    This uses an older version of Trinity that uses an older
    assembly algorithm that works better with highly diverse
    viral genomes.
'''

import logging, tools, util.file
import os, os.path, subprocess, tempfile

tool_version = "2011-11-26"
url = "http://sourceforge.net/projects/trinityrnaseq/files/" \
    + "trinityrnaseq_r{}.tgz/download".format(tool_version)

log = logging.getLogger(__name__)


class TrinityTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        #Hm, compiling this is a little tricky... this installer needs more work
        if install_methods == None :
            install_methods = [tools.DownloadPackage(url,
                'trinityrnaseq_r{}/Trinity.pl'.format(tool_version),
                post_download_command='cd trinityrnaseq_r{}; make'.format(tool_version))]
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def execute(self, inFastq1, inFastq2, outFasta, min_contig_length=300):
        outdir = tempfile.mkdtemp(prefix='trinity-')
        cmd = [self.install_and_get_path(),
            '--CPU', '1',
            '--min_contig_length', str(min_contig_length),
            '--seqType', 'fq',
            '--left', inFastq1,
            '--right', inFastq2,
            '--output', outdir]
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)
        shutil.copyfile(os.path.join(outdir, 'Trinity.fasta'), outFasta)
        shutil.rmtree(outdir, ignore_errors=True)
