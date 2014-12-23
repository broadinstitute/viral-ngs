'''
    The Trinity RNA-SEQ assembler
    
    This uses an older version of Trinity that uses an older
    assembly algorithm that works better with highly diverse
    viral genomes.
'''

import logging, tools, util.file
import os, os.path, subprocess

tool_version = "2011-11-26"
url = "http://sourceforge.net/projects/trinityrnaseq/files/" \
    + "trinityrnaseq_r{}.tgz/download".format(tool_version)

log = logging.getLogger(__name__)


'''
Hm, compiling this is a little tricky...

class TrinityTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = [tools.DownloadPackage(url,
                'trinityrnaseq_r{}/Trinity.pl'.format(tool_version),
                post_download_command='cd trinityrnaseq_r{}; make'.format(tool_version))]
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def execute(self, inFastq1, inFastq2, outFasta):
        raise NotImplementedException()
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

        """
        shutil.rmtree(params.tmpd_trinity, ignore_errors=True)
        shell("reuse -q Java-1.6 && perl /idi/sabeti-scratch/kandersen/bin/trinity_old/Trinity.pl --CPU 1 --min_contig_length 300 --seqType fq --left {params.tmpf_subsamp[0]} --right {params.tmpf_subsamp[1]} --output {params.tmpd_trinity}")
        shutil.copyfile(params.tmpd_trinity+"/Trinity.fasta", output[0])
        map(os.unlink, params.tmpf_subsamp)
        shutil.rmtree(params.tmpd_trinity, ignore_errors=True)
        """
'''
    

