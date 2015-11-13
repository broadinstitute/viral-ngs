'''
    The MUMMER aligner
    http://mummer.sourceforge.net/
'''

import logging
import tools
import util.file
import os
import os.path
import subprocess

tool_version = '3.23'
url = 'http://iweb.dl.sourceforge.net/project/mummer/mummer/{ver}/MUMmer{ver}.tar.gz'

log = logging.getLogger(__name__)


class MummerTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [
                tools.DownloadPackage(url.format(ver=tool_version),
                                      'MUMmer{}'.format(tool_version),
                                      post_download_command='cd MUMmer{}; make -s'.format(tool_version),
                                      verifycmd='{}/MUMmer{}/mummer -h > /dev/null 2>&1'.format(
                                          util.file.get_build_path(), tool_version))
                ]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return tool_version

    def execute(self, refFasta, qryFastas):
        toolCmd = [os.path.join(self.install_and_get_path(), 'mummer'),
            refFasta] + qryFastas
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def nucmer(self, refFasta, qryFasta, outDelta):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = [os.path.join(self.install_and_get_path(), 'nucmer'),
            '--prefix={}'.format(outDelta), refFasta, qryFasta]
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def promer(self, refFasta, qryFasta, outDelta):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = [os.path.join(self.install_and_get_path(), 'promer'),
            '--prefix={}'.format(outDelta), refFasta, qryFasta]
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)
    
    def delta_filter(self, inDelta, outDelta):
        toolCmd = [os.path.join(self.install_and_get_path(), 'delta-filter'),
            '-q', inDelta]
        log.debug(' '.join(toolCmd))
        with open(outDelta, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)

    def show_tiling(self, inDelta, outTiling, outFasta,
            circular=False, min_pct_id=None, min_contig_len=None):
        opts = []
        if circular:
            opts.append('-c')
        if min_pct_id is not None:
            opts.append('-i')
            opts.append(str(min_pct_id))
        if min_contig_len is not None:
            opts.append('-l')
            opts.append(str(min_contig_len))
        toolCmd = [os.path.join(self.install_and_get_path(), 'show-tiling'),
            '-p', outFasta] + opts + [inDelta]
        log.debug(' '.join(toolCmd))
        with open(outTiling, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)
    
    def scaffold_contigs(self, refFasta, contigsFasta, outFasta,
            aligner='nucmer', circular=False,
            min_pct_id=0.6, min_contig_len=200):
        
        if aligner=='nucmer':
            aligner = self.nucmer
        elif aligner=='promer':
            aligner = self.promer
        else:
            raise NameError()
        delta_1 = util.file.mkstempfname('.delta')
        delta_2 = util.file.mkstempfname('.delta')
        tiling = util.file.mkstempfname('.tiling')
        aligner(refFasta, contigsFasta, delta_1)
        self.delta_filter(delta_1, delta_2)
        self.show_tiling(delta_2, tiling, outFasta, circular=circular,
            min_pct_id=min_pct_id, min_contig_len=min_contig_len)
        os.unlink(delta_1)
        os.unlink(delta_2)
        os.unlink(tiling)
        

