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
        toolCmd = ['/'.join(self.install_and_get_path(), 'mummer'), refFasta] + qryFastas
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def nucmer(self, refFasta, qryFasta, outDelta):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = ['/'.join(self.install_and_get_path(), 'nucmer'),
            '--prefix={}'.format(outDelta), refFasta, qryFasta]
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def promer(self, refFasta, qryFasta, outDelta):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = ['/'.join(self.install_and_get_path(), 'nucmer'),
            '--prefix={}'.format(outDelta), refFasta, qryFasta]
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)
        

