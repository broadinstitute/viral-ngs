'''
    The MUSCLE aligner
    http://www.drive5.com/muscle
'''

import logging
import tools
import util.file
import os
import os.path
import subprocess

TOOL_VERSION = '3.8.31'
TOOL_URL = 'http://www.drive5.com/muscle/downloads{ver}/muscle{ver}_{os}.tar.gz'

LOG = logging.getLogger(__name__)


class MuscleTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            muscle_os = get_muscle_os()
            if muscle_os != 'src':
                install_methods.append(
                    tools.DownloadPackage(TOOL_URL.format(ver=TOOL_VERSION,
                                                     os=muscle_os),
                                          'muscle{}_{}'.format(TOOL_VERSION, muscle_os),
                                          verifycmd='{}/muscle{}_{} -version > /dev/null 2>&1'.format(
                                              util.file.get_build_path(), TOOL_VERSION, muscle_os)))
            install_methods.append(
                tools.DownloadPackage(TOOL_URL.format(ver=TOOL_VERSION,
                                                 os=muscle_os),
                                      'muscle{}/src/muscle'.format(TOOL_VERSION),
                                      post_download_command='cd muscle{}/src; make -s'.format(TOOL_VERSION),
                                      verifycmd='{}/muscle{}/src/muscle -version > /dev/null 2>&1'.format(
                                          util.file.get_build_path(), TOOL_VERSION)))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, inFasta, outFasta,
                maxiters=None, maxhours=None, fmt='fasta', diags=None, quiet=True, logFile=None): # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), '-in', inFasta, '-out', outFasta]

        if fmt in ('html', 'msf', 'clw', 'clwstrict'):
            tool_cmd.append('-' + fmt)
        else:
            if fmt != 'fasta':
                raise Exception()
        if quiet:
            tool_cmd.append('-quiet')
        if diags:
            tool_cmd.append('-diags')
        if maxiters:
            tool_cmd.append('-maxiters {}'.format(maxiters))
        if maxhours:
            tool_cmd.append('-maxhours {}'.format(maxhours))
        if logFile:
            tool_cmd.append('-log {}'.format(logFile))

        LOG.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)


def get_muscle_os():
    uname = os.uname()
    if uname[4].startswith('x86') and uname[0] in ('Darwin', 'Linux'):
        return 'i86' + uname[0].lower() + uname[4][-2:]
    return 'src'
