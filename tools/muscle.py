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

tool_version = '3.8.31'
url = 'http://www.drive5.com/muscle/downloads{ver}/muscle{ver}_{os}.tar.gz'

log = logging.getLogger(__name__)


class MuscleTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            muscle_os = get_muscle_os()
            if muscle_os != 'src':
                install_methods.append(
                    tools.DownloadPackage(url.format(ver=tool_version,
                                                     os=muscle_os),
                                          'muscle{}_{}'.format(tool_version, muscle_os),
                                          verifycmd='{}/muscle{}_{} -version > /dev/null 2>&1'.format(
                                              util.file.get_build_path(), tool_version, muscle_os)))
            install_methods.append(
                tools.DownloadPackage(url.format(ver=tool_version,
                                                 os=muscle_os),
                                      'muscle{}/src/muscle'.format(tool_version),
                                      post_download_command='cd muscle{}/src; make -s'.format(tool_version),
                                      verifycmd='{}/muscle{}/src/muscle -version > /dev/null 2>&1'.format(
                                          util.file.get_build_path(), tool_version)))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return tool_version

    def execute(self, inFasta, outFasta,
                maxiters=None, maxhours=None, format='fasta', diags=None, quiet=True, logFile=None):
        toolCmd = [self.install_and_get_path(), '-in', inFasta, '-out', outFasta]

        if format in ('html', 'msf', 'clw', 'clwstrict'):
            toolCmd.append('-' + format)
        else:
            if format != 'fasta':
                raise Exception()
        if quiet:
            toolCmd.append('-quiet')
        if diags:
            toolCmd.append('-diags')
        if maxiters:
            toolCmd.append('-maxiters {}'.format(maxiters))
        if maxhours:
            toolCmd.append('-maxhours {}'.format(maxhours))
        if logFile:
            toolCmd.append('-log {}'.format(logFile))

        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)


def get_muscle_os():
    uname = os.uname()
    if uname[4].startswith('x86') and uname[0] in ('Darwin', 'Linux'):
        return 'i86' + uname[0].lower() + uname[4][-2:]
    return 'src'
