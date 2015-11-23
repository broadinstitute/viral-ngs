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

#tool_version = '3.8.31'
#url = 'http://www.drive5.com/muscle/downloads{ver}/muscle{ver}_{os}.tar.gz'
tool_version = '3.8.1551'
url = 'http://www.drive5.com/muscle/muscle_{os}_{ver}.tar.gz'

log = logging.getLogger(__name__)


class MuscleTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            '''
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
            '''
            install_methods.append(DownloadAndBuildMuscle(tool_version))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return tool_version

    def execute(self, inFasta, outFasta,
                maxiters=None, maxhours=None, fmt='fasta', diags=None, quiet=True, logFile=None):
        toolCmd = [self.install_and_get_path(), '-in', inFasta, '-out', outFasta]

        if fmt in ('html', 'msf', 'clw', 'clwstrict'):
            toolCmd.append('-' + fmt)
        else:
            if fmt != 'fasta':
                raise Exception()
        if quiet:
            toolCmd.append('-quiet')
        if diags:
            toolCmd.append('-diags')
        if maxiters:
            toolCmd.extend(('-maxiters', str(maxiters)))
        if maxhours:
            toolCmd.extend(('-maxhours', str(maxhours)))
        if logFile:
            toolCmd.extend(('-log', logFile))

        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

'''
def get_muscle_os():
    uname = os.uname()
    if uname[4].startswith('x86') and uname[0] in ('Darwin', 'Linux'):
        return 'i86' + uname[0].lower() + uname[4][-2:]
    return 'src'
'''

class DownloadAndBuildMuscle(tools.DownloadPackage):

    def __init__(self, ver):
        tools.DownloadPackage.__init__(self,
            url.format(ver=ver, os='src'),
            'muscle',
            destination_dir=os.path.join(
                util.file.get_build_path(),
                'muscle-{}'.format(ver)),
            verifycmd='{}/muscle-{}/muscle -version > /dev/null 2>&1'.format(
                util.file.get_build_path(), tool_version),
            post_download_command='make -s')

    def post_download(self):
        if os.uname()[0] == 'Darwin':
            # Makefile needs to be modified for MacOSX
            badFilePath = os.path.join(self.destination_dir, 'Makefile')
            os.rename(badFilePath, badFilePath + '.orig')
            with open(badFilePath, 'wt') as outf:
                with open(badFilePath + '.orig', 'rt') as inf:
                    for line in inf:
                        if line.strip() == 'LDLIBS = -lm -static':
                            line = '# ' + line
                        elif line.strip() == '# LDLIBS = -lm':
                            line = line[2:]
                        outf.write(line)
        tools.DownloadPackage.post_download(self)

