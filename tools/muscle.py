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

TOOL_NAME = "muscle"
TOOL_VERSION = '3.8.1551'
CONDA_TOOL_VERSION = '3.8.1551'
TOOL_URL = 'http://www.drive5.com/muscle/downloads{ver}/muscle{ver}_{os}.tar.gz'

UNRELEASED_URL = 'http://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz'

_log = logging.getLogger(__name__)


class MuscleTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=CONDA_TOOL_VERSION))
            install_methods.append(
                tools.DownloadPackage(
                    UNRELEASED_URL,
                    'muscle{}/src/muscle'.format(TOOL_VERSION),
                    post_download_command='cd muscle{}/src; make -s'.format(TOOL_VERSION),
                    verifycmd='{}/muscle{}/src/muscle -version > /dev/null 2>&1'.format(
                        util.file.get_build_path(), TOOL_VERSION
                    )
                )
            )
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    # pylint: disable=W0221
    def execute(
        self,
        inFasta,
        outFasta,
        maxiters=None,
        maxhours=None,
        fmt='fasta',
        diags=None,
        quiet=True,
        logFile=None
    ):
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
            tool_cmd.extend(('-maxiters', str(maxiters)))
        if maxhours:
            tool_cmd.extend(('-maxhours', str(maxhours)))
        if logFile:
            tool_cmd.extend(('-log', logFile))

        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
    # pylint: enable=W0221


class DownloadAndBuildMuscle(tools.DownloadPackage):

    def __init__(self, ver):
        tools.DownloadPackage.__init__(self,
            TOOL_URL.format(ver=ver, os='src'),
            'muscle',
            destination_dir=os.path.join(
                util.file.get_build_path(),
                'muscle-{}'.format(ver)),
            verifycmd='{}/muscle-{}/muscle -version > /dev/null 2>&1'.format(
                util.file.get_build_path(), TOOL_VERSION),
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

