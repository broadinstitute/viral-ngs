'''
    The MUSCLE aligner
    http://www.drive5.com/muscle
'''


import tools
import util.file
import util.misc

import os
import os.path
import subprocess
import logging

TOOL_NAME = "muscle"
TOOL_VERSION = '3.8.1551'
CONDA_TOOL_VERSION = '3.8.1551'

_log = logging.getLogger(__name__)


class MuscleTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=CONDA_TOOL_VERSION))
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

