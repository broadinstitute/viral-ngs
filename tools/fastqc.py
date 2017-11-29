'''
    FastQC
'''

import logging
import os
import os.path
import shutil
import subprocess
import sys

import tools
import tools.samtools
import util.file
import util.misc

TOOL_NAME = 'fastqc'
TOOL_VERSION = '0.11.5'

log = logging.getLogger(__name__)

class FastQC(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, inBam, outHtml):    # pylint: disable=W0221
        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # fastqc can't deal with empty input
            with open(outHtml, 'wt') as outf:
                outf.write("<html><body>Input BAM has zero reads.</body></html>\n")

        else:
            # run fastqc
            with util.file.tmp_dir() as out_dir:
                tool_cmd = [self.install_and_get_path(),
                    '-t', str(util.misc.sanitize_thread_count()),
                    '-o', out_dir,
                    inBam]
                log.debug(' '.join(tool_cmd))
                subprocess.check_call(tool_cmd, stdout=sys.stderr)
                expected_out = os.path.join(out_dir, os.path.basename(inBam)[:-4]) + "_fastqc.html"
                shutil.copyfile(expected_out, outHtml)
                log.debug("complete")
