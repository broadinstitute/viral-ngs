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

log = logging.getLogger(__name__)

class FastQC(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
           install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(FastQC, self).__init__(install_methods=install_methods)

    def version(self):
        return subprocess.check_output([self.install_and_get_path(), '-v']).decode('UTF-8').strip().split()[1]

    def execute(self, inBam, out_html, out_zip=None, threads=None):    # pylint: disable=W0221
        threads =  util.misc.sanitize_thread_count(threads)

        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # fastqc can't deal with empty input
            with open(out_html, 'wt') as outf:
                outf.write("<html><body>Input BAM has zero reads.</body></html>\n")
            if out_zip:
                util.file.touch(out_zip)

        else:
            # run fastqc
            with util.file.tmp_dir() as out_dir:

                # fastqc sets java "-Xmx" (max java heap size) to 250mb/thread
                tool_cmd = [self.install_and_get_path(),
                    '-t', str(threads),
                    '-o', out_dir,
                    inBam]
                log.debug(' '.join(tool_cmd))
                subprocess.check_call(tool_cmd, stdout=sys.stderr)
                expected_out = os.path.join(out_dir, os.path.basename(inBam)[:-4]) + "_fastqc.html"
                shutil.copyfile(expected_out, out_html)
                if out_zip:
                    expected_out_zip = os.path.join(out_dir, os.path.basename(inBam)[:-4]) + "_fastqc.zip"
                    shutil.copyfile(expected_out_zip, out_zip)
                log.debug("complete")
