'''
    FreeBayes variant caller — a Bayesian genetic variant detector.

    Replaces GATK3 UnifiedGenotyper for consensus refinement.
    Uses --report-monomorphic to emit all sites (equivalent to EMIT_ALL_SITES).
'''

import logging
import shutil
import subprocess

from ..core import Tool, PrexistingUnixCommand

_log = logging.getLogger(__name__)

TOOL_NAME = 'freebayes'


class FreeBayesTool(Tool):

    def __init__(self):
        install_methods = [PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        Tool.__init__(self, install_methods=install_methods)

    def call(self, inBam, refFasta, outVcf, options=None):
        """Call variants with FreeBayes, emitting all sites."""
        options = options or ["--min-base-quality", "15", "--pooled-continuous"]
        opts = [
            '-b', inBam,
            '-f', refFasta,
            '-v', outVcf,
            '--use-best-n-alleles', '0',
            '--report-monomorphic',
        ]
        tool_cmd = [self.install_and_get_path()] + opts + options
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
