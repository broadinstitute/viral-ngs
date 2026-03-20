'''
    FreeBayes variant caller — a Bayesian genetic variant detector.

    Replaces GATK3 UnifiedGenotyper for consensus refinement.
    Uses --report-monomorphic to emit all sites (equivalent to EMIT_ALL_SITES).
'''

import logging
import os
import shutil
import subprocess

import pysam

from ..core import Tool, PrexistingUnixCommand
from ..core.file import mkstempfname

_log = logging.getLogger(__name__)

TOOL_NAME = 'freebayes'


class FreeBayesTool(Tool):

    def __init__(self):
        install_methods = [PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        Tool.__init__(self, install_methods=install_methods)

    def call(self, inBam, refFasta, outVcf, options=None):
        """Call variants with FreeBayes, emitting all sites.

        If outVcf ends with .vcf.gz, the output is bgzipped and tabix-indexed
        (FreeBayes always writes plain VCF, so we compress and index afterward).
        """
        options = options or ["--min-base-quality", "15", "--pooled-continuous"]

        if outVcf.endswith('.vcf.gz'):
            plainVcf = mkstempfname('.vcf')
        else:
            plainVcf = outVcf

        opts = [
            '-b', inBam,
            '-f', refFasta,
            '-v', plainVcf,
            '--use-best-n-alleles', '0',
            '--report-monomorphic',
        ]
        tool_cmd = [self.install_and_get_path()] + opts + options
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

        if outVcf.endswith('.vcf.gz'):
            pysam.tabix_compress(plainVcf, outVcf, force=True)
            pysam.tabix_index(outVcf, preset='vcf', force=True)
            os.unlink(plainVcf)
