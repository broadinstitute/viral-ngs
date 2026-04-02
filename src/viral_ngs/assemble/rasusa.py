'''
    Rasusa — randomly subsample sequencing reads or alignments to a
    specified coverage depth.

    https://github.com/mbhall88/rasusa
'''

import logging
import shutil
import subprocess

from ..core import Tool, PrexistingUnixCommand

_log = logging.getLogger(__name__)

TOOL_NAME = 'rasusa'


class RasusaTool(Tool):

    def __init__(self):
        install_methods = [PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        Tool.__init__(self, install_methods=install_methods)

    def downsample_bam(self, inBam, outBam, coverage, seed=None):
        """Downsample aligned reads in a BAM to a target coverage depth
        using ``rasusa aln``.

        The input BAM must be coordinate-sorted.  The output BAM is
        **not** guaranteed to be sorted — the caller is responsible for
        sorting and indexing if downstream tools require it.

        Args:
            inBam:    Input BAM file (coordinate-sorted).
            outBam:   Output BAM file path.
            coverage: Target coverage depth (integer).
            seed:     Random seed for reproducibility (optional).
        """
        tool_cmd = [
            self.install_and_get_path(), 'aln',
            '--coverage', str(coverage),
            '-o', outBam,
            inBam,
        ]
        if seed is not None:
            tool_cmd.extend(['--seed', str(seed)])
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
