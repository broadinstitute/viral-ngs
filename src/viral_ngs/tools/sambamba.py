'''
    Sambamba - high-performance BAM processing tool.

    Sambamba is a tool for processing BAM files that provides multi-threaded
    implementations of common operations like sorting, indexing, merging,
    and duplicate marking. It is significantly faster than single-threaded
    alternatives like Picard for many operations.

    https://github.com/biod/sambamba
'''

import logging
import os
import re
import shutil
import subprocess

import tools
import util.misc

TOOL_NAME = 'sambamba'

log = logging.getLogger(__name__)


class SambambaTool(tools.Tool):
    """Tool wrapper for sambamba - high-performance SAM/BAM processing"""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [
                tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)
            ]
        super(SambambaTool, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        output = subprocess.check_output(
            [self.install_and_get_path(), '--version'],
            stderr=subprocess.STDOUT
        ).decode('UTF-8')
        # Parse version from output like "sambamba 1.0.1"
        match = re.search(r'sambamba\s+(\d+\.\d+\.\d+)', output)
        if match:
            self.tool_version = match.group(1)
        else:
            # Fallback: try to get first version-like string
            for line in output.split('\n'):
                if 'sambamba' in line.lower():
                    parts = line.split()
                    for part in parts:
                        if re.match(r'\d+\.\d+', part):
                            self.tool_version = part
                            return
            self.tool_version = 'unknown'

    def execute(self, command, args, stdout=None, stderr=None):
        """Execute a sambamba subcommand with arguments"""
        tool_cmd = [self.install_and_get_path(), command] + [str(a) for a in args]
        log.debug(' '.join(tool_cmd))

        stdout_fh = open(stdout, 'w') if stdout else None
        stderr_fh = open(stderr, 'w') if stderr else None

        try:
            subprocess.check_call(tool_cmd, stdout=stdout_fh, stderr=stderr_fh)
        finally:
            if stdout_fh:
                stdout_fh.close()
            if stderr_fh:
                stderr_fh.close()

    def sort(self, inBam, outBam, sort_order='coordinate', threads=None):
        """
        Sort a BAM file.

        Args:
            inBam: Input BAM file path
            outBam: Output sorted BAM file path
            sort_order: 'coordinate' (default) or 'queryname'
            threads: Number of threads to use (default: auto)
        """
        threads = util.misc.sanitize_thread_count(threads)

        args = []

        # Thread count
        args.extend(['-t', str(threads)])

        # Sort order: -n for queryname, default is coordinate
        if sort_order == 'queryname':
            args.append('-n')
        elif sort_order != 'coordinate':
            raise ValueError(f"Invalid sort_order: {sort_order}. Must be 'coordinate' or 'queryname'")

        # Output file
        args.extend(['-o', outBam])

        # Input file
        args.append(inBam)

        self.execute('sort', args)

    def index(self, inBam, threads=None):
        """
        Create an index (.bai) for a coordinate-sorted BAM file.

        Args:
            inBam: Input BAM file path (must be coordinate-sorted)
            threads: Number of threads to use (default: auto)
        """
        threads = util.misc.sanitize_thread_count(threads)

        args = ['-t', str(threads), inBam]

        self.execute('index', args)

    def merge(self, inBams, outBam, threads=None):
        """
        Merge multiple sorted BAM files into a single BAM file.

        Args:
            inBams: List of input BAM file paths (must be coordinate-sorted)
            outBam: Output merged BAM file path
            threads: Number of threads to use (default: auto)
        """
        if not inBams or len(inBams) < 2:
            raise ValueError("merge requires at least 2 input BAM files")

        threads = util.misc.sanitize_thread_count(threads)

        args = ['-t', str(threads), outBam] + list(inBams)

        self.execute('merge', args)

    def flagstat(self, inBam, threads=None):
        """
        Get alignment statistics from a BAM file.

        Args:
            inBam: Input BAM file path
            threads: Number of threads to use (default: auto)

        Returns:
            Dictionary with parsed flagstat statistics
        """
        threads = util.misc.sanitize_thread_count(threads)

        args = ['-t', str(threads), inBam]

        # Capture output
        tool_cmd = [self.install_and_get_path(), 'flagstat'] + [str(a) for a in args]
        log.debug(' '.join(tool_cmd))

        output = subprocess.check_output(tool_cmd).decode('UTF-8')

        return self._parse_flagstat(output)

    def _parse_flagstat(self, output):
        """
        Parse sambamba flagstat output into a dictionary.

        The output format is similar to samtools flagstat:
        123 + 0 in total (QC-passed reads + QC-failed reads)
        0 + 0 secondary
        ...
        """
        stats = {}

        for line in output.strip().split('\n'):
            if not line.strip():
                continue

            # Parse lines like "123 + 0 in total..."
            match = re.match(r'(\d+)\s*\+\s*(\d+)\s+(.*)', line)
            if match:
                passed = int(match.group(1))
                failed = int(match.group(2))
                description = match.group(3).strip()

                # Create clean key from description
                if 'in total' in description:
                    key = 'total'
                elif 'secondary' in description:
                    key = 'secondary'
                elif 'supplementary' in description:
                    key = 'supplementary'
                elif 'duplicates' in description:
                    key = 'duplicates'
                elif 'mapped' in description and 'mate' not in description:
                    key = 'mapped'
                elif 'paired in sequencing' in description:
                    key = 'paired'
                elif 'read1' in description:
                    key = 'read1'
                elif 'read2' in description:
                    key = 'read2'
                elif 'properly paired' in description:
                    key = 'properly_paired'
                elif 'with itself and mate mapped' in description:
                    key = 'both_mapped'
                elif 'singletons' in description:
                    key = 'singletons'
                elif 'mate mapped to a different chr' in description:
                    if 'mapQ>=5' in description:
                        key = 'mate_diff_chr_mapq5'
                    else:
                        key = 'mate_diff_chr'
                else:
                    # Fallback: use cleaned description as key
                    key = re.sub(r'[^a-z0-9]+', '_', description.lower()).strip('_')

                stats[key] = passed
                stats[f'{key}_qcfailed'] = failed

        return stats

    def markdup(self, inBam, outBam, remove_duplicates=False, threads=None,
                hash_table_size=None, overflow_list_size=None):
        """
        Mark or remove duplicate reads using sambamba markdup.

        Args:
            inBam: Input BAM file (must be coordinate-sorted)
            outBam: Output BAM file
            remove_duplicates: If True, remove duplicates; if False, mark them (default: False)
            threads: Number of threads to use (default: auto)
            hash_table_size: Size of hash table for finding read pairs (optional)
            overflow_list_size: Size of overflow list (optional)
        """
        threads = util.misc.sanitize_thread_count(threads)

        args = []

        # Remove duplicates flag
        if remove_duplicates:
            args.append('-r')

        # Thread count
        args.extend(['-t', str(threads)])

        # Optional tuning parameters
        if hash_table_size:
            args.extend(['--hash-table-size', str(hash_table_size)])
        if overflow_list_size:
            args.extend(['--overflow-list-size', str(overflow_list_size)])

        # Input and output files
        args.append(inBam)
        args.append(outBam)

        self.execute('markdup', args)
