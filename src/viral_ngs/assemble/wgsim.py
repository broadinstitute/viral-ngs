'''
    WGSIM - read simulator for next-generation sequencing
    https://github.com/lh3/wgsim
'''

__author__ = "dpark@broadinstitute.org"

import logging
import math
import tools
import tools.picard
import util.file
import util.misc
import os
import os.path
import shutil
import subprocess
import Bio.SeqIO
import Bio.Seq

TOOL_NAME = "wgsim"

_log = logging.getLogger(__name__)

class WgsimTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(WgsimTool, self).__init__(install_methods=install_methods)

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        # wgsim doesn't have a clean version flag, parse from help output
        try:
            output = subprocess.check_output([self.install_and_get_path()], stderr=subprocess.STDOUT).decode('UTF-8')
            # Look for "Version:" in output
            for line in output.split('\n'):
                if 'Version:' in line:
                    self.tool_version = line.split('Version:')[1].strip()
                    return
            self.tool_version = 'unknown'
        except subprocess.CalledProcessError as e:
            # wgsim returns non-zero when called with no args, but outputs help
            output = e.output.decode('UTF-8')
            for line in output.split('\n'):
                if 'Version:' in line:
                    self.tool_version = line.split('Version:')[1].strip()
                    return
            self.tool_version = 'unknown'

    def slice_fasta(self, in_fasta, out_fasta, seq_id=None, start=None, end=None):
        ''' Slice a fasta file to extract a specific sequence and optional coordinate range.

        Args:
            in_fasta: Input fasta file
            out_fasta: Output fasta file with sliced sequence
            seq_id: Sequence ID to extract (default: None). If None, all sequences are copied.
            start: Start coordinate (1-based, inclusive). None means start of sequence.
            end: End coordinate (1-based, inclusive). None means end of sequence.

        If all optional parameters are None, this is a no-op (output identical to input).
        '''
        # If no filtering/slicing requested, just copy the file
        if seq_id is None and start is None and end is None:
            shutil.copyfile(in_fasta, out_fasta)
            return

        found = False
        with open(out_fasta, 'w') as outf:
            for record in Bio.SeqIO.parse(in_fasta, 'fasta'):
                if seq_id is None or record.id == seq_id:
                    found = True
                    # Convert to 0-based for python slicing
                    start_idx = (start - 1) if start is not None else None
                    end_idx = end if end is not None else None
                    sliced_seq = record.seq[start_idx:end_idx]
                    sliced_record = Bio.SeqRecord.SeqRecord(
                        sliced_seq,
                        id=record.id,
                        description=record.description
                    )
                    Bio.SeqIO.write(sliced_record, outf, 'fasta')
                    if seq_id is not None:
                        break

        if seq_id is not None and not found:
            raise ValueError(f"Sequence ID '{seq_id}' not found in {in_fasta}")

    def coverage_to_read_pairs(self, coverage, sequence_length, read_length):
        ''' Convert coverage depth to number of read pairs needed.

        Args:
            coverage: Desired coverage depth (e.g., 20 for 20X)
            sequence_length: Length of the sequence to cover
            read_length: Length of each read

        Returns:
            Number of read pairs (rounded up to nearest integer)
        '''
        # Each read pair contributes 2 * read_length bases
        total_bases_needed = coverage * sequence_length
        read_pairs = total_bases_needed / (2 * read_length)
        return math.ceil(read_pairs)

    def fastqs_to_bam(self, in_fastq1, in_fastq2, out_bam, sample_name='sample', library_name='lib1',
                      platform='ILLUMINA', platform_unit='wgsim'):
        ''' Convert paired fastq files to unaligned BAM with read group headers.

        Args:
            in_fastq1: Input fastq file for read 1
            in_fastq2: Input fastq file for read 2
            out_bam: Output unaligned BAM file
            sample_name: Sample name for read group
            library_name: Library name for read group
            platform: Platform for read group (default: ILLUMINA)
            platform_unit: Platform unit for read group (default: wgsim)
        '''
        picard = tools.picard.FastqToSamTool()
        picard.execute(
            in_fastq1, in_fastq2, sample_name, out_bam,
            picardOptions=[
                'LIBRARY_NAME=' + library_name,
                'PLATFORM=' + platform,
                'PLATFORM_UNIT=' + platform_unit
            ]
        )

    def execute(self, in_fasta, out_fastq1, out_fastq2,
                read_length=150, outer_distance=500, error_rate=0.02, mutation_rate=0.001,
                indel_fraction=0.15, indel_extended_prob=0.3, num_read_pairs=None,
                random_seed=None, other_opts=()):
        ''' Execute wgsim to generate simulated paired-end reads.

        Args:
            in_fasta: Input reference fasta file
            out_fastq1: Output fastq file for read 1
            out_fastq2: Output fastq file for read 2
            read_length: Length of each read (default: 150 for Illumina)
            outer_distance: Outer distance between the two ends (default: 500)
            error_rate: Base error rate (default: 0.02)
            mutation_rate: Rate of mutations (default: 0.001)
            indel_fraction: Fraction of indels (default: 0.15)
            indel_extended_prob: Probability an indel is extended (default: 0.3)
            num_read_pairs: Number of read pairs to generate (required)
            random_seed: Random seed for reproducibility
            other_opts: Additional command line options
        '''

        if num_read_pairs is None:
            raise ValueError("num_read_pairs must be specified")

        tool_cmd = [
            self.install_and_get_path(),
            '-e', str(error_rate),
            '-d', str(outer_distance),
            '-N', str(num_read_pairs),
            '-1', str(read_length),
            '-2', str(read_length),
            '-r', str(mutation_rate),
            '-R', str(indel_fraction),
            '-X', str(indel_extended_prob),
        ]

        if random_seed is not None:
            tool_cmd.extend(['-S', str(random_seed)])

        tool_cmd.extend(other_opts)
        tool_cmd.extend([in_fasta, out_fastq1, out_fastq2])

        _log.debug('running wgsim: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
