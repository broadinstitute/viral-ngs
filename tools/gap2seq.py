'''
    Gap2Seq - assembly gap closing tool
'''

import logging
import os
import os.path
import subprocess
import shlex
import shutil
import tempfile
import time

import Bio.Seq

import tools
import tools.samtools
import util.file
import util.misc

TOOL_NAME = 'gap2seq'
TOOL_VERSION = '2.1'

log = logging.getLogger(__name__)

class Gap2SeqTool(tools.Tool):
    """Tool wrapper for the Gap2Seq gap-closing tool."""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, env='gap2seq_env', executable='Gap2Seq.sh')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + args
        log.info('b4 gap2seq: ' + str(os.listdir()))
        log.info('running gap2seq: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
        log.info('aft gap2seq: ' + str(os.listdir()))

    def gapfill(self, in_scaffold, inBam, out_scaffold, solid_kmer_thresholds=(3,2), kmer_sizes=(90, 80, 70, 60, 50, 40, 31),
                min_gap_to_close=4, gap2seq_opts='', mem_limit_gb=4, threads=1, time_limit_minutes=60):
        """Try to fill the gaps in the given scaffold, using the reads.

        Inputs:
            in_scaffold: a FASTA file containing the scaffold.  Each FASTA record corresponds to one
                segment (for multi-segment genomes).  Contigs within each segment are
                separated by Ns.  The exact number of Ns between contigs does not matter, as the length of the gap is one 
                of the things determined by the gap-filling tool.  (But see `min_gap_to_close`).
            inBam: reads to use for filling the gaps.  Only paired-end reads from the bam file are used, any unpaired reads
                are ignored.
           
        Outputs:
            out_scaffold: the input scaffold, with some of the gaps between contigs possibly filled.

        Params:
            solid_kmer_thresholds: kmers must appear at least this many times in the reads to be considered solid.
                We try gapfilling for all combinations of values of solid_kmer_thresholds and kmer_sizes.
            kmer_sizes: kmer sizes to use.  We try gapfilling for all combinations of values of solid_kmer_thresholds and kmer_sizes.
            min_gap_to_close: stop gap-closing if all gaps are no longer than this many Ns
            gap2seq_opts: extra command-line flags to pass to Gap2Seq
            mem_limit_gb: max memory to use, in gigabytes
            threads: number of threads to use
            time_limit_minutes: stop trying to close more gaps after this many minutes (currently this is a soft/advisory limit)
        
        """
        from os.path import abspath, join
        solid_kmer_thresholds = sorted(util.misc.make_seq(solid_kmer_thresholds), reverse=True)
        kmer_sizes = sorted(util.misc.make_seq(kmer_sizes), reverse=True)
        samtools = tools.samtools.SamtoolsTool()
        with samtools.bam2fq_tmp(inBam) as (reads1, reads2):
            done = False
            prev_scaffold = in_scaffold
            stop_time = time.time() + 60*time_limit_minutes
            with util.file.tmp_dir('_gap2seq_dir') as gap2seq_dir:
                for solid_kmer_threshold in solid_kmer_thresholds:
                    for kmer_size in kmer_sizes:
                        gap2seq_filled = join(gap2seq_dir, 'gap2seq-filled.s{}.k{}.fasta'.format(solid_kmer_threshold, kmer_size))
                        log.info('s={} k={} gap2seq_filled={}'.format(solid_kmer_threshold, kmer_size, gap2seq_filled))

                        # gap2seq (or rather the gatb library it uses) currently has a bug where tempfiles are left in the
                        # current working directory.  So we run it in its own temp dir, but then must give absolute paths
                        # for all files.

                        args = ['-scaffolds', abspath(prev_scaffold), '-filled', abspath(gap2seq_filled),
                                '-reads', ','.join((abspath(reads1), abspath(reads2))), '-solid', str(solid_kmer_threshold), 
                                '-k', str(kmer_size), '-all-upper',
                                '-verbose', '-nb-cores', str(threads), '-max-mem', str(mem_limit_gb)] + shlex.split(gap2seq_opts)
                        
                        with util.file.tmp_dir('_gap2seq_run_dir') as gap2seq_run_dir:
                            with util.file.tmp_chdir(gap2seq_run_dir):
                                self.execute(args=args)

                        prev_scaffold = gap2seq_filled
                        if not any('N'*min_gap_to_close in str(rec.seq) for rec in Bio.SeqIO.parse(gap2seq_filled, 'fasta')):
                            log.info('no gaps left, quittting gap2seq early')
                            done=True
                            break
                        if time.time() > stop_time:
                            log.info('Time limit for gap closing reached')
                            done=True
                            break
                    # end: for kmer_size in kmer_sizes
                    if done: break
                # end: for solid_kmer_threshold in solid_kmer_thresholds

                shutil.copyfile(src=prev_scaffold, dst=out_scaffold)

