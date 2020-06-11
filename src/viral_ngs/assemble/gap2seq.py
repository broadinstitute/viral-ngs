'''
    Gap2Seq - assembly gap closing tool
'''

import itertools
import functools
import operator
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

TOOL_NAME = 'Gap2Seq'

log = logging.getLogger(__name__)

class Gap2SeqTool(tools.Tool):
    """Tool wrapper for the Gap2Seq gap-closing tool."""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(Gap2SeqTool, self).__init__(install_methods=install_methods)

    def version(self):
        return None

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(args)
        log.debug('running gap2seq: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def _run_gap2seq(self, reads, scaffolds, filled, *args, **kwargs):
        # gap2seq (or rather the gatb library it uses) currently has a bug where tempfiles are left in the
        # current working directory.  So we run it in its own temp dir, but then must give absolute paths
        # for all files.

        abspath = os.path.abspath
        file_args = ('-s', abspath(scaffolds), '--filled', abspath(filled), 
                     '--reads', ','.join(map(abspath,reads)))
        more_args = functools.reduce(operator.concat, 
                                     [(('--' if len(arg) > 1 else '-') + arg.replace('_','-'), str(val))
                                      for arg, val in kwargs.items()], ())
        with util.file.tmp_dir('_gap2seq_run_dir') as gap2seq_run_dir:
            with util.file.pushd_popd(gap2seq_run_dir):
                self.execute(file_args+args+more_args)

    def gapfill(self, in_scaffold, in_bam, out_scaffold, solid_kmer_thresholds=(3,), kmer_sizes=(90, 80, 70, 60, 50, 40, 31),
                min_gap_to_close=4, gap2seq_opts='', mem_limit_gb=4.0, threads=None, time_soft_limit_minutes=60.0, random_seed=0):
        """Try to fill the gaps in the given scaffold, using the reads.

        Inputs:
            in_scaffold: a FASTA file containing the scaffold.  Each FASTA record corresponds to one
                segment (for multi-segment genomes).  Contigs within each segment are
                separated by Ns.  The exact number of Ns between contigs does not matter, as the length of the gap is one 
                of the things determined by the gap-filling tool.  (But see `min_gap_to_close`).
            in_bam: reads to use for filling the gaps.  Only paired-end reads from the bam file are used, any unpaired reads
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
            threads: number of threads to use; None means use all available cores.
            time_soft_limit_minutes: stop trying to close more gaps after this many minutes (currently this is a soft/advisory limit)
            random_seed: random seed for choosing random paths (0 to use current time)
        
        """
        solid_kmer_thresholds = sorted(util.misc.make_seq(solid_kmer_thresholds), reverse=True)
        kmer_sizes = sorted(util.misc.make_seq(kmer_sizes), reverse=True)
        stop_time = time.time() + 60*time_soft_limit_minutes
        threads = util.misc.sanitize_thread_count(threads, tool_max_cores_value=0)
        util.misc.chk(out_scaffold != in_scaffold)
        with tools.samtools.SamtoolsTool().bam2fq_tmp(in_bam) as reads, util.file.tmp_dir('_gap2seq_dir') as gap2seq_dir:

            # We call Gap2Seq for a range of parameter combinations.  Output of each call is input to the next call, so
            # each call only deals with gaps not closed by prior calls.  We first try to close using higher-quality kmers,
            # and if that fails try with lower-quality ones.
            prev_scaffold = in_scaffold
            for solid_kmer_threshold, kmer_size in itertools.product(solid_kmer_thresholds, kmer_sizes):

                if not any('N'*min_gap_to_close in str(rec.seq) for rec in Bio.SeqIO.parse(prev_scaffold, 'fasta')):
                    log.info('no gaps left, quittting gap2seq early')
                    break
                if time.time() > stop_time:
                    log.info('Time limit for gap closing reached')
                    break

                filled_scaffold = os.path.join(gap2seq_dir, 'gap2seq-filled.s{}.k{}.fasta'.format(solid_kmer_threshold, kmer_size))
                self._run_gap2seq(reads, prev_scaffold, filled_scaffold,
                                  *(['--all-upper']+shlex.split(gap2seq_opts)),
                                  solid=solid_kmer_threshold, k=kmer_size, threads=threads,
                                  max_mem=mem_limit_gb, randseed=random_seed)

                prev_scaffold = filled_scaffold

            Bio.SeqIO.convert(prev_scaffold, 'fasta', out_scaffold, 'fasta')

