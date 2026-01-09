'''
    Tool wrapper for the KMC kmer counter
    ( http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about )
'''

import logging
import os
import os.path
import sys
import subprocess
import shlex
import argparse
import shutil
import re
import collections

import Bio.SeqIO

import read_utils
import tools
import tools.samtools
import util.file
import util.misc


TOOL_NAME = 'kmc'
TOOL_VERSION = '3.1.1rc1'

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

_chk = util.misc.chk

DEFAULT_KMER_SIZE = 25
DEFAULT_COUNTER_CAP = 255         # default for "cap counter values at" params

class KmcTool(tools.Tool):
    '''Tool wrapper for KMC kmer counter'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=False)]

        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def _kmer_db_name(self, kmer_db):
        """Return the kmer database path, given either the db path of the file name of the .kmc_pre or .kmc_suf file"""
        base, ext = os.path.splitext(kmer_db)
        return base if ext in ('.kmc_pre', '.kmc_suf') else kmer_db

    def is_kmer_db(self, kmer_db):
        """Quickly check that the given kmer database exists"""
        kmer_db = self._kmer_db_name(kmer_db)
        return os.path.isfile(kmer_db + '.kmc_pre') and os.path.isfile(kmer_db + '.kmc_suf')

    def _get_file_format_opt(self, fname):
        """Get the KMC command-line option to specify file format"""
        file_type = util.file.uncompressed_file_type(fname)
        if file_type in ('.fasta', '.fa'):
            return 'm'
        elif file_type in ('.fastq', '.fq'):
            return 'q'
        elif file_type == '.bam':
            return 'bam'
        else:
            raise RuntimeError('Unknown seq file format: {}'.format(file_type))

    def build_kmer_db(self, seq_files, kmer_db, kmer_size=DEFAULT_KMER_SIZE, min_occs=1, max_occs=util.misc.MAX_INT32,
                      counter_cap=DEFAULT_COUNTER_CAP, single_strand=False, mem_limit_gb=8, mem_limit_laxness=0,
                      threads=None, kmc_opts=''):
        """Build a database of kmers occurring in the given sequence files.

        Inputs:
          seq_files: filename(s) of files from which to gather kmers.  Files can be fasta, fastq or bam; fasta and fastq
             may be compressed with gzip or bzip2.

        Outputs:
          kmer_db: kmer database path, with or without the .kmc_pre/.kmc_suf extension.

        Params:
          kmer_size: kmer size

          min_occs: ignore kmers occurring fewer than this many times
          max_occs: ignore kmers occurring more than this many times

          counter_cap: when writing kmer counts to the database, cap the values at this number

          single_strand: if True, do not include kmers from reverse complements of input sequences

          mem_limit_gb: max RAM to use when building kmer database.  Might be exceeded unless `mem_limit_strict`
             is True.
          mem_limit_laxness: how strictly to adhere to `mem_limit_gb`: 0=strictly, higher values permit more laxness

          kmc_opts: any additional kmc flags
        """
        seq_files = util.misc.make_seq(seq_files)
        kmer_db = self._kmer_db_name(kmer_db)
        threads = util.misc.sanitize_thread_count(threads)
        with util.file.tmp_dir(suffix='_kmc_db') as t_dir, \
             util.file.tempfname(suffix='_kmc_seq_files') as seq_file_list:
            util.file.dump_file(seq_file_list, '\n'.join(seq_files))
            args = ['-v']

            input_fmt_opts_set = set(map(self._get_file_format_opt, seq_files))
            _chk(len(input_fmt_opts_set) == 1, "All input sequence files must be in the same format")
            input_fmt_opt = list(input_fmt_opts_set)[0]

            args += '-f{} -k{} -ci{} -cx{} -cs{} -m{} -t{}'.format(input_fmt_opt, kmer_size, min_occs, max_occs,
                                                                   counter_cap, mem_limit_gb, threads).split()
            if single_strand:
                args += ['-b']
            if mem_limit_laxness == 0:
                args += ['-sm']  # strict memory mode
            elif mem_limit_laxness == 2:
                args += ['-r'] # RAM-only mode

            if kmc_opts:
                args += shlex.split(kmc_opts)

            args += ['@'+seq_file_list, kmer_db, t_dir]
            tool_cmd = [self.install_and_get_path()] + args
            _log.info('Building kmer database with command: %s', ' '.join(tool_cmd))
            subprocess.check_call(tool_cmd)
            _chk(os.path.isfile(kmer_db+'.kmc_pre') and os.path.isfile(kmer_db+'.kmc_suf'),
                'kmer database files not created: {}'.format(kmer_db))

    def execute(self, args, threads=None, return_output=False):  # pylint: disable=arguments-differ
        """Run kmc_tools with the given args.  If `return_output` is True, returns the stdout output."""
        threads = util.misc.sanitize_thread_count(threads)
        tool_cmd = [self.install_and_get_path()+'_tools'] + ['-t{}'.format(threads)] + list(map(str, args))
        _log.info('Running kmc_tools command: %s', ' '.join(tool_cmd))

        result = None
        if not return_output:
            subprocess.check_call(tool_cmd)
        else:
            result = subprocess.check_output(tool_cmd).decode('utf-8')
        _log.info('Done running kmc_tools command: %s', ' '.join(tool_cmd))
        return result


    def dump_kmer_counts(self, kmer_db, out_kmers, min_occs=1, max_occs=util.misc.MAX_INT32, threads=None):
        """Dump the kmers from the database, with their counts, to a text file.  Kmers with counts
        below `min_occs` or above `max_occs` are ignored."""
        self.execute('transform {} -ci{} -cx{} dump {}'.format(self._kmer_db_name(kmer_db),
                                                               min_occs, max_occs, out_kmers).split(),
                     threads=threads)
        _chk(os.path.isfile(out_kmers), 'dump_kmer_counts: output file not created')

    def read_kmer_counts(self, kmer_counts_txt):
        """Read kmer counts from a file written by dump_kmer_counts, as collections.Counter"""
        counts = collections.Counter()
        with util.file.open_or_gzopen(kmer_counts_txt) as kmers_f:
            for line in kmers_f:
                kmer, count = line.strip().split()
                counts[kmer] = int(count)
        return counts

    def get_kmer_counts(self, kmer_db, **kwargs):
        """Extract and return the kmer counts from the kmer database as a dict.  `kwargs` are passed through
        to read_kmer_counts()."""
        with util.file.tempfname(suffix='_kmer_cnts.txt') as kmer_counts_file:
            self.dump_kmer_counts(kmer_db, out_kmers=kmer_counts_file, **kwargs)
            return self.read_kmer_counts(kmer_counts_file)

    def get_kmer_db_info(self, kmer_db):
        """Return params of a kmer db.
        This functionality is not documented in KMC docs but is supported
        ( https://github.com/refresh-bio/KMC/issues/83 )

        Returns: an argparse.Namespace() with attributes kmer_size, min_occs, max_occs,
           counter_size_bytes, and total_kmers.
        """
        output = self.execute(['info', self._kmer_db_name(kmer_db)], return_output=True, threads=1)
        db_info = dict(re.split('\\s+:\\s+', line.strip()) for line in output.strip().split('\n'))
        def fix_val(v):
            return v=='yes' if v in ('yes', 'no') else util.misc.as_type(v, (int, str))
        db_info = {k: fix_val(v) for k, v in db_info.items()}

        return argparse.Namespace(kmer_size=db_info['k'], single_strand=not db_info['both strands'],
                                  min_occs=db_info['cutoff min'],
                                  max_occs=db_info['cutoff max'],
                                  counter_size_bytes=int(db_info['counter size'].split()[0]),
                                  total_kmers=db_info['total k-mers'])

    def filter_reads(self, kmer_db, in_reads, out_reads, db_min_occs=1, db_max_occs=util.misc.MAX_INT32,
                     read_min_occs=0, read_max_occs=util.misc.MAX_INT32,
                     read_min_occs_frac=0.0, read_max_occs_frac=1.0, hard_mask=False, threads=None):
        """Filter reads based on their kmer contents.

        Note that "occurrence of a kmer" means "occurrence of the kmer or its reverse complement" if kmer_db was built
        with single_strand==False.

        Inputs:
          kmer_db: the kmc kmer database
          in_reads: the reads to filter.  can be a .fasta or .fastq or .bam; fasta or fastq can be compressed
             with gzip or bzip2.  If a .bam, a read pair is kept if either mate passes the filter.

        Outputs:
          out_reads: file to which filtered reads are written.  type is determined from extension,
             same types as above are supported.

        Params:
          db_min_occs: only consider database kmers with at least this count
          db_max_occs: only consider database kmers with at most this count

          read_min_occs: only keep reads with at least this many occurrences of kmers from database.
          read_max_occs: only keep reads with no more than this many occurrence of kmers from the database.
          read_min_occs_frac: only keep reads with at least this many occurrences of kmers from database,
             interpreted as a fraction of read length in kmers
          read_max_occs_frac: only keep reads with no more than this many occurrence of kmers from the database.
             interpreted as a fraction of read length in kmers.

          (Note: minimal read kmer content can be given as absolute counts or fraction of read length, but not both).

          hard_mask: if True, in the output reads, kmers not passing the filter are replaced by Ns
          threads: use this many threads
        """
        _log.debug('FILTER_READS: locals=%s dbinfo=%s', locals(), self.get_kmer_db_info(kmer_db))

        abs_thresholds = (read_min_occs, read_max_occs) != (0, util.misc.MAX_INT32)
        rel_thresholds = (read_min_occs_frac, read_max_occs_frac) != (0., 1.)
        if not (abs_thresholds or rel_thresholds): abs_thresholds = True
        _chk(not (abs_thresholds and rel_thresholds),
             "Mixing absolute and relative thresholds for kmer content not currently supported")
        _chk(0 <= read_min_occs <= read_max_occs <= util.misc.MAX_INT32,
             'Invalid kmer contents thresholds')
        _chk(0.0 <= read_min_occs_frac <= read_max_occs_frac <= 1.0,
             'Invalid kmer contents thresholds')

        in_reads_type = util.file.uncompressed_file_type(in_reads)
        _in_reads = in_reads
        _out_reads = out_reads
        with util.file.tmp_dir(suffix='_kmcfilt') as t_dir:
            if in_reads_type in ('.fa', '.fasta'):
                # kmc_tools filter currently requires fasta files to be in fasta-2line format
                # https://github.com/refresh-bio/KMC/issues/57
                _in_reads = os.path.join(t_dir, 'in_reads.fasta')
                Bio.SeqIO.convert(in_reads, 'fasta', _in_reads, 'fasta-2line')
            if in_reads_type == '.bam':
                # kmc_tools filter currently does not support .bam files
                # https://github.com/refresh-bio/KMC/issues/66
                _in_reads = os.path.join(t_dir, 'in_reads.fasta')
                tools.samtools.SamtoolsTool().bam2fa(in_reads, _in_reads)
                _out_reads = os.path.join(t_dir, 'out_reads.fasta')

                # TODO: if db is single-strand, reverse-complement read2's?

            in_reads_fmt = 'q' if in_reads_type in ('.fq', '.fastq') else 'a'

            if os.path.getsize(self._kmer_db_name(kmer_db)+'.kmc_suf') == 8:
                _chk(self.get_kmer_db_info(kmer_db).total_kmers == 0,
                     'Something odd: did kmc file format change?')
                # kmc has a bug filtering on empty kmer databases:
                # https://github.com/refresh-bio/KMC/issues/86
                if read_min_occs > 0 or read_min_occs_frac > 0.0:
                    util.file.make_empty(_out_reads)
                else:
                    shutil.copyfile(_in_reads, _out_reads)
            else:
                thresholds_opt = ' -ci{} -cx{}'.format(*((read_min_occs, read_max_occs) if abs_thresholds else \
                                                         (read_min_occs_frac, read_max_occs_frac)))
                self.execute('filter {} {} -ci{} -cx{} {} -f{} {} {}'.format('-hm' if hard_mask else '',
                                                                             self._kmer_db_name(kmer_db),
                                                                             db_min_occs, db_max_occs,
                                                                             _in_reads, in_reads_fmt,
                                                                             thresholds_opt,
                                                                             _out_reads).split(),
                             threads=threads)

            if in_reads_type == '.bam':
                _chk(out_reads.endswith('.bam'), 'output from .bam to non-.bam not yet supported')
                passing_read_names = os.path.join(t_dir, 'passing_read_names.txt')
                read_utils.fasta_read_names(_out_reads, passing_read_names)
                # Load passing read names into ReadIdStore and filter BAM (keep matching reads)
                db_path = os.path.join(t_dir, 'read_ids.db')
                with read_utils.ReadIdStore(db_path) as store:
                    store.add_from_readlist(passing_read_names)
                    store.filter_bam_by_ids(in_reads, out_reads, include=True)
        # end: with util.file.tmp_dir(suffix='kmcfilt') as t_dir
    # end: def filter_reads(self, kmer_db, in_reads, out_reads, db_min_occs=1, db_max_occs=util.misc.MAX_INT32, ...)

    def kmers_binary_op(self, op, kmer_db1, kmer_db2, kmer_db_out,
                        result_min_occs=1, result_max_occs=util.misc.MAX_INT32,
                        result_counter_cap=DEFAULT_COUNTER_CAP,
                        threads=None):
        """Perform a simple binary operation on two kmer sets"""
        kmer_db1, kmer_db2, kmer_db_out = map(self._kmer_db_name, (kmer_db1, kmer_db2, kmer_db_out))
        # db1_min_occs, db1_max_occs, db2_min_occs, db2_max_occs, db_out_min_occs, db_out_max_occs,
        self.execute(['simple', kmer_db1, kmer_db2, op, kmer_db_out,
                      '-ci{}'.format(result_min_occs),
                      '-cx{}'.format(result_max_occs),
                      '-cs{}'.format(result_counter_cap)], threads=threads)
        _chk(self.is_kmer_db(kmer_db_out), 'kmer_binary_op: output not created')

    def set_kmer_counts(self, kmer_db_in, value, kmer_db_out, threads=None):
        """Create a copy of the kmer database with all counts set to specified value"""
        _chk(1 <= value <= util.misc.MAX_INT32, 'can only set kmer counts to a positive 32-bit int')
        kmer_db_in, kmer_db_out = map(self._kmer_db_name, (kmer_db_in, kmer_db_out))
        # db1_min_occs, db1_max_occs, db2_min_occs, db2_max_occs, db_out_min_occs, db_out_max_occs,
        self.execute(['transform', kmer_db_in, 'set_counts', value, kmer_db_out], threads=threads)
        _chk(self.is_kmer_db(kmer_db_out), 'set_kmer_counts: output not created')

# end: class KmcTool(tools.Tool)
