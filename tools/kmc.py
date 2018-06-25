'''
    Tool wrapper for the KMC kmer counter 
    ( http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about )
'''

import logging
import os
import os.path
import subprocess
import shutil
import shlex

import read_utils
import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

import Bio.SeqIO

TOOL_NAME = 'kmc'
TOOL_VERSION = '3.1.0'

log = logging.getLogger(__name__)


DEFAULT_KMER_SIZE=27
DEFAULT_COUNTER_CAP=255         # default for "cap counter values at" params

class KmcTool(tools.Tool):
    '''Tool wrapper for KMC kmer counter'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='kmc',
                                                  verifycmd='kmc -h > /dev/null')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    @staticmethod
    def _kmer_db_name(kmer_db):
        """Return the kmer database path, given either the db path of the file name of the .kmc_pre or .kmc_suf file"""
        base, ext = os.path.splitext(kmer_db)
        return base if ext in ('.kmc_pre', '.kmc_suf') else kmer_db

    @staticmethod
    def _get_file_format_opt(fname):
        """Get the KMC command-line option to specify file format"""
        file_type = util.file.uncompressed_file_type(fname)
        if file_type in ('.fasta', '.fa'): return 'm'
        elif file_type in ('.fastq', '.fq'): return 'q'
        elif file_type == '.bam': return 'bam'
        else: raise RuntimeError('Unknown seq file format: {}'.format(file_type))

    def build_kmer_db(self, seq_files, kmer_db, kmer_size=DEFAULT_KMER_SIZE, min_occs=None, max_occs=None, 
                      counter_cap=DEFAULT_COUNTER_CAP, mem_limit_gb=8, threads=None, kmc_opts=''):
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
          kmc_opts: any additional kmc flags

          
        """
        if min_occs is None: min_occs = 1
        if max_occs is None: max_occs = util.misc.MAX_INT32
        seq_files = util.misc.make_seq(seq_files)
        kmer_db = self._kmer_db_name(kmer_db)
        threads = util.misc.sanitize_thread_count(threads)
        with util.file.tmp_dir(suffix='kmcdb') as t_dir, util.file.tempfname(suffix='kmcfiles') as seq_file_list:
            util.file.dump_file(seq_file_list, '\n'.join(seq_files))
            args = ['-v']

            input_fmt_opts = set(map(self._get_file_format_opt, seq_files))
            assert len(input_fmt_opts) == 1, "All input files must be of the same format"
            input_fmt_opt = list(input_fmt_opts)[0]

            args += '-f{} -k{} -ci{} -cx{} -cs{} -m{} -r -t{} @{} {}'.format(input_fmt_opt, kmer_size, min_occs, max_occs, counter_cap, mem_limit_gb, threads,
                                                                              seq_file_list, kmer_db).split()
            if kmc_opts: args += shlex.split(kmc_opts)
            args += [t_dir]
            tool_cmd = [self.install_and_get_path()] + args
            log.info('Building kmer database with command: ' + ' '.join(tool_cmd))
            subprocess.check_call(tool_cmd)
            assert os.path.isfile(kmer_db+'.kmc_pre') and os.path.isfile(kmer_db+'.kmc_suf'), \
                'kmer database files not created: {}'.format(kmer_db)

    def execute(self, args, threads=None):
        """Run kmc_tools with the given args"""
        threads = util.misc.sanitize_thread_count(threads)
        tool_cmd = [self.install_and_get_path()+'_tools'] + ['-v', '-t{}'.format(threads)] + list(map(str, args))
        log.info('Running kmc_tools command: ' + ' '.join(tool_cmd))
        print('Running kmc_tools command: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)


    def dump_kmers(self, kmer_db, out_kmers, min_occs=None, max_occs=None, threads=None):
        """Dump the kmers from the database to a text file"""
        if min_occs is None: min_occs = 1
        if max_occs is None: max_occs = util.misc.MAX_INT32
        self.execute('transform {} -ci{} -cx{} dump -s {}'.format(self._kmer_db_name(kmer_db), min_occs, max_occs, out_kmers).split(), threads=threads)
        assert os.path.isfile(out_kmers)


    def filter_reads(self, kmer_db, in_reads, out_reads, db_min_occs=None, db_max_occs=None, 
                     read_min_occs=None, read_max_occs=None, hard_mask=False, threads=None):
        """Filter reads based on their kmer contents.

        Note that "occurrence of a kmer" means "occurrence of the kmer or its reverse complement".

        Inputs:
          kmer_db: the kmc kmer database
          in_reads: the reads to filter.  can be a .fasta or .fastq or .bam; fasta or fastq can be compressed with gzip or bzip2.
             If a .bam, a read pair is kept if either mate passes the filter.

        Outputs:
          out_reads: file to which filtered reads are written.  type is determined from extension, same types as above are supported.

        Params:
          db_min_occs: only consider database kmers with at least this count
          db_max_occs: only consider database kmers with at most this count
          read_min_occs: only keep reads with at least this many occurrences of kmers from database.  If `as_perc` is True, interpreted as percent
             of read length
          read_max_occs: only keep reads with no more than this many occurrence of kmers from the database
          hard_mask: if True, in the output reads, kmers not passing the filter are replaced by Ns
        """

        if db_min_occs is None: db_min_occs=1
        if db_max_occs is None: db_max_occs=util.misc.MAX_INT32

        if read_min_occs is None and read_max_occs is None:
            read_min_occs, read_max_occs = 1, util.misc.MAX_INT32
        elif read_min_occs is None and read_max_occs is not None:
            read_min_occs = 0.0 if isinstance(read_max_occs, float) else 1
        elif read_max_occs is None and read_min_occs is not None:
            read_max_occs = 1.0 if isinstance(read_min_occs, float) else util.misc.MAX_INT32

        assert type(read_min_occs) == type(read_max_occs), 'read_min_occs and read_max_occs must be specified the same way (as kmer count or fraction of read length)'
        assert read_min_occs <= read_max_occs, 'vals are {} {}'.format(read_min_occs, read_max_occs)
        assert not isinstance(read_min_occs, float) or 0.0 <= read_min_occs <= read_max_occs <= 1.0

        in_reads_type = util.file.uncompressed_file_type(in_reads)
        _in_reads = in_reads
        _out_reads = out_reads
        with util.file.tmp_dir(suffix='kmcfilt') as t_dir:
            if in_reads_type in ('.fa', '.fasta'):
                # kmc_tools filter currently requires fasta files to be in fasta-2line format
                # https://github.com/refresh-bio/KMC/issues/57
                _in_reads = os.path.join(t_dir, 'in_reads.fasta')
                with util.file.open_or_gzopen(in_reads, 'rt'):
                    Bio.SeqIO.convert(in_reads, 'fasta', _in_reads, 'fasta-2line')
            if in_reads_type == '.bam':
                # kmc_tools filter currently does not support .bam files
                # https://github.com/refresh-bio/KMC/issues/66
                _in_reads = os.path.join(t_dir, 'in_reads.fasta')
                tools.samtools.SamtoolsTool().bam2fa(in_reads, _in_reads, add_mate_num=True)
                _out_reads = os.path.join(t_dir, 'out_reads.fasta')

            in_reads_fmt = 'q' if in_reads_type in ('.fq', '.fastq') else 'a'

            self.execute('filter {} {} -ci{} -cx{} {} -f{} -ci{} -cx{} {}'.format('-hm' if hard_mask else '',
                                                                                  self._kmer_db_name(kmer_db), db_min_occs, db_max_occs, _in_reads, in_reads_fmt,
                                                                                  read_min_occs, read_max_occs, _out_reads).split(), threads=threads)

            if in_reads_type == '.bam':
                assert out_reads.endswith('.bam')
                passed_read_names = os.path.join(t_dir, 'passed_read_names.txt')
                read_utils.fasta_read_names(_out_reads, passed_read_names)
                shutil.copyfile(passed_read_names, '/tmp/passed.txt')
                tools.picard.FilterSamReadsTool().execute(inBam=in_reads, exclude=False, readList=passed_read_names, outBam=out_reads)
        # end: with util.file.tmp_dir(suffix='kmcfilt') as t_dir
    # end: def filter_reads(self, kmer_db, in_reads, out_reads, db_min_occs=1, db_max_occs=util.misc.MAX_INT32, read_min_occs=None, read_max_occs=None, threads=None)

    def kmers_binary_op(self, op, kmer_db1, kmer_db2, kmer_db_out, threads=None):
        """Perform a simple binary operation on two kmer sets"""
        kmer_db1, kmer_db2, kmer_db_out = map(self._kmer_db_name, (kmer_db1, kmer_db2, kmer_db_out))
        # db1_min_occs, db1_max_occs, db2_min_occs, db2_max_occs, db_out_min_occs, db_out_max_occs,
        self.execute(['simple', kmer_db1, kmer_db2, op, kmer_db_out], threads=threads)

    def set_kmer_counts(self, kmer_db_in, value, kmer_db_out, threads=None):
        """Create a copy of the kmer database with all counts set to specified value"""
        assert 1 <= value <= util.misc.MAX_INT32, 'can only set kmer counts to a positive 32-bit int'
        kmer_db_in, kmer_db_out = map(self._kmer_db_name, (kmer_db_in, kmer_db_out))
        # db1_min_occs, db1_max_occs, db2_min_occs, db2_max_occs, db_out_min_occs, db_out_max_occs,
        self.execute(['transform', kmer_db_in, 'set_counts', value, kmer_db_out], threads=threads)

# end: class KmcTool(tools.Tool)

