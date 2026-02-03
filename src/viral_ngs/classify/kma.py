'''
KMA k-mer alignment tool
'''
import logging
import os
import os.path
import shutil
import subprocess

from viral_ngs import core
from viral_ngs.core import picard
from viral_ngs.core import samtools
from viral_ngs.core import file
from viral_ngs.core import misc
from builtins import super

log = logging.getLogger(__name__)

class KMA(core.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(core.PrexistingUnixCommand(shutil.which('kma'), require_executability=False))
        super(KMA, self).__init__(install_methods=install_methods)

    def version(self):
        '''Parse version from kma -v output'''
        try:
            result = subprocess.run(['kma', '-v'], capture_output=True, text=True, check=True)
            version_line = result.stdout.strip() or result.stderr.strip()
            if version_line:
                parts = version_line.split()
                if len(parts) >= 2:
                    return parts[1]
            return 'unknown'
        except (subprocess.CalledProcessError, FileNotFoundError, IndexError):
            return 'unknown'

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())

    def execute(self, command, args=None, options=None):
        '''Run a kma command.

        Args:
          command: Command to run (e.g., 'kma', 'kma_index').
          args: List of positional args.
          options: Dict of keyword options.
        '''
        options = options or {}
        args = args or []

        cmd = [command]

        for key, value in options.items():
            if value is None:
                cmd.append(key)
            else:
                cmd.extend([key, str(value)])

        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        subprocess.check_call(cmd)

    def build(self, ref_fasta, db_prefix, num_threads=None):
        '''Build a KMA database from a reference FASTA file.

        Args:
          ref_fasta: Input reference FASTA file.
          db_prefix: Output database prefix.
          num_threads: Number of threads to use.
        '''
        opts = {
            '-i': ref_fasta,
            '-o': db_prefix
        }
        if num_threads:
            opts['-t'] = misc.sanitize_thread_count(num_threads)

        self.execute('kma_index', options=opts)

    def classify(self, in_bam, db, out_prefix, num_threads=None):
        '''Classify reads from BAM file using KMA.

        Args:
          in_bam: Input unaligned reads in BAM format.
          db: KMA database prefix (without file extensions).
          out_prefix: Output prefix for KMA result files.
          num_threads: Number of threads to use.
        '''
        if samtools.SamtoolsTool().isEmpty(in_bam):
            log.warning("Input BAM is empty, skipping KMA classification")
            with open(out_prefix + '.res', 'wt') as outf:
                pass
            return

        tmp_fastq1 = file.mkstempfname('.1.fastq')
        tmp_fastq2 = file.mkstempfname('.2.fastq')
        tmp_fastq3 = file.mkstempfname('.s.fastq')

        try:
            picard = picard.SamToFastqTool()
            picard_opts = {
                'CLIPPING_ATTRIBUTE': picard.SamToFastqTool.illumina_clipping_attribute,
                'CLIPPING_ACTION': 'X'
            }
            picard.execute(in_bam, tmp_fastq1, tmp_fastq2, outFastq0=tmp_fastq3,
                           picardOptions=picard.PicardTools.dict_to_picard_opts(picard_opts),
                           JVMmemory=picard.jvmMemDefault)

            opts = {
                '-t_db': db,
                '-o': out_prefix
            }
            if num_threads:
                opts['-t'] = misc.sanitize_thread_count(num_threads)

            # SamToFastq outputs paired reads to tmp_fastq1/2 and unpaired to tmp_fastq3.
            # For true paired-end data, tmp_fastq2 will be larger than tmp_fastq3.
            # For single-end data, tmp_fastq2 is empty/tiny and tmp_fastq3 contains all reads.
            if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
                log.warning("running in single-end read mode")
                opts['-i'] = tmp_fastq3
            else:
                opts['-ipe'] = ' '.join([tmp_fastq1, tmp_fastq2])

            self.execute('kma', options=opts)

        finally:
            if os.path.exists(tmp_fastq1):
                os.unlink(tmp_fastq1)
            if os.path.exists(tmp_fastq2):
                os.unlink(tmp_fastq2)
            if os.path.exists(tmp_fastq3):
                os.unlink(tmp_fastq3)
