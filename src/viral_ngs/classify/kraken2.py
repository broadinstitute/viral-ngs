'''
KRAKEN metagenomics classifier
'''
import collections
import itertools
import logging
import os
import os.path
import shutil
import subprocess
import sys
import tempfile

import tools
import tools.picard
import tools.samtools
import util.file
import util.misc
from builtins import super

log = logging.getLogger(__name__)

class Kraken2(tools.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(tools.PrexistingUnixCommand('kraken2', verifycmd=False, require_executability=False))
        super(Kraken2, self).__init__(install_methods=install_methods)

    def version(self):
        return '2'

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())

    def build(self, db, standard_libraries=(), custom_libraries=(), taxdump_out=None,
        k=None, l=None, s=None, protein=False, max_db_size=None, num_threads=None):
        '''Create a kraken database.

        Args:
          db: Kraken database directory to build. Must have library/ and
            taxonomy/ subdirectories to build from.
          standard_libraries: a list of strings from the set of:
             "archaea", "bacteria", "plasmid",
             "viral", "human", "fungi", "plant", "protozoa",
             "nr", "nt", "env_nr", "env_nt", "UniVec",
             "UniVec_Core"
          *args: List of input filenames to process.
        '''

        self.execute('kraken2-build', db, None, options={
            '--threads':util.misc.sanitize_thread_count(num_threads),
            '--download-taxonomy': None
            })

        # add standard libraries:
        for lib in standard_libraries:
            self.execute('kraken2-build', db, None, options={
                '--threads':util.misc.sanitize_thread_count(num_threads),
                '--download-library': lib
                })
        for lib in custom_libraries:
            self.execute('kraken2-build', db, None, options={
                '--threads':util.misc.sanitize_thread_count(num_threads),
                '--add-to-library': lib
                })

        # build db
        build_opts = {
            '--build': None,
            '--threads': util.misc.sanitize_thread_count(num_threads)
        }
        if k:
            build_opts['--kmer-len'] = k
        if l:
            build_opts['--minimizer-len'] = l
        if s:
            build_opts['--minimizer-spaces'] = s
        if protein:
            build_opts['--protein'] = None
        if max_db_size:
            build_opts['--max-db-size'] = max_db_size
        self.execute('kraken2-build', db, None, options=build_opts)

        # grab the taxdump for krona!
        if taxdump_out:
            shutil.copyfile(os.path.join(db, 'taxonomy', 'taxdump.tar.gz'), taxdump_out)

        # clean db
        self.execute('kraken2-build', db, None, options={
                '--threads':util.misc.sanitize_thread_count(num_threads),
                '--clean': None
                })

    def inspect(self, db, output, num_threads=None):
        with open(output, 'wt') as outf:
            subprocess.check_call('kraken2-inspect',
                '--db', db,
                '--threads', util.misc.sanitize_thread_count(num_threads),
                stdout=outf)

    def execute(self, command, db, output, args=None, options=None):
        '''Run a kraken-* command.

        Args:
          db: Kraken database directory.
          output: Output file of command.
          args: List of positional args.
          options: List of keyword options.
        '''
        options = options or {}

        if output:
            options['--output'] = output
        args = args or []

        cmd = [command, '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        subprocess.check_call(cmd)

    def pipeline(self, db, in_bams, out_reports=None, out_reads=None,
                 min_base_qual=None, confidence=None, num_threads=None):

        assert out_reads is not None or out_reports is not None
        out_reports = out_reports or []
        out_reads = out_reads or []

        for in_bam, out_read, out_report in itertools.zip_longest(in_bams, out_reads, out_reports):
            self.classify(in_bam, db, out_reads=out_read, out_report=out_report,
                min_base_qual=min_base_qual, confidence=confidence, num_threads=None)

    def classify(self, in_bam, db, out_reads=None, out_report=None,
                 confidence=None, min_base_qual=None, num_threads=None):
        """Classify input reads (bam)

        Args:
          in_bam: unaligned reads
          db: Kraken built database directory.
          out_reads: Output file of command.
        """

        if tools.samtools.SamtoolsTool().isEmpty(in_bam):
            # kraken cannot deal with empty input
            if out_reads:
                with open(out_reads, 'wt') as outf:
                    pass
            if out_report:
                with open(out_report, 'wt') as outf:
                    pass
            return

        opts = {
            '--threads': util.misc.sanitize_thread_count(num_threads)
        }
        if out_report:
            opts['--report'] = out_report
        if not out_reads:
            out_reads = '-' # in kraken2, this suppresses normal output
        if min_base_qual:
            opts['--minimum-base-quality'] = min_base_qual
        if confidence:
            opts['--confidence'] = confidence

        tmp_fastq1 = util.file.mkstempfname('.1.fastq.gz')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq.gz')
        tmp_fastq3 = util.file.mkstempfname('.s.fastq.gz')
        # Do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(in_bam, tmp_fastq1, tmp_fastq2, outFastq0=tmp_fastq3,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        if out_report:
            opts['--report'] = out_report
        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
            log.warn("running in single-end read mode!")
            res = self.execute('kraken2', db, out_reads, args=[tmp_fastq3], options=opts)
        else:
            opts['--paired'] = None
            res = self.execute('kraken2', db, out_reads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)
        os.unlink(tmp_fastq3)

