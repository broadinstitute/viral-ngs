'''
KRAKEN metagenomics classifier
'''
from __future__ import print_function
import itertools
import logging
import os
import os.path
import shlex
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

TOOL_NAME = 'kraken-all'
TOOL_VERSION = '0.10.6_fork2'

log = logging.getLogger(__name__)

class Kraken(tools.Tool):

    BINS = ['kraken', 'kraken-build', 'kraken-filter', 'kraken-mpa-report', 'kraken-report', 'kraken-translate']

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "kraken"
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION, channel='broad-viral'))
        super(Kraken, self).__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())

    def build(self, db, options=None, option_string=None):
        '''Create a kraken database.

        Args:
          db: Kraken database directory to build. Must have library/ and
            taxonomy/ subdirectories to build from.
          *args: List of input filenames to process.
        '''
        self.execute('kraken-build', db, db, options=options,
                     option_string=option_string)

    def pipeline(self, inBam, db, outReport=None, outReads=None, filterThreshold=None, numThreads=None):
        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # kraken cannot deal with empty input
            with open(outReads, 'rt') as outf:
                pass
            return

        with util.file.fifo(2) as (fastq1_pipe, fastq2_pipe):
            # do not convert this to samtools bam2fq unless we can figure out how to replicate
            # the clipping functionality of Picard SamToFastq
            picard = tools.picard.SamToFastqTool()
            picard_opts = {
                'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
                'CLIPPING_ACTION': 'X'
            }
            picard.execute(inBam, fastq1_pipe, fastq2_pipe,
                           picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                           JVMmemory=picard.jvmMemDefault,
                           background=True)

            if numThreads is None:
                numThreads = 10000000
            opts = {
                '--threads': min(int(numThreads), util.misc.available_cpu_count()),
            }

            kraken_bin = os.path.join(self.libexec, 'kraken')
            cmd = '''export KRAKEN_DEFAULT_DB={kraken_db}; {{ {kraken} --paired --fastq-input --threads {threads} {fastq1} {fastq2} 2>&1 1>&3 3>&- | sed '/Processed [0-9]* sequences/d'; }} \
            3>&1 1>&2'''.format(
                kraken_db=db,
                kraken=kraken_bin,
                threads=numThreads,
                fastq1=fastq1_pipe,
                fastq2=fastq2_pipe)

            if outReads is not None:
                cmd += '| tee >(gzip > {kraken_reads})'.format(kraken_reads=outReads)

            if filterThreshold is not None:

                kraken_filter_bin = os.path.join(self.libexec, 'kraken-filter')
                cmd += '| {kraken_filter} --threshold {filterThreshold}'.format(
                    kraken_filter=kraken_filter_bin,
                    filterThreshold=filterThreshold)

            if outReport is not None:
                kraken_report_bin = os.path.join(self.libexec, 'kraken-report')
                cmd += '| {kraken_report} > {outReport}'.format(
                    kraken_report=kraken_report_bin,
                    outReport=outReport)
            subprocess.check_call(cmd, shell=True, executable='/bin/bash')


    def classify(self, inBam, db, outReads, numThreads=None):
        """Classify input reads (bam)

        Args:
          inBam: unaligned reads
          db: Kraken built database directory.
          outReads: Output file of command.
        """
        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # kraken cannot deal with empty input
            with open(outReads, 'rt') as outf:
                pass
            return
        tmp_fastq1 = util.file.mkstempfname('.1.fastq.gz')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq.gz')
        # do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(inBam, tmp_fastq1, tmp_fastq2,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        if numThreads is None:
            numThreads = 10000000
        opts = {
            '--threads': min(int(numThreads), util.misc.available_cpu_count()),
            '--fastq-input': None,
            '--gzip-compressed': None,
        }
        if os.path.getsize(tmp_fastq2) < 50:
            res = self.execute('kraken', db, outReads, args=[tmp_fastq1], options=opts)
        else:
            opts['--paired'] = None
            res = self.execute('kraken', db, outReads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)

    def filter(self, inReads, db, outReads, filterThreshold):
        """Filter Kraken hits
        """
        self.execute('kraken-filter', db, outReads, args=[inReads],
                            options={'--threshold': filterThreshold})

    def report(self, inReads, db, outReport):
        """Convert Kraken read-based output to summary reports
        """
        self.execute('kraken-report', db, outReport, args=[inReads])

    def execute(self, command, db, output, args=None, options=None,
                option_string=None):
        '''Run a kraken-* command.

        Args:
          db: Kraken database directory.
          output: Output file of command.
          args: List of positional args.
          options: List of keyword options.
          option_string: Raw strip command line options.
        '''
        assert command in Kraken.BINS, 'Kraken command is unknown'
        options = options or {}

        if command == 'kraken':
            options['--output'] = output
        option_string = option_string or ''
        args = args or []

        cmd = [os.path.join(self.libexec, command), '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(shlex.split(option_string))
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        if command == 'kraken':
            subprocess.check_call(cmd)
        elif command == 'kraken-build':
            jellyfish_path = Jellyfish().install_and_get_path()
            env = os.environ.copy()
            env['PATH'] = ':'.join([os.path.dirname(jellyfish_path), env['PATH']])
            subprocess.check_call(cmd, env=env)
        else:
            with util.file.open_or_gzopen(output, 'w') as of:
                util.misc.run(cmd, stdout=of, stderr=subprocess.PIPE, check=True)


class Jellyfish(Kraken):
    """ Tool wrapper for Jellyfish (installed by kraken-all metapackage) """
    subtool_name = 'jellyfish'
