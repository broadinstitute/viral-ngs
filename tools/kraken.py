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

TOOL_NAME = 'kraken'
TOOL_VERSION = '1.0.0_fork4'

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
        options['--threads'] = util.misc.sanitize_thread_count(options.get('--threads'))
        self.execute('kraken-build', db, db, options=options,
                     option_string=option_string)

    def _db_opts(self, db, threads):
        '''Determine kraken command line options based on db string'''

        env = os.environ.copy()
        def s3_psub(path):
            cmd = 'aws s3 cp {} -'.format(path)
            if path.endswith('.bz2'):
                cmd += ' | lbzip2 -n {threads} -d'.format(threads=threads)
            elif path.endswith('.gz'):
                cmd += ' | pigz -p {threads} -d'.format(threads=threads)
            elif path.endswith('.lz4'):
                cmd += ' | lz4 -d'
            cmd = '<({cmd})'.format(cmd=cmd)
            return cmd

        if db.startswith('s3://'):
            import boto3
            import yaml
            s3 = boto3.resource('s3')
            path = db[5:]
            bucket_name, db_dir = path.split('/', 1)
            obj = s3.Object(bucket_name, '/'.join([db_dir, 'config.yaml']))
            db_config = yaml.load(obj.get()['Body'].read())

            db_opts = (' --db-pipe --db-index {index_f} --index-size {index_size} --db-file {db_f} --db-size {db_size} '
                       '--taxonomy {nodes}'.format(
                           index_f=s3_psub(db_config['db_index']),
                           index_size=db_config['index_size'],
                           db_f=s3_psub(db_config['db_file']),
                           db_size=db_config['db_size'],
                           nodes=s3_psub(db_config['taxonomy_nodes'])))
            tax_filter_opts = ' --taxonomy-nodes {nodes}'.format(
                nodes=s3_psub(db_config['taxonomy_nodes']))
            tax_report_opts = tax_filter_opts + ' --taxonomy-names {names}'.format(
                names=s3_psub(db_config['taxonomy_names']))
            env['KRAKEN_DEFAULT_DB'] = '.'
        else:
            env['KRAKEN_DEFAULT_DB'] = db
            db_opts = ''
            tax_filter_opts = ''
            tax_report_opts = ''
        return db_opts, env, tax_filter_opts, tax_report_opts

    def pipeline(self, db, inBams, outReports=None, outReads=None,
                 lockMemory=None, filterThreshold=None, numThreads=None):
        assert outReads is not None or outReports is not None

        n_bams = len(inBams)
        # 2n for paired fastq, 1n for kraken output
        n_pipes = n_bams * 3
        if outReports and len(outReports) != n_bams:
            raise Exception("--outReports specified with {} output files, which does not match the number of input bams ({})".format(len(outReports), n_bams))
        if outReads and len(outReads) != n_bams:
            raise Exception("--outReads specified with {} output files, which does not match the number of input bams ({})".format(len(outReads), n_bams))
        threads = util.misc.sanitize_thread_count(numThreads)

        with util.file.fifo(n_pipes) as pipes:
            fastq_pipes = pipes[:n_bams * 2]
            kraken_output_pipes = pipes[n_bams * 2:]

            kraken_bin = os.path.join(self.libexec, 'kraken')
            opts = ''
            if lockMemory:
                opts += ' --lock-memory'

            db_opts, env, tax_filter_opts, tax_report_opts = self._db_opts(db, threads)
            opts += db_opts

            cmd = '''set -ex -o pipefail; {kraken}{opts} --paired --fastq-input --threads {threads} {outputs} {fastqs}'''.format(
                kraken=kraken_bin,
                opts=opts,
                threads=threads,
                outputs=' '.join('--output {}'.format(x) for x in kraken_output_pipes),
                fastqs=' '.join(fastq_pipes))
            log.debug('Calling kraken command line: %s', cmd)
            subprocess.Popen(cmd, shell=True, executable='/bin/bash', env=env)

            for i, in_bam in enumerate(inBams):
                cmd = 'cat {kraken_output}'.format(kraken_output=kraken_output_pipes[i])

                if outReads:
                    if outReports:
                        cmd += ' | tee >(gzip --best > {kraken_reads})'
                    else:
                        cmd += ' | gzip --best > {kraken_reads}'

                    cmd = cmd.format(kraken_reads=outReads[i])

                if outReports:
                    if filterThreshold is not None:

                        kraken_filter_bin = os.path.join(self.libexec, 'kraken-filter')
                        cmd += ' | {kraken_filter}{tax_opts} --threshold {filterThreshold}'.format(
                            kraken_filter=kraken_filter_bin,
                            tax_opts=tax_filter_opts,
                            filterThreshold=filterThreshold)

                    kraken_report_bin = os.path.join(self.libexec, 'kraken-report')
                    cmd += ' | {kraken_report}{tax_opts} > {outReport}'.format(
                        kraken_report=kraken_report_bin,
                        tax_opts=tax_report_opts,
                        outReport=outReports[i])

                # do not convert this to samtools bam2fq unless we can figure out how to replicate
                # the clipping functionality of Picard SamToFastq
                picard = tools.picard.SamToFastqTool()
                picard_opts = {
                    'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
                    'CLIPPING_ACTION': 'X'
                }
                picard.execute(in_bam, fastq_pipes[i*2], fastq_pipes[i*2 + 1],
                               picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                               JVMmemory=picard.jvmMemDefault, background=True)

                log.debug('Calling kraken output command line: %s', cmd)
                subprocess.check_call(cmd, shell=True, executable='/bin/bash', env=env)


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

        opts = {
            '--threads': util.misc.sanitize_thread_count(numThreads),
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
