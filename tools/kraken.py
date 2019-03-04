'''
KRAKEN metagenomics classifier
'''
from __future__ import print_function
import collections
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

KRAKEN_VERSION = '1.0.0_fork3'
KRAKENUNIQ_VERSION = '0.5.7_yesimon'


log = logging.getLogger(__name__)


class Kraken(tools.Tool):

    BINS = {
        'classify': 'kraken',
        'build': 'kraken-build',
        'filter': 'kraken-filter',
        'report': 'kraken-report'}

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "kraken"
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage('kraken', executable=self.subtool_name, version=KRAKEN_VERSION, channel='broad-viral'))
        super(Kraken, self).__init__(install_methods=install_methods)

    def version(self):
        return KRAKEN_VERSION

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
        self.execute(self.BINS['build'], db, db, options=options,
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
            db_config = yaml.load(obj.get()['Body'].read(), loader=yaml.FullLoader)

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

            kraken_bin = 'kraken'
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
                        cmd += ' | tee >(pigz --best > {kraken_reads})'
                    else:
                        cmd += ' | pigz --best > {kraken_reads}'

                    cmd = cmd.format(kraken_reads=outReads[i])

                if outReports:
                    if filterThreshold is not None:

                        kraken_filter_bin = 'kraken-filter'
                        cmd += ' | {kraken_filter}{tax_opts} --threshold {filterThreshold}'.format(
                            kraken_filter=kraken_filter_bin,
                            tax_opts=tax_filter_opts,
                            filterThreshold=filterThreshold)

                    kraken_report_bin = 'kraken-report'
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
                bam2fq_ps = picard.execute(in_bam, fastq_pipes[i*2], fastq_pipes[i*2 + 1],
                    picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                    JVMmemory=picard.jvmMemDefault, background=True)

                log.debug('Calling kraken output command line: %s', cmd)
                subprocess.check_call(cmd, shell=True, executable='/bin/bash', env=env)

                if bam2fq_ps.poll():
                    raise subprocess.CalledProcessError(bam2fq_ps.returncode, "SamToFastqTool().execute({})".format(in_bam))


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
        # Detect if input bam was paired by checking fastq 2
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
        self.execute(self.BINS['filter'], db, outReads, args=[inReads],
                            options={'--threshold': filterThreshold})

    def report(self, inReads, db, outReport):
        """Convert Kraken read-based output to summary reports
        """
        self.execute(self.BINS['report'], db, outReport, args=[inReads])

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
        options = options or {}

        if command == self.BINS['classify']:
            if output:
                options['--output'] = output
            elif 'krakenuniq' in command:
                options['--output'] = 'off'
        option_string = option_string or ''
        args = args or []

        cmd = [command, '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(shlex.split(option_string))
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        if command == self.BINS['classify']:
            subprocess.check_call(cmd)
        elif command == self.BINS['build']:
            subprocess.check_call(cmd)
        else:
            with util.file.open_or_gzopen(output, 'w') as of:
                util.misc.run(cmd, stdout=of, stderr=subprocess.PIPE, check=True)


@tools.skip_install_test()
class Jellyfish(Kraken):
    """ Tool wrapper for Jellyfish (installed by kraken-all metapackage) """
    subtool_name = 'jellyfish'


class KrakenUniq(Kraken):

    BINS = {
        'classify': 'krakenuniq',
        'build': 'krakenuniq-build',
        'filter': 'krakenuniq-filter',
        'report': 'krakenuniq-report'}

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, 'subtool_name') else 'krakenuniq'
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage('krakenuniq', executable=self.subtool_name, version=KRAKENUNIQ_VERSION, channel='broad-viral'))
        super(KrakenUniq, self).__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def pipeline(self, db, in_bams, out_reports=None, out_reads=None,
                 filter_threshold=None, num_threads=None):

        try:
            from itertools import zip_longest
        except:  # Python 2 compat
            from itertools import izip_longest as zip_longest
        assert out_reads is not None or out_reports is not None
        out_reports = out_reports or []
        out_reads = out_reads or []

        for in_bam, out_read, out_report in zip_longest(in_bams, out_reads, out_reports):
            self.classify(in_bam, db, out_reads=out_read, out_report=out_report, num_threads=None)

    def classify(self, in_bam, db, out_reads=None, out_report=None, num_threads=None):
        """Classify input reads (bam)

        Args:
          in_bam: unaligned reads
          db: Kraken built database directory.
          outReads: Output file of command.
        """
        tmp_fastq1 = util.file.mkstempfname('.1.fastq.gz')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq.gz')
        # Do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(in_bam, tmp_fastq1, tmp_fastq2,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        opts = {
            '--threads': util.misc.sanitize_thread_count(num_threads),
            '--fastq-input': None,
            '--gzip-compressed': None,
            '--preload': None
        }
        if out_report:
            opts['--report-file'] = out_report
        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < 50:
            res = self.execute(self.BINS['classify'], db, out_reads, args=[tmp_fastq1], options=opts)
        else:
            opts['--paired'] = None
            res = self.execute(self.BINS['classify'], db, out_reads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)
        if out_report:
            with open(out_report, 'rt+') as f:
                lines = [line.strip() for line in f.readlines() if not line.startswith('#')]
                lines = [line for line in lines if line]
                if not lines:
                    f.seek(f.tell() - 1, os.SEEK_SET)
                    print('\t'.join(['%', 'reads', 'taxReads', 'kmers', 'dup', 'cov', 'taxID', 'rank', 'taxName']), file=f)
                    print('\t'.join(['100.00', '0', '0', '0', '0', 'NA', '0', 'no rank', 'unclassified']), file=f)

    def read_report(self, report_fn):
        report = collections.Counter()
        with open(report_fn) as f:
            for line in f:
                if line.startswith('#') or line.startswith('%'):
                    continue
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                percent = float(parts[0])
                cum_reads = int(parts[1])
                tax_reads = int(parts[2])
                tax_kmers = int(parts[3])
                if parts[5] == 'NA':  # unclassified
                    cov = 0
                else:
                    cov = float(parts[5])
                tax_id = int(parts[6])
                rank = parts[7]
                name = parts[8]
                report[tax_id] = (tax_reads, tax_kmers)
        return report



class Kraken2(tools.Tool):

    BINS = {
        'classify': 'kraken2',
        'build': 'kraken2-build'
    }

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "kraken2"
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage('kraken2', executable=self.subtool_name, version=KRAKEN_VERSION, channel='broad-viral'))
        super(Kraken2, self).__init__(install_methods=install_methods)

    def version(self):
        return KRAKEN2_VERSION

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
        self.execute(self.BINS['build'], db, db, options=options,
                     option_string=option_string)

    def pipeline(self, db, inBams, outReports=None, outReads=None,
                 lockMemory=None, filterThreshold=None, numThreads=None):
        assert outReads is not None or outReports is not None


        try:
            from itertools import zip_longest
        except:  # Python 2 compat
            from itertools import izip_longest as zip_longest
        assert out_reads is not None or out_reports is not None
        out_reports = out_reports or []
        out_reads = out_reads or []

        for in_bam, out_read, out_report in zip_longest(in_bams, out_reads, out_reports):
            self.classify(in_bam, db, out_reads=out_read, out_report=out_report, num_threads=None)


    def classify(self, in_bam, db, out_reads=None, out_report=None, num_threads=None):
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
        picard.execute(in_bam, tmp_fastq1, tmp_fastq2,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        opts = {
            '--threads': util.misc.sanitize_thread_count(num_threads),
            '--fastq-input': None,
            '--gzip-compressed': None,
        }
        if out_report:
            opts['--report'] = out_report
        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < 50:
            res = self.execute('kraken2', db, out_reads, args=[tmp_fastq1], options=opts)
        else:
            opts['--paired'] = None
            res = self.execute('kraken2', db, out_reads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)

    def filter(self, inReads, db, outReads, filterThreshold):
        """Filter Kraken hits
        """
        self.execute(self.BINS['filter'], db, outReads, args=[inReads],
                            options={'--threshold': filterThreshold})

    def report(self, in_reads, db, outReport):
        """Convert Kraken read-based output to summary reports
        """
        self.execute(self.BINS['report'], db, outReport, args=[inReads])

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
        options = options or {}

        if command == self.BINS['classify']:
            if output:
                options['--output'] = output
        option_string = option_string or ''
        args = args or []

        cmd = [command, '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(shlex.split(option_string))
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        if command == self.BINS['classify']:
            subprocess.check_call(cmd)
        elif command == self.BINS['build']:
            subprocess.check_call(cmd)
        else:
            with util.file.open_or_gzopen(output, 'w') as of:
                util.misc.run(cmd, stdout=of, stderr=subprocess.PIPE, check=True)
