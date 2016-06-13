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
import tools
import tools.picard
import util.file
import util.misc
from builtins import super

TOOL_NAME = "kraken-all"
TOOL_VERSION = '0.10.6_eaf8fb68'

log = logging.getLogger(__name__)

@tools.skip_install_test(condition=tools.is_osx())
class Kraken(tools.Tool):

    BINS = ['kraken', 'kraken-build', 'kraken-filter', 'kraken-mpa-report', 'kraken-report', 'kraken-translate']

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "kraken"
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
        super().__init__(install_methods=install_methods)

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
        return self.execute('kraken-build', db, db, options=options,
                            option_string=option_string)
    
    def classify(self, inBam, db, outReads, numThreads=None):
        """Classify input reads (bam)

        Args:
          inBam: unaligned reads
          db: Kraken built database directory.
          outReads: Output file of command.
        """
        # add hardening code here to defend against div-by-zero bugs
        tmp_fastq1 = util.file.mkstempfname('.1.fastq')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq')
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
            '--paired': None,
            '--threads': min(int(numThreads), util.misc.available_cpu_count()),
        }
        res = self.execute('kraken', db, outReads, args=[tmp_fastq1, tmp_fastq2], options=opts)
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)
        return res
        
    def filter(self, inReads, db, outReads, filterThreshold):
        """Filter Kraken hits
        """
        # add hardening code here to defend against div-by-zero bugs
        self.execute('kraken-filter', db, outReads, args=[inReads],
                            options={'--threshold': filterThreshold})
    
    def report(self, inReads, db, outReport):
        """Convert Kraken read-based output to summary reports
        """
        # add hardening code here to defend against div-by-zero bugs
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

        jellyfish_path = Jellyfish().install_and_get_path()
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(os.path.dirname(jellyfish_path),
                                     env['PATH'])
        cmd = [os.path.join(self.libexec, command), '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(shlex.split(option_string))
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        if command in ('kraken', 'kraken-build'):
            return util.misc.run_and_print(cmd, env=env, check=True)
        else:
            with util.file.open_or_gzopen(output, 'w') as of:
                res = util.misc.run(cmd, stdout=of, stderr=subprocess.PIPE,
                                    check=True)
                if res.returncode != 0:
                    print(res.stderr.decode('utf-8'), file=sys.stderr)
            return res

@tools.skip_install_test(condition=tools.is_osx())
class Jellyfish(Kraken):
    """ Tool wrapper for Jellyfish (installed by kraken-all metapackage) """
    subtool_name = 'jellyfish'
