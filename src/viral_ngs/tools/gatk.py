'''
    GATK genotyping toolkit from the Broad Institute

    This software has different licenses depending on use cases.
    As such, we do not have an auto-downloader. The user must have GATK
    pre-installed on their own and available in $GATK_PATH.
'''

import tools
import tools.picard
import tools.samtools
import util.file
import util.misc

import logging
import os
import os.path
import shutil
import subprocess
import tempfile

_log = logging.getLogger(__name__)

TOOL_NAME = 'gatk3'

class GATKTool(tools.Tool):
    jvmMemDefault = '2g'

    def __init__(self, path=None):
        self.tool_version = None
        install_methods = []
        for jarpath in [path, os.environ.get('GATK_PATH')]:
            if jarpath:
                if not jarpath.endswith('.jar'):
                    jarpath = os.path.join(jarpath, 'GenomeAnalysisTK.jar')
                install_methods.append(
                    tools.PrexistingUnixCommand(
                        jarpath,
                        verifycmd='java -jar %s --version &> /dev/null' % jarpath,
                        verifycode=0,
                        require_executability=False
                    )
                )
        install_methods.append(tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True))
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, command, gatkOptions=None, JVMmemory=None):    # pylint: disable=W0221
        gatkOptions = gatkOptions or []

        if not JVMmemory:
            JVMmemory = self.jvmMemDefault

        # the conda version wraps the jar file with a shell script
        if self.install_and_get_path().endswith(".jar"):
            tool_cmd = [
                'java', '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), '-jar', self.install_and_get_path(),
                '-T', command
            ] + list(map(str, gatkOptions))
        else:
            tool_cmd = [
                self.install_and_get_path(), '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), '-T', command
            ] + list(map(str, gatkOptions))

        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    @staticmethod
    def dict_to_gatk_opts(options):
        return ["%s=%s" % (k, v) for k, v in options.items()]

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        return subprocess.check_output([self.install_and_get_path(), '--version']).decode('UTF-8').strip()

    def ug(self, inBam, refFasta, outVcf, options=None, JVMmemory=None, threads=None):
        options = options or ["--min_base_quality_score", 15, "-ploidy", 4]

        threads = util.misc.sanitize_thread_count(threads)
        opts = [
            '-I', inBam,
            '-R', refFasta,
            '-o', outVcf,
            '-glm', 'BOTH',
            '--baq', 'OFF',
            '--useOriginalQualities',
            '-out_mode', 'EMIT_ALL_SITES',
            '-dt', 'NONE',
            '--num_threads', threads,
            '-A', 'AlleleBalance',
        ]
        if tuple(map(int, self.version().split('-')[0].split('.'))) < (3, 7):
            opts += ['-stand_call_conf', 0, '-stand_emit_conf', 0]  # deprecated in 3.7+

        self.execute('UnifiedGenotyper', opts + options, JVMmemory=JVMmemory)

    def local_realign(self, inBam, refFasta, outBam, JVMmemory=None, threads=None):
        threads = util.misc.sanitize_thread_count(threads)
        _log.debug("Running local realign with %s threads", threads)
        intervals = util.file.mkstempfname('.intervals')
        opts = ['-I', inBam, '-R', refFasta, '-o', intervals, '--num_threads', threads]
        self.execute('RealignerTargetCreator', opts, JVMmemory=JVMmemory)
        opts = ['-I', inBam, '-R', refFasta, '-targetIntervals', intervals, '-o', outBam]
        self.execute('IndelRealigner', opts, JVMmemory=JVMmemory)
        os.unlink(intervals)
