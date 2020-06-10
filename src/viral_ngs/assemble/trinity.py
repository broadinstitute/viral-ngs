'''
    The Trinity RNA-SEQ assembler

    This uses an older version of Trinity that uses an older
    assembly algorithm that works better with highly diverse
    viral genomes.
'''

import contextlib
import logging
import os
import os.path
import platform
import resource
import subprocess
import tempfile
import shutil
import sys
import tools
import util.misc

TOOL_NAME = "Trinity"

log = logging.getLogger(__name__)


@contextlib.contextmanager
def unlimited_stack():
    '''Set the ulimit on stack size to be infinity.

    OS X has a fixed hard stack size limit of 64 MB, so we're not setting it here.
    '''
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    if sys.platform.startswith('linux'):
        new_soft, new_hard = (resource.RLIM_INFINITY, resource.RLIM_INFINITY)
        try:
            resource.setrlimit(resource.RLIMIT_STACK, (new_soft, new_hard))
            yield
            resource.setrlimit(resource.RLIMIT_STACK, (soft, hard))
        except (ValueError, OSError) as e:
            log.warning('Error raising stacksize to unlimited: %s', str(e))
            yield
    else:
        yield


class TrinityTool(tools.Tool):
    jvm_mem_default = '4g'

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(TrinityTool, self).__init__(install_methods=install_methods)

    def version(self):
        return None

    def execute(self,
                inFastq1,
                inFastq2,
                outFasta,
                min_contig_length=300,
                JVMmemory=None,
                threads=None):    # pylint: disable=W0221
        if JVMmemory is None:
            JVMmemory = self.jvm_mem_default
        outdir = tempfile.mkdtemp(prefix='trinity-')
        util.misc.sanitize_thread_count(threads)
        cmd = [
            self.install_and_get_path(), '--CPU', '{}'.format(util.misc.sanitize_thread_count(threads)), '--bflyHeapSpace', JVMmemory.upper(),
            '--min_contig_length', str(min_contig_length), '--seqType', 'fq', '--left', inFastq1, '--right', inFastq2,
            '--output', outdir
        ]
        log.debug(' '.join(cmd))

        #
        # Fix some quirks of the Trinity.pl script
        #
        trinity_env = dict(os.environ)
        # Ensure OSTYPE is set
        if 'OSTYPE' not in trinity_env:
            trinity_env['OSTYPE'] = platform.system().lower()
        # Ensure _JAVA_OPTIONS is not set: OpenJDK java always prints a message to stdout that it picked up
        # _JAVA_OPTIONS, which confuses Trinity.pl's attempt to get the java version.
        # Note that the Java heap options will still be passed to trinity via the --bflyHeapSpace parameter.
        if '_JAVA_OPTIONS' in trinity_env:
            del trinity_env['_JAVA_OPTIONS']

        with unlimited_stack():
            subprocess.check_call(cmd, env=trinity_env)

        shutil.copyfile(os.path.join(outdir, 'Trinity.fasta'), outFasta)
        shutil.rmtree(outdir, ignore_errors=True)
