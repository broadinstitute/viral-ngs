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
import resource
import subprocess
import tempfile
import shutil
import sys
import tools

TOOL_NAME = "trinity"
TOOL_VERSION = "2011-11-26"
CONDA_TOOL_VERSION = "date.2011_11_26" # conda versions cannot have hyphens...

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
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable="Trinity", version=CONDA_TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self,
                inFastq1,
                inFastq2,
                outFasta,
                min_contig_length=300,
                JVMmemory=None,
                threads=1):    # pylint: disable=W0221
        if JVMmemory is None:
            JVMmemory = self.jvm_mem_default
        outdir = tempfile.mkdtemp(prefix='trinity-')
        if int(threads) < 1:
            threads = 1
        cmd = [
            self.install_and_get_path(), '--CPU', '{}'.format(int(threads)), '--bflyHeapSpace', JVMmemory.upper(),
            '--min_contig_length', str(min_contig_length), '--seqType', 'fq', '--left', inFastq1, '--right', inFastq2,
            '--output', outdir
        ]
        log.debug(' '.join(cmd))
        with unlimited_stack():
            subprocess.check_call(cmd)
        shutil.copyfile(os.path.join(outdir, 'Trinity.fasta'), outFasta)
        shutil.rmtree(outdir, ignore_errors=True)
