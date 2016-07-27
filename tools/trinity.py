'''
    The Trinity RNA-SEQ assembler

    This uses an older version of Trinity that uses an older
    assembly algorithm that works better with highly diverse
    viral genomes.
'''

import logging
import os
import os.path
import subprocess
import tempfile
import shutil
import tools

TOOL_NAME = "trinity"
TOOL_VERSION = "2011-11-26"
CONDA_TOOL_VERSION = "date.2011_11_26" # conda versions cannot have hyphens...

log = logging.getLogger(__name__)


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
        subprocess.check_call(cmd)
        shutil.copyfile(os.path.join(outdir, 'Trinity.fasta'), outFasta)
        shutil.rmtree(outdir, ignore_errors=True)

