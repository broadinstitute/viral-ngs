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
CONDA_TOOL_VERSION = "2.1.1"
TRINITY_VERSION = "trinityrnaseq_r{}".format(TOOL_VERSION)
url = "http://sourceforge.net/projects/trinityrnaseq/files/{}.tgz".format(TRINITY_VERSION)

log = logging.getLogger(__name__)


class TrinityTool(tools.Tool):
    jvm_mem_default = '4g'

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append( tools.CondaPackage(TOOL_NAME, version=CONDA_TOOL_VERSION) )
            install_methods.append(DownloadAndBuildTrinity(url, TRINITY_VERSION + '/Trinity.pl'))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, inFastq1, inFastq2, outFasta, min_contig_length=300, JVMmemory=None, threads=1): # pylint: disable=W0221
        if JVMmemory is None:
            JVMmemory = self.jvm_mem_default
        outdir = tempfile.mkdtemp(prefix='trinity-')
        if int(threads) < 1:
            threads = 1
        cmd = [self.install_and_get_path(), '--CPU', '{}'.format(int(threads)), '--bflyHeapSpace', JVMmemory.upper(),
               '--min_contig_length', str(min_contig_length), '--seqType', 'fq', '--left', inFastq1, '--right',
               inFastq2, '--output', outdir]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)
        shutil.copyfile(os.path.join(outdir, 'Trinity.fasta'), outFasta)
        shutil.rmtree(outdir, ignore_errors=True)


class DownloadAndBuildTrinity(tools.DownloadPackage):

    def post_download(self):
        trinity_dir = os.path.join(self.destination_dir, TRINITY_VERSION)
        if TOOL_VERSION == "2011-11-26":
            # Chrysalis doesn't compile. Need to add an include file.
            badFilePath = os.path.join(trinity_dir, 'Chrysalis', 'analysis', 'RunButterfly.cc')
            os.rename(badFilePath, badFilePath + '.orig')
            with open(badFilePath, 'wt') as outf:
                outf.write('#include <unistd.h>\n')
                with open(badFilePath + '.orig', 'rt') as inf:
                    for line in inf:
                        outf.write(line)

            # Trinity.pl insists on Java 1.6, but Java >= 1.6 is fine
            badFilePath = os.path.join(trinity_dir, 'Trinity.pl')
            os.rename(badFilePath, badFilePath + '.orig')
            with open(badFilePath, 'wt') as outf:
                with open(badFilePath + '.orig', 'rt') as inf:
                    for line in inf:
                        if line.startswith('unless ($java_version =~ /java version'):
                            outf.write(r'$java_version =~ /java version "1\.(\d+)\./;')
                            outf.write('\nunless ($1 >= 6) {\n')
                        else:
                            outf.write(line)
            shutil.copymode(badFilePath + '.orig', badFilePath)

        # Now we can make:
        os.system('cd "{}" && make -s'.format(trinity_dir))
        shutil.rmtree(os.path.join(trinity_dir, 'sample_data'), ignore_errors=True)

    def verify_install(self):
        if not tools.DownloadPackage.verify_install(self):
            return False
        # Verify that chrysalis and inchworm were built
        trinity_dir = os.path.join(self.destination_dir, TRINITY_VERSION)
        chrysalisPath = os.path.join(trinity_dir, 'Chrysalis', 'Chrysalis')
        inchwormPath = os.path.join(trinity_dir, 'Inchworm', 'src', 'inchworm')
        for path in [chrysalisPath, inchwormPath]:
            if not os.access(path, (os.X_OK | os.R_OK)):
                log.debug('%s was not built.', path)
                self.installed = False
        self.installed = True
        return self.installed
