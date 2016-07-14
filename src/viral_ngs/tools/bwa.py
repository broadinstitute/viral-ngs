'''
    The BWA aligner.

'''

import logging
import os
import os.path
import subprocess
import shutil

import tools
import tools.samtools
import util.file
import util.misc

TOOL_NAME = 'bwa'
TOOL_VERSION = '0.7.15'

log = logging.getLogger(__name__)


class Bwa(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, args, stdout=None):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        subprocess.check_call(tool_cmd, stdout=stdout)
        if stdout:
            stdout.close()

    def index(self, inFasta, algorithm=None):
        cmd = []
        if algorithm is not None:
            if algorithm not in ('is', 'bwtsw'):
                raise NameError(algorithm + " is not a recognized algorithm")
            cmd.extend(('-a', algorithm))
        cmd.append(inFasta)
        self.execute('index', cmd)

    def mem(self, inReads, refDb, outAlign, opts=None, threads=None):
        opts = [] if not opts else opts

        threads = threads or util.misc.available_cpu_count()
        samtools = tools.samtools.SamtoolsTool()
        fq1 = util.file.mkstempfname('.1.fastq')
        fq2 = util.file.mkstempfname('.2.fastq')
        aln_sam = util.file.mkstempfname('.sam')
        aln_sam_sorted = util.file.mkstempfname('sorted.sam')
        samtools.bam2fq(inReads, fq1, fq2)
        self.execute('mem', opts + ['-t', str(threads), refDb, fq1, fq2], stdout=aln_sam)
        os.unlink(fq1)
        os.unlink(fq2)
        samtools.sort(aln_sam, aln_sam_sorted)
        os.unlink(aln_sam)
        # cannot index sam files; only do so if a bam is desired
        if outAlign.endswith(".bam") or outAlign.endswith(".cram"):
            # convert sam -> bam
            samtools.view(["-b"], aln_sam_sorted, outAlign)
            samtools.index(outAlign)
        elif outAlign.endswith(".sam"):
            shutil.copyfile(aln_sam_sorted, outAlign)
        os.unlink(aln_sam_sorted)


