"tools.Tool for prinseq."

import logging
import os.path
import shutil
import subprocess
import tools
import util.file

TOOL_NAME = "prinseq"
TOOL_VERSION = '0.20.4'

log = logging.getLogger(__name__)

class PrinseqTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable="prinseq-lite.pl", version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)

    def rmdup_fastq_single(self, inFastq, outFastq):
        ''' remove duplicate reads and reads with multiple Ns '''

        if os.path.getsize(inFastq) == 0:
            # prinseq-lite fails on empty file input so handle this scenario specially
            log.debug("prinseq input is empty")
            shutil.copyfile(inFastq, outFastq)
        else:
            outprefix = util.file.mkstempfname("-prinseq-rmdup-out")
            cmd = [
                'perl', self.install_and_get_path(), '-verbose',
                '-ns_max_n', '1', '-derep', '1',
                '-fastq', inFastq,
                '-out_bad', 'null', '-line_width', '0', '-out_good', outprefix
            ]
            log.debug(' '.join(cmd))
            subprocess.check_call(cmd)
            shutil.copyfile(outprefix+'.fastq', outFastq)
            os.unlink(outprefix+'.fastq')

    def rmdup_fastq_paired(self, inFastq1, inFastq2, outFastq1, outFastq2, purgeUnmated=False):
        ''' remove duplicate reads and reads with multiple Ns '''

        if not os.path.exists(inFastq2) or os.path.getsize(inFastq2) == 0:
            # input is single-end
            log.debug("prinseq input is unpaired")
            util.file.touch(outFastq2)
            if purgeUnmated:
                # purgeUnmated on unpaired input implies no output
                util.file.touch(outFastq1)
            else:
                self.rmdup_fastq_single(inFastq1, outFastq1)
        elif os.path.getsize(inFastq1) == 0:
            # prinseq-lite fails on empty file input so handle this scenario specially
            log.debug("prinseq input is empty")
            util.file.touch(outFastq1)
            util.file.touch(outFastq2)
        else:
            outprefix = util.file.mkstempfname("-prinseq-rmdup-out")
            cmd = [
                'perl', self.install_and_get_path(), '-verbose',
                '-ns_max_n', '1', '-derep', '1',
                '-fastq', inFastq1, '-fastq2', inFastq2,
                '-out_bad', 'null', '-line_width', '0', '-out_good', outprefix
            ]
            log.debug(' '.join(cmd))
            subprocess.check_call(cmd)
            for fn in (outprefix+'_1.fastq', outprefix+'_1_singletons.fastq',
                    outprefix+'_2.fastq', outprefix+'_2_singletons.fastq'):
                # if any of the four output files is empty, prinseq doesn't create it at all, which is bad for us
                util.file.touch(fn)
            if purgeUnmated:
                shutil.copyfile(outprefix+'_1.fastq', outFastq1)
                shutil.copyfile(outprefix+'_2.fastq', outFastq2)
            else:
                util.file.cat(outFastq1, (outprefix+'_1.fastq', outprefix+'_1_singletons.fastq'))
                util.file.cat(outFastq2, (outprefix+'_2.fastq', outprefix+'_2_singletons.fastq'))
            os.unlink(outprefix+'_1.fastq')
            os.unlink(outprefix+'_2.fastq')
            os.unlink(outprefix+'_1_singletons.fastq')
            os.unlink(outprefix+'_2_singletons.fastq')
