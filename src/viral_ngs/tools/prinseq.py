"tools.Tool for prinseq."

import logging
import os.path
import shutil
import subprocess
import tools
import util.file

TOOL_NAME = "prinseq-lite.pl"

log = logging.getLogger(__name__)

class PrinseqTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(PrinseqTool, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '-version']).decode('UTF-8').strip().split()[1]

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

    def rmdup_fastq_paired(self, inFastq1, inFastq2, outFastq1, outFastq2, includeUnmated=False, unpairedOutFastq1=None, unpairedOutFastq2=None):
        ''' remove duplicate reads and reads with multiple Ns '''

        if not os.path.exists(inFastq2) or os.path.getsize(inFastq2) == 0:
            # input is single-end
            log.debug("prinseq input is unpaired")
            util.file.touch(outFastq2)
            if not includeUnmated:
                # excluding unmated on unpaired input implies no output
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
            
            # if the user desires the unmated reads to be included in the output fastqs
            if includeUnmated:
                util.file.cat(outFastq1, (outprefix+'_1.fastq', outprefix+'_1_singletons.fastq'))
                util.file.cat(outFastq2, (outprefix+'_2.fastq', outprefix+'_2_singletons.fastq'))
            else:
                shutil.copyfile(outprefix+'_1.fastq', outFastq1)
                shutil.copyfile(outprefix+'_2.fastq', outFastq2)

            # if a path is specified for the unmated reads
            if unpairedOutFastq1 is not None:
                shutil.copyfile(outprefix+'_1_singletons.fastq', unpairedOutFastq1)

            if unpairedOutFastq2 is not None:
                shutil.copyfile(outprefix+'_2_singletons.fastq', unpairedOutFastq2)

            os.unlink(outprefix+'_1.fastq')
            os.unlink(outprefix+'_2.fastq')
            os.unlink(outprefix+'_1_singletons.fastq')
            os.unlink(outprefix+'_2_singletons.fastq')
