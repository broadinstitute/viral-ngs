"tools.Tool for prinseq."

import logging
import os.path
import shutil
import subprocess
import tools

TOOL_NAME = "prinseq"
TOOL_VERSION = '0.20.4'

log = logging.getLogger(__name__)

class PrinseqTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable="prinseq-lite.pl", version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)

    def rmdup_fastq_single(inFastq, outFastq):
	    # remove duplicate reads and reads with multiple Ns
	    if not outFastq.endswith('.fastq'):
	        raise Exception()
	    if os.path.getsize(inFastq) == 0:
	        # prinseq-lite fails on empty file input (which can happen in real life
	        # if no reads match the refDbs) so handle this scenario specially
	        log.info("output is empty: no reads in input match refDb")
	        shutil.copyfile(inFastq, outFastq)
	    else:
	        cmd = [
	            'perl', self.install_and_get_path(), '-ns_max_n', '1', '-derep', '1', '-fastq',
	            inFastq, '-out_bad', 'null', '-line_width', '0', '-out_good', outFastq[:-6]
	        ]
	        log.debug(' '.join(cmd))
	        subprocess.check_call(cmd)

    def rmdup_fastq_paired(inFastq1, inFastq2, outFastq1, outFastq2):
        self.rmdup_fastq_single(inFastq1, outFastq1)
        self.rmdup_fastq_single(inFastq2, outFastq2)
        return
	    # remove duplicate reads and reads with multiple Ns
	    if not outFastq.endswith('.fastq'):
	        raise Exception()
	    if os.path.getsize(inFastq) == 0:
	        # prinseq-lite fails on empty file input (which can happen in real life
	        # if no reads match the refDbs) so handle this scenario specially
	        log.info("output is empty: no reads in input match refDb")
	        shutil.copyfile(inFastq, outFastq)
	    else:
	        cmd = [
	            'perl', self.install_and_get_path(), '-ns_max_n', '1', '-derep', '1', '-fastq',
	            inFastq, '-out_bad', 'null', '-line_width', '0', '-out_good', outFastq[:-6]
	        ]
	        log.debug(' '.join(cmd))
	        subprocess.check_call(cmd)
