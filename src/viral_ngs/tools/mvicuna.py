"tools.Tool for mvicuna."

import logging
import os
import subprocess
import shutil

import tools
import util.file

# BroadUnixPath = '/gsap/garage-viral/viral/analysis/xyang/programs'\
#                 '/M-Vicuna/bin/mvicuna'

TOOL_NAME = "mvicuna"
TOOL_VERSION = "1.0"

_log = logging.getLogger(__name__)


class MvicunaTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            path = _get_mvicuna_path()
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
            install_methods.append(tools.PrexistingUnixCommand(path))
        tools.Tool.__init__(self, install_methods=install_methods)

    def rmdup(self, inPair, outPair, outUnpaired=None):
        """
        Run mvicuna's duplicate removal operation on paired-end input reads in
            fastq format, producing various outputs in fastq format.
        inPair, pairedOutPair are pairs of file names,
            while unpairedOut is a single file name.
        Notes on weird behaviors of M-Vicuna DupRm:
            For some reason, it requires you to specify both -opfq and -drm_op.
            The -drm_op pair (here, tmp1OutPair) is where the output is initially written.
            Then M-Vicuna renames the -drm_op files to the -opfq filenames (here, tmp2OutPair).
            So obviously, we use throwaway (tempfile) names for -drm_op and use the real
            desired output file names in -opfq.
            The problem is that M-Vicuna uses a rename/move operating system function, which
            means your tempfiles cannot be on a different file system than the final output
            files. This is our typical use case (local disks for tempfile.tempdir and
            network file systems for final output). Hence, our wrapper sets -opfq to yet
            another set of temp file names, and we move the final output files ourselves
            using shutil.move (which is capable of crossing file system boundaries).
        """
        if not outUnpaired:
            outUnpaired = util.file.mkstempfname(suffix='.unpaired.fastq')
        tmp1OutPair = (
            util.file.mkstempfname(suffix='.tmp1out.1.fastq'), util.file.mkstempfname(suffix='.tmp1out.2.fastq')
        )
        tmp2OutPair = (
            util.file.mkstempfname(suffix='.tmp2out.1.fastq'), util.file.mkstempfname(suffix='.tmp2out.2.fastq')
        )
        cmdline = [
            self.install_and_get_path(), '-ipfq', ','.join(inPair), '-opfq', ','.join(tmp2OutPair), '-osfq',
            outUnpaired, '-drm_op', ','.join(tmp1OutPair), '-tasks', 'DupRm'
        ]
        _log.debug(' '.join(cmdline))
        subprocess.check_call(cmdline)
        for tmpfname, outfname in zip(tmp2OutPair, outPair):
            shutil.copyfile(tmpfname, outfname)

    def rmdup_single(self, inFastq, outFastq):
        """
        Run mvicuna's duplicate removal operation on single-end input reads in
            fastq format, producing output in fastq format.
        Notes on weird behaviors of M-Vicuna DupRm:
            For some reason, it requires you to specify both -opfq and -drm_op even for
            single-end input.

            Actually, more interestingly, this doesn't work at all.
        """
        _log.warn("MVicuna duplicate removal doesn't work for single-end reads. Copying input to output and allowing all reads to pass.")
        shutil.copyfile(inFastq, outFastq)
        '''
        tmp1OutPair = (
            util.file.mkstempfname(suffix='.tmp1out.1.fastq'), util.file.mkstempfname(suffix='.tmp1out.2.fastq')
        )
        tmp2OutPair = (
            util.file.mkstempfname(suffix='.tmp2out.1.fastq'), util.file.mkstempfname(suffix='.tmp2out.2.fastq')
        )
        cmdline = [
            self.install_and_get_path(), '-isfq', inFastq, '-osfq', outFastq,
            '-opfq', ','.join(tmp2OutPair), '-drm_op', ','.join(tmp1OutPair),
            '-tasks', 'DupRm'
        ]
        _log.debug(' '.join(cmdline))
        subprocess.check_call(cmdline)
        '''


def _get_mvicuna_path():
    uname = os.uname()
    if uname[0] == 'Darwin':
        osName = 'MacOSX'
    elif uname[0] == 'Linux' and uname[4].endswith('64'):
        osName = 'linux64'
    else:
        _log.debug('mvicuna not implemented for OS %s %s', uname[0], uname[4])
        return ''
    binaries_path = util.file.get_binaries_path()
    return os.path.join(binaries_path, 'mvicuna', osName, 'mvicuna')
