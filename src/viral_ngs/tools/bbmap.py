'''
    Tool wrapper for the BBMap aligner and related tools.
'''

import logging
import os
import os.path
import shutil
import subprocess

import util.file
import tools
import tools.samtools
import tools.picard

TOOL_NAME = 'bbmap.sh'

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class BBMapTool(tools.Tool):
    '''Tool wrapper for the BBMap aligner and related tools.'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(BBMapTool, self).__init__(install_methods=install_methods)

    def version(self):
        return subprocess.check_output([os.path.join(os.path.dirname(self.install_and_get_path()), 'bbversion.sh')]).decode('UTF-8').strip()

    def execute(self, tool, **kwargs):  # pylint: disable=arguments-differ
        tool_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(tool_dir, tool)] + \
                   ['{}={}'.format('in' if arg=='in_' else arg,
                                   (val is True and 't') or (val is False and 'f') or val)
                    for arg, val in kwargs.items()]
        _log.debug('Running BBMap tool: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def align(self, inBam, refFasta, outBam, min_qual=0, nodisk=True, JVMmemory=None, **kwargs):
        with tools.samtools.SamtoolsTool().bam2fq_tmp(inBam) as (in1, in2), \
             util.file.tmp_dir('_bbmap_align') as t_dir:
            tmp_bam = os.path.join(t_dir, 'bbmap_out.bam')
            self.execute(tool='bbmap.sh', in1=in1, in2=in2, ref=refFasta, out=tmp_bam, nodisk=nodisk, **kwargs)
            
            # Samtools filter (optional)
            if min_qual:
                tmp_bam2 = os.path.join(tdir, 'bbmap.filtered.bam')
                cmd = [samtools.install_and_get_path(), 'view', '-b', '-S', '-1', '-q', str(min_qual), tmp_bam]
                _log.debug('%s > %s', ' '.join(cmd), tmp_bam2)
                with open(tmp_bam2, 'wb') as outf:
                    util.misc.run_and_save(cmd, outf=outf)
                os.unlink(tmp_bam)
                tmp_bam = tmp_bam2

            # Picard SortSam
            sorter = tools.picard.SortSamTool()
            sorter.execute(
                tmp_bam,
                outBam,
                sort_order='coordinate',
                picardOptions=['CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'],
                JVMmemory=JVMmemory
            )

