"tools.Tool for trimmomatic."

import logging
import os
import subprocess
import tools
import util.file
import util.misc

TOOL_NAME = "trimmomatic"
TOOL_VERSION = "0.36"

_log = logging.getLogger(__name__)

class TrimmomaticTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self,
        inFastq1,
        inFastq2,
        pairedOutFastq1,
        pairedOutFastq2,
        clipFasta,
        unpairedOutFastq1=None,
        unpairedOutFastq2=None,
        leading_q_cutoff=15,
        trailing_q_cutoff=15,
        minlength_to_keep=30,
        sliding_window_size=4,
        sliding_window_q_cutoff=25
    ):
        '''Trim read sequences with Trimmomatic.'''
        trimmomaticPath = self.install_and_get_path()
        unpairedFastq1 = unpairedOutFastq1 or util.file.mkstempfname()
        unpairedFastq2 = unpairedOutFastq2 or util.file.mkstempfname()
        javaCmd = [trimmomaticPath]

        if inFastq2 is None:
            # Unpaired reads
            javaCmd.extend([
                    'SE', '-phred33',
                    inFastq1, unpairedFastq1
                ])
        else:
            # Paired reads
            javaCmd.extend([
                    'PE', '-phred33',
                    inFastq1, inFastq2,
                    pairedOutFastq1, unpairedFastq1, pairedOutFastq2, unpairedFastq2
                ])

        # all the options
        javaCmd.extend(
            [
                'LEADING:{leading_q_cutoff}'.format(leading_q_cutoff=leading_q_cutoff),
                'TRAILING:{trailing_q_cutoff}'.format(trailing_q_cutoff=trailing_q_cutoff),
                'SLIDINGWINDOW:{sliding_window_size}:{sliding_window_q_cutoff}'.format(
                    sliding_window_size=sliding_window_size,
                    sliding_window_q_cutoff=sliding_window_q_cutoff,
                ), 
                'MINLEN:{minlength_to_keep}'.format(minlength_to_keep=minlength_to_keep),
                'ILLUMINACLIP:{clipFasta}:2:30:12'.format(clipFasta=clipFasta)
            ]
        )

        _log.debug(' '.join(javaCmd))
        util.misc.run_and_print(javaCmd, check=True)

        if not unpairedOutFastq1:
            os.unlink(unpairedFastq1)
        if not unpairedOutFastq2:
            os.unlink(unpairedFastq2)
