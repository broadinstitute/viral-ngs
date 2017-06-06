'''
    CLARK - metagenomics tool
'''

import logging
import os
import os.path
import subprocess
import shutil
import tempfile
import random

import tools
import util.file
import util.misc

TOOL_NAME = 'clark'
TOOL_VERSION = '1.2.3'

log = logging.getLogger(__name__)

class ClarkTool(tools.Tool):
    """Tool wrapper for the CLARK metagenomics tool"""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION,
                                              executable='CLARK', env='clark_env')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args, stdout=None):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        subprocess.check_call(tool_cmd, stdout=stdout)
        if stdout:
            stdout.close()







                      
                      

    

