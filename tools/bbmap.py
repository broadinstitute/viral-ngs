'''
    Tool wrapper for the BBMap aligner and related tools.
'''

import logging
import os
import os.path
import subprocess

import tools

TOOL_NAME = 'bbmap'
TOOL_VERSION = '38.20'

log = logging.getLogger(__name__)

class BBMapTool(tools.Tool):
    '''Tool wrapper for the BBMap aligner and related tools.'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='bbmap.sh')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, tool, args, **kwargs):    # pylint: disable=W0221
        tool_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(tool_dir, tool)] + list(map(str, args)) + \
                   [arg if val is True else '{}={}'.format(arg, val) for arg, val in kwargs.items()]
        log.debug('Running BBMap tool: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
