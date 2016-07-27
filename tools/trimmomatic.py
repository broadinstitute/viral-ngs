"tools.Tool for trimmomatic."

import subprocess
import tools

TOOL_NAME = "trimmomatic"
TOOL_VERSION = "0.36"


class TrimmomaticTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()] + args
        subprocess.check_call(cmd)
