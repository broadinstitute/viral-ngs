"tools.Tool for prinseq."

import tools

TOOL_NAME = "prinseq"
TOOL_VERSION = '0.20.4'

class PrinseqTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable="prinseq-lite.pl", version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)
