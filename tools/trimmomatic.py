"tools.Tool for trimmomatic."

import tools
import util.misc

TOOL_NAME = "trimmomatic"
TOOL_VERSION = "0.35"

TOOL_URL = 'http://www.usadellab.org/cms/uploads/supplementary/' \
                 'Trimmomatic/Trimmomatic-0.32.zip'


class TrimmomaticTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
            install_methods.append(
                tools.DownloadPackage(
                    TOOL_URL,
                    'Trimmomatic-0.32/trimmomatic-0.32.jar',
                    require_executability=False
                )
            )
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, *args):
        util.misc.run_and_print(self.exec_path, args, check=True)
