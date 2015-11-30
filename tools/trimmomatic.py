"tools.Tool for trimmomatic."

import tools

TOOL_URL = 'http://www.usadellab.org/cms/uploads/supplementary/' \
                 'Trimmomatic/Trimmomatic-0.32.zip'


class TrimmomaticTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.DownloadPackage(TOOL_URL,
                                                         'Trimmomatic-0.32/trimmomatic-0.32.jar',
                                                         require_executability=False))
        tools.Tool.__init__(self, install_methods=install_methods)
