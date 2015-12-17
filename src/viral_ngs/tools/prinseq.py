"tools.Tool for prinseq."

import tools

TOOL_NAME = "prinseq"
TOOL_VERSION = '0.20.4'
TOOL_URL = 'http://sourceforge.net/projects/prinseq/files/standalone/' \
      'prinseq-lite-{ver}.tar.gz'.format(ver=TOOL_VERSION)


class PrinseqTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []

            install_methods.append(tools.CondaPackage(TOOL_NAME, executable="prinseq-lite.pl", version=TOOL_VERSION))

            target_rel_path = 'prinseq-lite-{ver}/prinseq-lite.pl'.format(ver=TOOL_VERSION)
            install_methods.append(
                tools.DownloadPackage(
                    TOOL_URL,
                    target_rel_path,
                    post_download_command='chmod +x {}'.format(target_rel_path),
                    require_executability=True
                )
            )
        tools.Tool.__init__(self, install_methods=install_methods)
