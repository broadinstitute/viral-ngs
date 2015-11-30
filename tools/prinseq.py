"tools.Tool for prinseq."

import tools

tool_version = '0.20.4'
url = 'http://sourceforge.net/projects/prinseq/files/standalone/' \
      'prinseq-lite-{ver}.tar.gz'.format(ver=tool_version)


class PrinseqTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            target_rel_path = 'prinseq-lite-{ver}/prinseq-lite.pl'.format(ver=tool_version)
            install_methods.append(
                tools.DownloadPackage(url,
                                      target_rel_path,
                                      post_download_command='chmod +x {}'.format(target_rel_path),
                                      require_executability=True))
        tools.Tool.__init__(self, install_methods=install_methods)
