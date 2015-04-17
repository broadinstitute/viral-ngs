"tools.Tool for prinseq."

import tools, util.file

url = 'http://sourceforge.net/projects/prinseq/files/standalone/' \
      'prinseq-lite-0.19.3.tar.gz'

class PrinseqTool(tools.Tool) :
    def __init__(self, install_methods = None):
        if install_methods == None:
            install_methods = []
            target_rel_path = 'prinseq-lite-0.19.3/prinseq-lite.pl'
            install_methods.append(
                tools.DownloadPackage(url, target_rel_path,
                    post_download_command='chmod +x {}'.format(target_rel_path),
                    require_executability = True))
        tools.Tool.__init__(self, install_methods = install_methods)



