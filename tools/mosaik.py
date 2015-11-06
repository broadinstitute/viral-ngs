'''
    The MOSAIK aligner
'''

import logging
import os
import os.path
import tools
import util.file

tool_version = '2.2.26'
commit_hash = 'e04c806bb1410cf1dbd1534991c46d696aec6723'
url = 'https://github.com/wanpinglee/MOSAIK/archive/{commit_hash}.zip'

log = logging.getLogger(__name__)


class MosaikTool(tools.Tool):

    def __init__(self):
        os.environ['BLD_PLATFORM'] = get_build_env()
        install_methods = []
        destination_dir = os.path.join(util.file.get_build_path(), 'mosaik-{}'.format(commit_hash))
        install_methods.append(
            DownloadAndBuildMosaik(url.format(commit_hash=commit_hash,
                                              os='source'),
                                   os.path.join(destination_dir, 'MOSAIK-{}'.format(commit_hash), 'bin', 'MosaikAligner'),
                                   destination_dir))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return commit_hash

    def get_networkFile(self):
        # this is the directory to return
        mosaikDir = os.path.join(util.file.get_build_path(), 'mosaik-{}'.format(commit_hash), 'MOSAIK-{}'.format(commit_hash) , 'src', 'networkFile')
        if not os.path.isdir(mosaikDir):
            # if it doesn't exist, run just the download-unpack portion of the
            #     source installer
            self.get_install_methods()[-1].download()
        assert os.path.isdir(mosaikDir)
        return mosaikDir


class DownloadAndBuildMosaik(tools.DownloadPackage):

    def post_download(self):
        mosaikDir = os.path.join(self.destination_dir, 'MOSAIK-{}/src'.format(commit_hash))

        incFilePath = os.path.join(mosaikDir, 'includes', get_build_env()+".inc")
        os.rename(incFilePath, incFilePath + '.orig')        
        with open(incFilePath + '.orig') as inf:
            makeText = inf.read()
        with open(incFilePath, 'wt') as outf:
            outf.write(makeText.replace('-static', ''))

        # Now we can make:
        os.system('cd "{}" && make -s'.format(mosaikDir))

def get_build_env():
    if os.uname()[0] == 'Darwin':
        if os.uname()[4].endswith('_64'):
            return 'macosx64'
        else:
            return 'macosx'
    else:
        return "linux"