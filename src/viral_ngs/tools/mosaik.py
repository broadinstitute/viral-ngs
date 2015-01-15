'''
    The MOSAIK aligner
'''

import logging, tools, util.file
import os, os.path, subprocess

tool_version = '2.1.33'
url = 'https://mosaik-aligner.googlecode.com/files/MOSAIK-{ver}-{os}.tar'

log = logging.getLogger(__name__)

class MosaikTool(tools.Tool) :
    def __init__(self) :
        if os.uname()[0] == 'Darwin':
            if os.uname()[4].endswith('_64'):
                os.environ['BLD_PLATFORM'] = 'macosx64'
            else:
                os.environ['BLD_PLATFORM'] = 'macosx'
        install_methods = []
        install_methods.append(tools.DownloadPackage(url.format(ver=tool_version, os='source'),
            os.path.join('bin', 'MosaikAligner'),
            post_download_command='cd MOSAIK-{}-source; make -s'.format(tool_version)))
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def get_networkFile(self):
        # this is the directory to return
        dir = os.path.join(util.file.get_build_path(),
            'MOSAIK-{}-source'.format(tool_version),
            'networkFile')
        if not os.path.isdir(dir):
            # if it doesn't exist, run just the download-unpack portion of the source installer
            self.get_install_methods()[-1].download()
        return dir
