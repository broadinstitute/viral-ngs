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
        tool_dir = 'MOSAIK-{ver}-{os}'.format(ver=tool_version, os='source')
        install_methods = [
            tools.DownloadPackage(url.format(ver=tool_version, os='source'),
                os.path.join(tool_dir, 'bin', 'MosaikAligner'),
                post_download_command='cd {}; make -s'.format(tool_dir))]
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def get_networkFile(self):
        dir = os.path.dirname(os.path.dirname(self.install_and_get_path()))
        return os.path.join(dir, 'networkFile')
