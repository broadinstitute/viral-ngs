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
        # set env variable BLD_PLATFORM=macosx or macosx64 (default is linux)
        tool_dir = 'MOSAIK-{ver}-{os}'.format(ver=tool_version, os='source')
        self.tool_dir = os.path.join(util.file.get_build_path(), tool_dir)
        install_methods = [
            tools.DownloadPackage(url.format(ver=tool_version, os='source'),
                os.path.join(tool_dir, 'bin', 'MosaikAligner'),
                post_download_command='cd {}; make'.format(tool_dir))
        ]
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def get_networkFile(self):
        return os.path.join(self.tool_dir, 'networkFile')
