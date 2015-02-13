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
        destination_dir = os.path.join(
            util.file.get_build_path(), 'mosaik-{}'.format(tool_version))
        install_methods.append(
            DownloadAndBuildMosaik(url.format(ver=tool_version, os='source'),
                os.path.join(destination_dir, 'bin', 'MosaikAligner'),
                destination_dir))
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def get_networkFile(self):
        # this is the directory to return
        dir = os.path.join(util.file.get_build_path(),
            'mosaik-{}'.format(tool_version),
            'MOSAIK-{}-source'.format(tool_version),
            'networkFile')
        if not os.path.isdir(dir):
            # if it doesn't exist, run just the download-unpack portion of the
            #     source installer
            self.get_install_methods()[-1].download()
        assert os.path.isdir(dir)
        return dir

class DownloadAndBuildMosaik(tools.DownloadPackage) :
    def post_download(self) :
        mosaikDir = os.path.join(self.destination_dir,
                                 'MOSAIK-{}-source'.format(tool_version))
        if tool_version == "2.1.33" :
            # In this version, obsolete LDFLAGS breaks make. Remove it
            makeFilePath = os.path.join(mosaikDir, 'Makefile')
            os.rename(makeFilePath, makeFilePath + '.orig')
            with open(makeFilePath + '.orig') as inf :
                makeText = inf.read()
            with open(makeFilePath, 'wt') as outf :
                outf.write(makeText.replace('export LDFLAGS = -Wl',
                                            '#export LDFLAGS = -Wl'))
    
        # Now we can make:
        os.system('cd "{}" && make -s'.format(mosaikDir))
