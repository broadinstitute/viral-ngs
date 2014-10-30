import tools, util.file
import os

url = 'http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2'

class SamtoolsTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = []
            install_methods.append(
                tools.DownloadPackage(url, 'samtools-0.1.19/samtools',
                                      post_download_command = 'cd samtools-0.1.19; make'))
            #path = '/idi/sabeti-data/software/samtools/samtools-0.1.19/samtools',
            #install_methods.append(tools.PrexistingUnixCommand(path))
        tools.Tool.__init__(self, install_methods = install_methods)

