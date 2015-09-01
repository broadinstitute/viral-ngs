"tools.Tool for bmtagger.sh."

import tools
import util.file
import os
import logging
from tools import urlretrieve
log = logging.getLogger(__name__)


class BmtaggerTools(tools.Tool):
    '''
    "Abstract" base class for bmtagger.sh, bmfilter, extract_fullseq, srprism.
    Subclasses must define class member subtoolName.

    Note: bmtagger calls blastn so that must be installed somewhere in $PATH.

    WARNING: bmtagger.sh does not work with the version of getopt that ships
    with Mac OS X. This can be worked around by installing linux getopt
    using fink and assuring that /sw/bin comes before /usr/bin in $PATH.

    '''

    # subtoolName must be defined in subclass

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(DownloadBmtagger(self.subtoolName))
        tools.Tool.__init__(self, install_methods=install_methods)


class BmtaggerShTool(BmtaggerTools):
    subtoolName = 'bmtagger.sh'


class BmfilterTool(BmtaggerTools):
    subtoolName = 'bmfilter'


class Extract_fullseqTool(BmtaggerTools):
    subtoolName = 'extract_fullseq'


class SrprismTool(BmtaggerTools):
    subtoolName = 'srprism'


class DownloadBmtagger(tools.InstallMethod):
    executables = ['bmtagger.sh', 'bmfilter', 'extract_fullseq', 'srprism']

    def __init__(self, subtoolName):
        self.installed = False
        self.targetDir = os.path.join(util.file.get_build_path(), 'bmtagger')
        self.targetpath = os.path.join(self.targetDir, subtoolName)
        tools.InstallMethod.__init__(self)

    def is_installed(self):
        return self.installed

    def executable_path(self):
        return self.installed and self.targetpath or None

    def verify_install(self):
        self.installed = all(os.access(os.path.join(self.targetDir, executable), (os.X_OK | os.R_OK))
                             for executable in self.executables)
        return self.installed

    def _attempt_install(self):
        if self.verify_install():
            return
        util.file.mkdir_p(self.targetDir)
        urlBase = 'ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/'
        uname = os.uname()
        if uname[0] == 'Darwin':
            urlBase += 'mac-os/'
        elif uname[0] != 'Linux' or not uname[4].endswith('64'):
            log.debug('OS {} not implemented'.format(uname[0]))
            return
        for executable in self.executables:
            path = os.path.join(self.targetDir, executable)
            url = urlBase + executable
            log.info('Downloading from %s ...', url)
            urlretrieve(url, path)
            os.system('chmod +x ' + path)
        self.verify_install()
