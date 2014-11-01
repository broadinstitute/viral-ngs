"tools.Tool for bmtagger.sh."

import tools, util.file
import os, logging
from tools import urlretrieve
log = logging.getLogger(__name__)


class BmtaggerTool(tools.Tool) :
    '''
    Install bmtagger.sh, bmfilter, extract_fullseq, and srprism.
    
    Note: bmtagger calls blastn so that must be installed somewhere in $PATH.
    
    WARNING: bmtagger.sh does not work with the version of getopt that ships
    with Mac OS X. This can be worked around by installing linux getopt
    using fink and assuring that /sw/bin comes before /usr/bin in $PATH.

    '''
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = []
            install_methods.append(DownloadBmtagger())
            #bmtaggerBroadUnixPath = \
            #   '/idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh'
            #install_methods.append(tools.PrexistingUnixCommand(
            #    bmtaggerBroadUnixPath, require_executability=False))
        tools.Tool.__init__(self, install_methods = install_methods)
    def install_and_get_related_path(self, otherTool) :
        '''
        Get path to one of the associated tools we download with bmtagger:
            bmfilter, extract_fullseq, srprism
        '''
        self.install()
        return os.path.join(os.path.dirname(self.executable_path()), otherTool)

class DownloadBmtagger(tools.InstallMethod) :
    executables = ['bmtagger.sh', 'bmfilter', 'extract_fullseq', 'srprism']
    def __init__(self) :
        self.installed = False
        self.targetDir = os.path.join(util.file.get_build_path(), 'bmtagger')
        self.targetpath = os.path.join(self.targetDir, 'bmtagger.sh')
        tools.InstallMethod.__init__(self)
    def is_installed(self):
        return self.installed
    def executable_path(self) :
        return self.installed and self.targetpath or None
    def verify_install(self) :
        self.installed = all(os.access(os.path.join(self.targetDir, executable),
                                       (os.X_OK | os.R_OK))
                             for executable in self.executables)
        return self.installed
    def _attempt_install(self) :
        if self.verify_install() :
            return
        util.file.mkdir_p(self.targetDir)
        urlBase = 'ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/'
        uname = os.uname()
        if uname[0] == 'Darwin' :
            urlBase += 'mac-os/'
        elif uname[0] != 'Linux' or not uname[4].endswith('64') :
            log.debug('OS {} not implemented'.format(uname[0]))
            return
        for executable in self.executables :
            path = os.path.join(self.targetDir, executable)
            url = urlBase + executable
            log.info('Downloading from {} ...'.format(url))
            urlretrieve(url, path)
            os.system('chmod +x ' + path)
        self.verify_install()
