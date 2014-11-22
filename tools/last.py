"Tools in the 'last' suite."

import tools, util.file
import os

#lastBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/last'

class LastTools(tools.Tool) :
    """
    "Abstract" base class for tools in the 'last' suite.
    Subclasses must define class members subtoolName #and subtoolNameOnBroad.
    """
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = []
            install_methods.append(DownloadAndBuildLast(self.subtoolName))
            #Version of last on broad is old, database-incompatible with newer
            # one, so don't use it and always load the newer version
            #path = os.path.join(lastBroadUnixPath, self.subtoolNameOnBroad)
            #install_methods.append(tools.PrexistingUnixCommand(path))
        tools.Tool.__init__(self, install_methods = install_methods)

class DownloadAndBuildLast(tools.DownloadPackage) :
    lastWithVersion = 'last-490'
    def __init__(self, subtoolName) :
        url = 'http://last.cbrc.jp/{}.zip'.format(self.lastWithVersion)
        target_rel_path = os.path.join(self.lastWithVersion, 'bin', subtoolName)
        tools.DownloadPackage.__init__(self, url, target_rel_path)
    def post_download(self) :
        path = os.path.join(self.destination_dir, self.lastWithVersion)
        os.system('cd {}; make; make install prefix=.'.format(path))

        # maf_convert doesn't run in Python 3.x. Fix it (in place).
        binpath = os.path.join(path, 'bin')
        mafConvertPath = os.path.join(binpath, 'maf-convert')
        os.system('2to3 {mafConvertPath} -W'.format(**locals()))
        # Still more fixes needed, the first for 3.x, the second for 2.7
        fileContents = open(mafConvertPath).read()
        fileContents = fileContents.replace('string.maketrans("", "")', 'None')
        fileContents = fileContents.replace(
            '#! /usr/bin/env python',
            '#! /usr/bin/env python\nfrom __future__ import print_function')
        open(mafConvertPath, 'w').write(fileContents)
    def verify_install(self) :
        'Default checks + verify python 2.7/3.x compatibility fixes were done'
        if not tools.DownloadPackage.verify_install(self) :
            return False
        mafConvertPath = os.path.join(self.destination_dir,
            self.lastWithVersion, 'bin', 'maf-convert')
        return 'print_function' in open(mafConvertPath).read()

class Lastal(LastTools) :
    subtoolName = 'lastal'
    subtoolNameOnBroad = 'lastal'

class MafSort(LastTools) :
    subtoolName = 'maf-sort'
    subtoolNameOnBroad = 'scripts/maf-sort.sh'

class Lastdb(LastTools) :
    subtoolName = 'lastdb'
    subtoolNameOnBroad = 'lastdb'

class MafConvert(LastTools) :
    subtoolName = 'maf-convert'
    subtoolNameOnBroad = 'scripts/maf-convert.py'


