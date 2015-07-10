"Tools in the 'last' suite."

# built-ins
import os, logging, subprocess

# within this module
import tools

log = logging.getLogger(__name__)

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
        os.system('cd {}; make -s; make -s install prefix=.'.format(path))

        # maf_convert doesn't run in Python 3.x. Fix it (in place).
        binpath = os.path.join(path, 'bin')
        mafConvertPath = os.path.join(binpath, 'maf-convert')
        os.system('2to3 {mafConvertPath} -w --no-diffs'.format(**locals()))
        # Still more fixes needed: the first for 3.x, the second for 2.7
        with open(mafConvertPath, 'rt') as inf:
            fileContents = inf.read()
            fileContents = fileContents.replace('string.maketrans("", "")', 'None')
            fileContents = fileContents.replace(
                '#! /usr/bin/env python',
                '#! /usr/bin/env python\nfrom __future__ import print_function')
        with open(mafConvertPath, 'wt') as outf:
            outf.write(fileContents)
    def verify_install(self) :
        'Default checks + verify python 2.7/3.x compatibility fixes were done'
        if not tools.DownloadPackage.verify_install(self) :
            return False
        mafConvertPath = os.path.join(self.destination_dir,
            self.lastWithVersion, 'bin', 'maf-convert')
        if not os.access(mafConvertPath, os.X_OK | os.R_OK) :
            return False
        with open(mafConvertPath, 'rt') as inf:
            return 'print_function' in inf.read()

class Lastal(LastTools) :
    subtoolName = 'lastal'
    subtoolNameOnBroad = 'lastal'

class MafSort(LastTools) :
    subtoolName = 'maf-sort'
    subtoolNameOnBroad = 'scripts/maf-sort.sh'

class Lastdb(LastTools) :
    subtoolName = 'lastdb'
    subtoolNameOnBroad = 'lastdb'

    def execute(self, inputFasta, outputDirectory, outputFilePrefix):
        # get the path to the binary
        toolCmd = [self.install_and_get_path()]
        
        # if the output directory (and its parents) do not exist, create them
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)

        # store the cwd because we will be changing it to the file destination
        cwdBeforeLastdb = os.getcwd()

        # lastdb writes files to the current working directory, so we need to set 
        # it to the desired output location
        os.chdir(os.path.realpath(outputDirectory))

        # append the prefix given to files created by lastdb
        toolCmd.append(outputFilePrefix)

        # append the input filepath
        toolCmd.append(inputFasta)

        #execute the lastdb command
        log.debug(" ".join(toolCmd))
        subprocess.check_call(toolCmd)

        # restore cwd
        os.chdir(cwdBeforeLastdb)

class MafConvert(LastTools) :
    subtoolName = 'maf-convert'
    subtoolNameOnBroad = 'scripts/maf-convert.py'


