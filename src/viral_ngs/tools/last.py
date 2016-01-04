"Tools in the 'last' suite."

# built-ins
import os
import logging
import subprocess

# within this module
import tools

LOG = logging.getLogger(__name__)

TOOL_NAME = "last"
TOOL_VERSION = "638"


class LastTools(tools.Tool):
    """
    "Abstract" base class for tools in the 'last' suite.
    Subclasses must define class members subtool_name #and subtool_name_on_broad.
    """

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else None
        self.subtool_name_on_broad = self.subtool_name_on_broad if hasattr(self, "subtool_name_on_broad") else None
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
            install_methods.append(DownloadAndBuildLast(self.subtool_name))
            # Version of last on broad is old, database-incompatible with newer
            # one, so don't use it and always load the newer version
            #path = os.path.join(lastBroadUnixPath, self.subtool_name_on_broad)
            # install_methods.append(tools.PrexistingUnixCommand(path))
        tools.Tool.__init__(self, install_methods=install_methods)


class DownloadAndBuildLast(tools.DownloadPackage):
    """ class for platform-specific last download and install """
    last_with_version = 'last-490'

    def __init__(self, subtool_name):
        url = 'http://last.cbrc.jp/{}.zip'.format(self.last_with_version)
        target_rel_path = os.path.join(self.last_with_version, 'bin', subtool_name)
        tools.DownloadPackage.__init__(self, url, target_rel_path)

    def post_download(self):
        path = os.path.join(self.destination_dir, self.last_with_version)
        os.system('cd {}; make -s; make -s install prefix=.'.format(path))

        # maf_convert doesn't run in Python 3.x. Fix it (in place).
        binpath = os.path.join(path, 'bin')
        maf_convert_path = os.path.join(binpath, 'maf-convert')
        os.system('2to3 {maf_convert_path} -w --no-diffs'.format(**locals()))
        # Still more fixes needed: the first for 3.x, the second for 2.7
        with open(maf_convert_path, 'rt') as inf:
            file_contents = inf.read()
            file_contents = file_contents.replace('string.maketrans("", "")', 'None')
            file_contents = file_contents.replace(
                '#! /usr/bin/env python', '#! /usr/bin/env python\nfrom __future__ import print_function'
            )
        with open(maf_convert_path, 'wt') as outf:
            outf.write(file_contents)

    def verify_install(self):
        'Default checks + verify python 2.7/3.x compatibility fixes were done'
        if not tools.DownloadPackage.verify_install(self):
            return False
        maf_convert_path = os.path.join(self.destination_dir, self.last_with_version, 'bin', 'maf-convert')
        if not os.access(maf_convert_path, os.X_OK | os.R_OK):
            return False
        with open(maf_convert_path, 'rt') as inf:
            return 'print_function' in inf.read()


class Lastal(LastTools):
    """ wrapper for lastal subtool """
    subtool_name = 'lastal'
    subtool_name_on_broad = 'lastal'


class MafSort(LastTools):
    """ wrapper for maf-sort subtool """
    subtool_name = 'maf-sort'
    subtool_name_on_broad = 'scripts/maf-sort.sh'


class Lastdb(LastTools):
    """ wrapper for lastdb subtool """
    subtool_name = 'lastdb'
    subtool_name_on_broad = 'lastdb'

    def execute(self, inputFasta, outputDirectory, outputFilePrefix):    # pylint: disable=W0221
        # get the path to the binary
        tool_cmd = [self.install_and_get_path()]

        # if the output directory (and its parents) do not exist, create them
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)

        # store the cwd because we will be changing it to the file destination
        cwd_before_lastdb = os.getcwd()

        # append the prefix given to files created by lastdb
        tool_cmd.append(outputFilePrefix)

        # append the input filepath
        tool_cmd.append(os.path.realpath(inputFasta))

        # lastdb writes files to the current working directory, so we need to set
        # it to the desired output location
        os.chdir(os.path.realpath(outputDirectory))

        # execute the lastdb command
        LOG.debug(" ".join(tool_cmd))
        subprocess.check_call(tool_cmd)

        # restore cwd
        os.chdir(cwd_before_lastdb)


class MafConvert(LastTools):
    """ wrapper for maf-convert subtool """
    subtool_name = 'maf-convert'
    subtool_name_on_broad = 'scripts/maf-convert.py'
