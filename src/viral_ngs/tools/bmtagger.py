"tools.Tool for bmtagger.sh."

import tools
import util.file
import os
import logging
from tools import urlretrieve
_log = logging.getLogger(__name__)

TOOL_NAME = "bmtagger"
TOOL_VERSION = "3.101"


class BmtaggerTools(tools.Tool):
    '''
    "Abstract" base class for bmtagger.sh, bmfilter, extract_fullseq, srprism.
    Subclasses must define class member subtool_name.

    Note: bmtagger calls blastn so that must be installed somewhere in $PATH.

    WARNING: bmtagger.sh does not work with the version of getopt that ships
    with Mac OS X. This can be worked around by installing linux getopt
    using fink and assuring that /sw/bin comes before /usr/bin in $PATH.

    '''

    # subtool_name must be defined in subclass

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "bmtagger.sh"
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
            install_methods.append(DownloadBmtagger(self.subtool_name))
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        util.misc.run_and_print(cmd, check=True)


class BmtaggerShTool(BmtaggerTools):
    """ tool wrapper for bmtagger """
    subtool_name = 'bmtagger.sh'


class BmfilterTool(BmtaggerTools):
    """ tool wrapper for bmfilter """
    subtool_name = 'bmfilter'


class BmtoolTool(BmtaggerTools):
    """ tool wrapper for bmtool """
    subtool_name = 'bmtool'

    def build_database(self, fasta_files, bitmask_file_path, max_amig=0, word_size=8):
        """ builds a bmtool database (*.bitmask file) """

        input_fasta = ""

        # we can pass in a string containing a fasta file path
        # or a list of strings
        if 'basestring' not in globals():
           basestring = str
        if isinstance(fasta_files, basestring):
            fasta_files = [fasta_files]
        elif isinstance(fasta_files, list):
            pass
        else:
            raise TypeError("fasta_files was not a single fasta file, nor a list of fasta files") # or something along that line

        # if more than one fasta file is specified, join them
        # otherwise if only one is specified, just use it
        if len(fasta_files) > 1:
            input_fasta = util.file.mkstempfname("fasta")
            util.file.cat(input_fasta, fasta_files)
        elif len(fasta_files) == 1:
            input_fasta = fasta_files[0]
        else:
            raise IOError("No fasta file provided")

        args = ['-d', input_fasta, '-o', bitmask_file_path, '-A', str(max_amig), '-w', str(word_size)]
        self.execute(*args)

        return bitmask_file_path

class ExtractFullseqTool(BmtaggerTools):
    """ tool wrapper for extract_fullseq """
    subtool_name = 'extract_fullseq'


class SrprismTool(BmtaggerTools):
    """ tool wrapper for srprism """
    subtool_name = 'srprism'

    def build_database(self, fasta_files, database_prefix_path):
        """ builds a srprism database """

        input_fasta = ""

        # we can pass in a string containing a fasta file path
        # or a list of strings
        if 'basestring' not in globals():
           basestring = str
        if isinstance(fasta_files, basestring):
            fasta_files = [fasta_files]
        elif isinstance(fasta_files, list):
            pass
        else:
            raise TypeError("fasta_files was not a single fasta file, nor a list of fasta files") # or something along that line

        # if more than one fasta file is specified, join them
        # otherwise if only one is specified, just use it
        if len(fasta_files) > 1:
            input_fasta = util.file.mkstempfname("fasta")
            util.file.cat(input_fasta, fasta_files)
        elif len(fasta_files) == 1:
            input_fasta = fasta_files[0]
        else:
            raise IOError("No fasta file provided")

        args = ['mkindex', '-i', input_fasta, '-o', database_prefix_path]
        self.execute(*args)

        return database_prefix_path


class DownloadBmtagger(tools.InstallMethod):
    """ InstallMethod class for downloading platform-specific bmtagger """
    executables = ['bmtagger.sh', 'bmfilter', 'extract_fullseq', 'srprism']

    def __init__(self, subtool_name):
        self.installed = False
        self.target_dir = os.path.join(util.file.get_build_path(), 'bmtagger')
        self.target_path = os.path.join(self.target_dir, subtool_name)
        tools.InstallMethod.__init__(self)

    def is_installed(self):
        return self.installed

    def executable_path(self):
        return self.installed and self.target_path or None

    def verify_install(self):
        """ confirms that the tools are present and executable """
        self.installed = all(
            os.access(
                os.path.join(self.target_dir, executable), (os.X_OK | os.R_OK)) for executable in self.executables
        )
        return self.installed

    def _attempt_install(self):
        if self.verify_install():
            return
        util.file.mkdir_p(self.target_dir)
        url_base = 'ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/'
        uname = os.uname()
        if uname[0] == 'Darwin':
            url_base += 'mac-os/'
        elif uname[0] != 'Linux' or not uname[4].endswith('64'):
            _log.debug('OS %s not implemented', uname[0])
            return
        for executable in self.executables:
            path = os.path.join(self.target_dir, executable)
            url = url_base + executable
            _log.info('Downloading from %s ...', url)
            urlretrieve(url, path)
            os.system('chmod +x ' + path)
        self.verify_install()
