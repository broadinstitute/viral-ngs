"tools.Tool for bmtagger.sh."

import tools
import util.file
import os
import logging
import subprocess
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
            install_methods = [tools.PrexistingUnixCommand(shutil.which(self.subtool_name), require_executability=False)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        subprocess.check_call(cmd)

    def silent_execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        with open(os.devnull, 'w') as fnull:
            subprocess.check_call(cmd, stderr=fnull)


class BmtaggerShTool(BmtaggerTools):
    """ tool wrapper for bmtagger """
    subtool_name = 'bmtagger.sh'


class BmfilterTool(BmtaggerTools):
    """ tool wrapper for bmfilter """
    subtool_name = 'bmfilter'


class BmtoolTool(BmtaggerTools):
    """ tool wrapper for bmtool """
    subtool_name = 'bmtool'

    def build_database(self, fasta_files, bitmask_file_path, max_ambig=0, word_size=18):
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

        args = ['-d', input_fasta, '-o', bitmask_file_path, '-A', str(max_ambig), '-w', str(word_size)]
        self.silent_execute(*args)

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
