"Tools in the 'last' suite."

# built-ins
import os
import logging
import subprocess

# within this module
import util.file
import tools

_log = logging.getLogger(__name__)

TOOL_NAME = "last"
TOOL_VERSION = "719"


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
        tools.Tool.__init__(self, install_methods=install_methods)


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

    def build_database(self, fasta_files, database_prefix_path): # pylint: disable=W0221
        output_file_prefix = os.path.basename(database_prefix_path)
        output_directory = os.path.dirname(database_prefix_path)


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

        self.execute(input_fasta, output_directory, output_file_prefix)    

        return database_prefix_path


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
        _log.debug(" ".join(tool_cmd))
        subprocess.check_call(tool_cmd)

        # restore cwd
        os.chdir(cwd_before_lastdb)


class MafConvert(LastTools):
    """ wrapper for maf-convert subtool """
    subtool_name = 'maf-convert'
    subtool_name_on_broad = 'scripts/maf-convert.py'
