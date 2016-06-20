"Tools in the blast+ suite."
import tools
import os
import util.misc

TOOL_NAME = "blast"
TOOL_VERSION = "2.2.31"


class BlastTools(tools.Tool):
    """'Abstract' base class for tools in the blast+ suite.
       Subclasses must define class member subtool_name."""

    def __init__(self, install_methods=None):
        unwanted = [
            'blast_formatter', 'blastdb_aliastool', 'blastdbcheck', 'blastdbcmd', 'convert2blastmask', 'deltablast',
            'legacy_blast.pl', 'makembindex', 'makeprofiledb', 'psiblast', 'rpsblast', 'rpstblastn', 'segmasker',
            'tblastn', 'tblastx', 'update_blastdb.pl', 'windowmasker'
        ]
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "blastn"
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
        super(BlastTools, self).__init__(install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        util.misc.run_and_print(cmd, buffered=True, check=True)


class BlastnTool(BlastTools):
    """ Tool wrapper for blastn """
    subtool_name = 'blastn'


class MakeblastdbTool(BlastTools):
    """ Tool wrapper for makeblastdb """
    subtool_name = 'makeblastdb'

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

        args = ['-dbtype', 'nucl', '-in', input_fasta, '-out', database_prefix_path]
        self.execute(*args)

        return database_prefix_path
