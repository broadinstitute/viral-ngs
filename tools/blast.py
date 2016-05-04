"Tools in the blast+ suite."
import tools
import os

URL_PREFIX = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables' \
            '/blast+/2.2.29/ncbi-blast-2.2.29+-'

TOOL_NAME = "blast"
TOOL_VERSION = "2.2.31"


def get_url():
    """ creates the download url for this tool """
    uname = os.uname()
    if uname[0] == 'Darwin':
        os_str = 'universal-macosx'
    elif uname[0] == 'Linux':
        if uname[4].endswith('64'):
            os_str = 'x64-linux'
        else:
            os_str = 'ia32-linux'
    else:
        raise NotImplementedError('OS {} not implemented'.format(uname[0]))
    return URL_PREFIX + os_str + '.tar.gz'


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
            target_rel_path = 'ncbi-blast-2.2.29+/bin/' + self.subtool_name
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
            install_methods.append(
                tools.DownloadPackage(
                    get_url(),
                    target_rel_path,
                    post_download_command=' '.join(
                        ['rm'] + [
                            'ncbi-blast-2.2.29+/bin/' + f for f in unwanted
                        ]
                    ),
                    post_download_ret=None
                )
            )
        #tools.Tool.__init__(self, install_methods=install_methods)
        super(BlastTools, self).__init__(install_methods=install_methods)

    def execute(self, *args):
        util.misc.run_and_print(self.exec_path, args, check=True)


class BlastnTool(BlastTools):
    """ Tool wrapper for blastn """
    subtool_name = 'blastn'


class MakeblastdbTool(BlastTools):
    """ Tool wrapper for makeblastdb """
    subtool_name = 'makeblastdb'
