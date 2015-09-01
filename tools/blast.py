"Tools in the blast+ suite."
import tools
import os

urlPrefix = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables' \
            '/blast+/2.2.29/ncbi-blast-2.2.29+-'


def get_url():
    uname = os.uname()
    if uname[0] == 'Darwin':
        osStr = 'universal-macosx'
    elif uname[0] == 'Linux':
        if uname[4].endswith('64'):
            osStr = 'x64-linux'
        else:
            osStr = 'ia32-linux'
    else:
        raise NotImplementedError('OS {} not implemented'.format(uname[0]))
    return urlPrefix + osStr + '.tar.gz'


class BlastTools(tools.Tool):
    """'Abstract' base class for tools in the blast+ suite.
       Subclasses must define class member subtoolName."""

    def __init__(self, install_methods=None):
        unwanted = ['blast_formatter', 'blastdb_aliastool', 'blastdbcheck', 'blastdbcmd', 'convert2blastmask',
                    'deltablast', 'legacy_blast.pl', 'makembindex', 'makeprofiledb', 'psiblast', 'rpsblast',
                    'rpstblastn', 'segmasker', 'tblastn', 'tblastx', 'update_blastdb.pl', 'windowmasker']
        if install_methods is None:
            target_rel_path = 'ncbi-blast-2.2.29+/bin/' + self.subtoolName
            install_methods = [tools.DownloadPackage(get_url(),
                                                     target_rel_path,
                                                     post_download_command=' '.join(
                                                         ['rm'] + ['ncbi-blast-2.2.29+/bin/' + f for f in unwanted]),
                                                     post_download_ret=None)]
        tools.Tool.__init__(self, install_methods=install_methods)


class BlastnTool(BlastTools):
    subtoolName = 'blastn'


class MakeblastdbTool(BlastTools):
    subtoolName = 'makeblastdb'
