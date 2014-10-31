"Tools in the blast+ suite."
import tools
import os

urlPrefix = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables' \
            '/blast+/2.2.29/ncbi-blast-2.2.29+-'

def get_url() :
    uname = os.uname()
    if uname[0] == 'Darwin' :
        osStr = 'universal-macosx'
    elif uname[0] == 'Linux' :
        if uname[4].endswith('64') :
            osStr = 'x64-linux'
        else :
            osStr = 'ia32-linux'
    else :
        raise NotImplementedError('OS {} not implemented'.format(uname[0]))
    return urlPrefix + osStr + '.tar.gz'

class BlastTools(tools.Tool) :
    """'Abstract' base class for tools in the blast+ suite.
       Subclasses must define class member subtoolName."""
    def __init__(self, install_methods = None) :
        if install_methods == None :
            target_rel_path = 'ncbi-blast-2.2.29+/bin/' + self.subtoolName
            install_methods = [tools.DownloadPackage(get_url(), target_rel_path)]
        tools.Tool.__init__(self, install_methods = install_methods)

class BlastnTool(BlastTools) :
    subtoolName = 'blastn'
