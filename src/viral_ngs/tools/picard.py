"Tools in the 'picard' suite."
import tools

url = 'https://github.com/broadinstitute/picard/releases/download/1.123/picard-tools-1.123.zip'
#         Note: Version 1.123 is latest as of 2014-10-22
picardBroadUnix = '/seq/software/picard/current/bin'

class PicardTool(tools.Tool) :
    """'Abstract' base class for tools in the picard suite.
       Subclasses must define class member subtoolName=."""
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = []
            target_rel_path = 'picard-tools-1.123/{}'.format(self.subtoolName)
            install_methods.append(tools.DownloadPackage(url, target_rel_path,
                                                         require_executability=False))
            #broadUnixPath = picardBroadUnix + '/' + subtoolName
            #install_methods.append(tools.PrexistingUnixCommand(broadUnixPath,
            #                                                  require_executability=False))
        tools.Tool.__init__(self, install_methods = install_methods)

class MarkDuplicatesTool(PicardTool) :
    subtoolName = 'MarkDuplicates.jar'

class SamToFastqTool(PicardTool) :
    subtoolName = 'SamToFastq.jar'

class SortSamTool(PicardTool) :
    subtoolName = 'SortSam.jar'
