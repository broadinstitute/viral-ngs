"tools.Tool for bmtagger.sh."

import tools, util.file
import os

bmtaggerBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/bmtagger/bmtagger.sh'
bmtaggerURL = 'ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmtagger.sh'
#    Note that as of 2014-10-21, the version of bmtagger.sh in the mac-os directory
#    is the same as the top level, so no need to distinguish target os.

class BmtaggerTool(tools.Tool) :
	def __init__(self, install_methods = None) :
		if install_methods == None :
			install_methods = []
			install_methods.append(tools.PrexistingUnixCommand(bmtaggerBroadUnixPath,
															   require_executability=False))
			install_methods.append(tools.DownloadScript(bmtaggerURL, 'bmtagger.sh'))
		tools.Tool.__init__(self, install_methods = install_methods)
