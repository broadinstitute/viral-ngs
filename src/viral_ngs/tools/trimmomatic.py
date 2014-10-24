"tools.Tool for trimmomatic."

import tools, util.file
import os, tempfile


trimmomaticBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/trimmomatic/trimmomatic-0.32.jar'

trimmomaticURL = 'http://www.usadellab.org/cms/uploads/supplementary/' \
				 'Trimmomatic/Trimmomatic-0.32.zip'

class TrimmomaticTool(tools.Tool) :
	def __init__(self, install_methods = None) :
		if install_methods == None :
			install_methods = []
			install_methods.append(tools.DownloadPackage(trimmomaticURL,
														 'Trimmomatic-0.32/trimmomatic-0.32.jar',
														 require_executability=False))
			#install_methods.append(tools.PrexistingUnixCommand(trimmomaticBroadUnixPath,
			#												   require_executability=False))
		tools.Tool.__init__(self, install_methods = install_methods)
