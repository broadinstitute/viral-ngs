"tools.Tool for prinseq."

import os
import tools, scripts

# BroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl'

class PrinseqTool(tools.Tool) :
	def __init__(self, install_methods = None):
		if install_methods == None:
			path = os.path.join(scripts.get_scripts_path(), 'prinseq-lite.pl')
			install_methods = [tools.PrexistingUnixCommand(path)]
		tools.Tool.__init__(self, install_methods = install_methods)
