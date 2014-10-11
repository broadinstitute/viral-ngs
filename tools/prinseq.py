"tools.Tool for prinseq."

import tools

prinseqBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl'

class PrinseqTool(tools.Tool) :
	def __init__(self, install_methods = None):
		if install_methods == None:
			install_methods = [tools.PrexistingUnixCommand(prinseqBroadUnixPath)]
		super(PrinseqTool, self).__init__(install_methods = install_methods)
