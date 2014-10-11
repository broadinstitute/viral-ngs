"Tools in the 'last' suite."

import tools, os

lastBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/last'

class LastTools(tools.Tool) :
	"""Abstract base class for tools in the 'last' suite.
	   Subclasses must define a class member subtoolName."""
	def __init__(self, install_methods = None):
		if install_methods == None:
			path = os.path.join(lastBroadUnixPath, self.subtoolName)
			install_methods = [tools.PrexistingUnixCommand(path)]
		super(LastTools, self).__init__(install_methods = install_methods)

class Lastal(LastTools) :
	subtoolName = 'lastal'

class MafSort(LastTools) :
	subtoolName = 'maf-sort'

class MafConvert(LastTools) :
	subtoolName = 'maf-convert'
