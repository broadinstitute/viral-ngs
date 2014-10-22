"Tools in the 'last' suite."

import tools, util.file
import os

#lastBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/last'

class LastTools(tools.Tool) :
	"""Abstract base class for tools in the 'last' suite.
	   Subclasses must define class members subtoolName and subtoolNameOnBroad."""
	def __init__(self, install_methods = None) :
		if install_methods == None :
			install_methods = []
			install_methods.append(DownloadAndBuildLast(self.subtoolName))
			#Version of last on broad is old, database-incompatible with newer one,
			#   so don't use it and always load the newer version
			#path = os.path.join(lastBroadUnixPath, self.subtoolNameOnBroad)
			#install_methods.append(tools.PrexistingUnixCommand(path))
		tools.Tool.__init__(self, install_methods = install_methods)

class DownloadAndBuildLast(tools.DownloadPackage) :
	lastWithVersion = 'last-490'
	def __init__(self, subtoolName) :
		url = 'http://last.cbrc.jp/{}.zip'.format(self.lastWithVersion)
		target_rel_path = os.path.join(self.lastWithVersion, 'bin', subtoolName)
		tools.DownloadPackage.__init__(self, url, target_rel_path)
	def post_download(self) :
		path = os.path.join(self.destination_dir, self.lastWithVersion)
		os.system('cd {}; make; make install prefix=.'.format(path))

class Lastal(LastTools) :
	subtoolName = 'lastal'
	subtoolNameOnBroad = 'lastal'

class MafSort(LastTools) :
	subtoolName = 'maf-sort'
	subtoolNameOnBroad = 'scripts/maf-sort.sh'

class Lastdb(LastTools) :
	subtoolName = 'lastdb'
	subtoolNameOnBroad = 'lastdb'

# As of version 490 of "last", maf-convert doesn't run in python 3.x.
# Workaround this by distributing our own copy that does.
"""
class MafConvert(LastTools) :
	subtoolName = 'maf-convert'
	subtoolNameOnBroad = 'scripts/maf-convert.py'
"""
class MafConvert(tools.Tool) :
	def __init__(self, install_methods = None):
		if install_methods == None:
			path = os.path.join(util.file.get_scripts_path(), 'maf-convert.last-490.2to3.py')
			install_methods = [tools.PrexistingUnixCommand(path)]
		tools.Tool.__init__(self, install_methods = install_methods)
