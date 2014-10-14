"Tools in the 'last' suite."

import tools
import os, tempfile

#lastBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/last'

class LastTools(tools.Tool) :
	"""Abstract base class for tools in the 'last' suite.
	   Subclasses must define class members subtoolName and subtoolNameOnBroad."""
	def __init__(self, install_methods = None):
		if install_methods == None:
			install_methods = []
			#Version of last on broad is old, database-incompatible with newer one,
			#   so don't use it and always load the newer version
			#path = os.path.join(lastBroadUnixPath, self.subtoolNameOnBroad)
			#install_methods.append(tools.PrexistingUnixCommand(path))
			install_methods.append(DownloadAndBuildLast(self.subtoolName))
		tools.Tool.__init__(self, install_methods = install_methods)

class DownloadAndBuildLast(tools.DownloadPackage) :
	lastWithVersion = 'last-490'
	def __init__(self, subtoolName) :
		url = 'http://last.cbrc.jp/{}.zip'.format(self.lastWithVersion)
		targetpath = os.path.join('build', self.lastWithVersion, 'bin', subtoolName)
		download_dir = tempfile.tempdir
		unpack_dir = 'build'
		tools.DownloadPackage.__init__(self, url = url, targetpath = targetpath,
									   download_dir = download_dir, unpack_dir = unpack_dir)
	def post_download(self) :
		self.unpack()
		path = os.path.join('build', self.lastWithVersion)
		os.system('cd {}; make; make install prefix=.'.format(path))

class Lastal(LastTools) :
	subtoolName = 'lastal'
	subtoolNameOnBroad = 'lastal'

class MafSort(LastTools) :
	subtoolName = 'maf-sort'
	subtoolNameOnBroad = os.path.join('scripts', 'maf-sort.sh')

class MafConvert(LastTools) :
	subtoolName = 'maf-convert'
	subtoolNameOnBroad = os.path.join('scripts', 'maf-convert.py')
