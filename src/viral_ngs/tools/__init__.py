'''Stuff here'''

__author__ = "dpark@broadinstitute.org"

import os, logging, tempfile
import util.file

try:
	# Python 3.x
	from urllib.request import urlretrieve
	from urllib.parse import urlparse
except ImportError:
	# Python 2.x
	from urllib import urlretrieve
	from urlparse import urlparse

__all__ = ['snpeff','prinseq','last','trimmomatic','bmtagger','samtools','picard', 'blast']
installed_tools = {}

log = logging.getLogger(__name__)

def get_tool_by_name(name):
	if name not in installed_tools:
		pass
		raise NotImplementedError
	return installed_tools[name]

class Tool(object):
	''' Base tool class that includes install machinery.
		
		TO DO: add something about dependencies..
	'''
	def __init__(self, install_methods=[]):
		self.install_methods = install_methods
		self.installed_method = None
		self.exec_path = None
	def is_installed(self):
		return (self.installed_method != None)
	def install(self):
		if not self.is_installed():
			for m in self.install_methods:
				if not m.is_attempted():
					m.attempt_install()
				if m.is_installed():
					self.installed_method = m
					self.exec_path = m.executable_path()
					installed_tools[self.__class__.__name__] = self
					break
	def get_install_methods(self):
		return self.install_methods
	def set_install_methods(self, methods):
		self.install_methods = methods
	def version(self):
		return None
	def executable_path(self):
		return self.exec_path
	def execute(self, args):
		assert not os.system(self.exec_path + ' ' + args)
	def install_and_get_path(self) :
		self.install()
		return self.executable_path()

class InstallMethod(object):
	''' Base class for installation methods for a given tool.
		None of these methods should ever fail/error. attempt_install should
		return silently regardless of the outcome (is_installed must be
		called to verify success or failure).
	'''
	def __init__(self):
		self.attempts = 0
	def is_attempted(self):
		return self.attempts
	def attempt_install(self):
		self.attempts += 1
	def is_installed(self):
		raise NotImplementedError
	def executable_path(self):
		raise NotImplementedError

class PrexistingUnixCommand(InstallMethod):
	''' This is an install method that tries to find whether an executable binary
		already exists for free on the unix file system--it doesn't actually try to
		install anything.
	'''
	def __init__(self, path, verifycmd=None, verifycode=0, require_executability=True):
		self.path = path
		self.verifycmd = verifycmd
		self.verifycode = verifycode
		self.attempted = False
		self.installed = False
		self.require_executability = require_executability
	def is_attempted(self):
		return self.attempted
	def attempt_install(self):
		self.attempted = True
		if os.access(self.path,
					 (os.X_OK | os.R_OK) if self.require_executability else os.R_OK):
			if self.verifycmd:
				self.installed = (os.system(self.verifycmd) == self.verifycode)
			else:	
				self.installed = True
		else:
			self.installed = False
	def is_installed(self):
		if not self.is_attempted():
			self.attempt_install()
		return self.installed
	def executable_path(self):
		return self.installed and self.path or None

class DownloadPackage(InstallMethod):
	''' This is an install method for downloading, unpacking, and post-processing
		    something straight from the source.
		target_rel_path is the path of the executable relative to destination_dir
		destination_dir defaults to the project build directory
		post_download_command will be executed if it isn't None, in destination_dir.
	'''
	def __init__(self, url, target_rel_path, destination_dir=None,
				 verifycmd=None, verifycode=0, require_executability=True,
				 post_download_command=None):
		if destination_dir == None :
			destination_dir = util.file.get_build_path()
		self.url = url
		self.targetpath = os.path.join(destination_dir, target_rel_path)
		self.destination_dir = destination_dir
		self.verifycmd = verifycmd
		self.verifycode = verifycode
		self.attempted = False
		self.installed = False
		self.require_executability = require_executability
		self.post_download_command = post_download_command
	def is_attempted(self):
		return self.attempted
	def is_installed(self):
		return self.installed
	def executable_path(self):
		return self.installed and self.targetpath or None
	def verify_install(self):
		if os.access(self.targetpath,
					 (os.X_OK | os.R_OK) if self.require_executability else os.R_OK):
			if self.verifycmd:
				log.debug("validating")
				self.installed = (os.system(self.verifycmd) == self.verifycode)
			else:	
				self.installed = True
		else:
			self.installed = False
		return self.installed
	def attempt_install(self):
		self.attempted = True
		if not self.verify_install():
			self.pre_download()
			self.download()
			self.post_download()
			self.verify_install()
	def pre_download(self):
		pass
	def download(self):
		log.debug("downloading")
		download_dir = tempfile.gettempdir()
		util.file.mkdir_p(download_dir)
		filepath = urlparse(self.url).path
		filename = filepath.split('/')[-1]
		urlretrieve(self.url, "%s/%s" % (download_dir,filename))
		self.download_file = filename
		self.unpack(download_dir)
	def post_download(self):
		if self.post_download_command:
			os.system('cd "{}"; {}'.format(self.destination_dir, self.post_download_command))
	def unpack(self, download_dir):
		log.debug("unpacking")
		util.file.mkdir_p(self.destination_dir)
		if self.download_file.endswith('.zip'):
			if os.system("unzip -o %s/%s -d %s > /dev/null" % (download_dir, self.download_file,
															   self.destination_dir)):
				return
			else:
				os.unlink("%s/%s" % (download_dir, self.download_file))
		elif (self.download_file.endswith('.tar.gz') or
			  self.download_file.endswith('.tar.bz2') or
			  self.download_file.endswith('.tgz')):
			compression_option = 'j' if self.download_file.endswith('.tar.bz2') else 'z'
			if os.system("tar -C %s -x%spf %s/%s" % (self.destination_dir,
													 compression_option,
													 download_dir,
													 self.download_file)):
				return
			else:
				os.unlink("%s/%s" % (download_dir, self.download_file))
		else :
			os.rename(os.path.join(download_dir, self.download_file),
					  os.path.join(self.destination_dir, self.download_file))
