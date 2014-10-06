'''Stuff here'''

__author__ = "dpark@broadinstitute.org"

import os

class Tool:
	''' Base tool class that includes install machinery.
	'''
	def __init__(self, install_methods=[]):
		self.install_methods = install_methods
		self.installed_method = None
		self.exec_path = None
	def is_installed(self):
		return (self.installed_method != None)
	def install(self):
		if not is_installed():
			for m in self.install_methods:
				if not m.is_attempted():
					m.attempt_install()
				if m.is_installed():
					self.installed_method = m
					self.exec_path = m.executable_path()
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

class InstallMethod:
	''' Base class for installation methods for a given tool
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
	''' This is a method that tries to find whether an executable binary already exists
	'''
	def __init__(self, path, verifycmd=None, verifycode=0):
		self.path = path
		self.verifycmd = verifycmd
		self.verifycode = verifycode
		self.attempted = False
		self.installed = False
	def is_attempted(self):
		return self.attempted
	def attempt_install(self):
		self.attempted = True
		if os.access(self.path, os.X_OK | os.R_OK):
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


