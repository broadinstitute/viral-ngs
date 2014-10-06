'''Stuff here'''

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import logging, re, os
import util.file, util.misc

log = logging.getLogger(__name__)


class Tool:
	''' 
	'''
	def __init__(self):
		pass
	def is_installed(self):
		pass
	def install(self):
		pass
	def get_install_methods(self):
		pass
	def set_install_methods(self, methods):
		pass
	def executable_path(self):
		pass
	def execute(self, args):
		pass

class InstallMethod:
	'''
	'''
	def __init__(self):
		pass
	def is_attempted(self):
		pass
	def attempt_install(self):
		pass
	def is_installed(self):
		pass
	def executable_path(self):
		pass


