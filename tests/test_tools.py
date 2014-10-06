# Unit tests for tools/__init__.py

__author__ = "dpark@broadinstitute.org"

import tools
from tools import *
import unittest, tempfile, shutil, os
import util.cmd

def set_tmpDir(name):
	proposed_prefix = ['tmp']
	if name:
		proposed_prefix.append(name)
	for e in ('LSB_JOBID','LSB_JOBINDEX'):
		if e in os.environ:
			proposed_prefix.append(os.environ[e])
	tempfile.tempdir = tempfile.mkdtemp(prefix='-'.join(proposed_prefix)+'-',
		dir=util.cmd.find_tmpDir())
def destroy_tmpDir():
	shutil.rmtree(tempfile.tempdir)

class TestToolsInstallation(unittest.TestCase):
	def setUp(self):
		set_tmpDir('TestToolsInstallation')
	def tearDown(self):
		destroy_tmpDir()
	def testAllToolInstallers(self):
		'''Load every tool's default chain of install methods and try them.'''
		for tool_class in tools.Tool.__subclasses__():
			t = tool_class()
			t.install()
			self.assertTrue(t.is_installed(), "installation of tool %s failed" % tool_class.__name__)
			
