# Unit tests for tools/__init__.py

__author__ = "dpark@broadinstitute.org"

import tools
from tools import *
import unittest
import tempfile
import shutil
import os
import logging
import util.cmd
import util.file
from test import TestCaseWithTmp
import operator
import nose

log = logging.getLogger(__name__)
util.cmd.setup_logger('INFO')

# nose cannot use test generators with subclasses of unittest.TestCase
# so we must have setUp() and tearDown() defined here
# and setup/destroy the tmp directory as appropriate
class TestToolsInstallation(object):

    def setUp(self):
        #super(TestToolsInstallation, self).setUp()
        util.file.set_tmp_dir(type(self).__name__)

    def tearDown(self):
        util.file.destroy_tmp_dir()

    @nose.with_setup(setUp, tearDown)
    def check_tool(self, tool_class):
        t = tool_class()
        t.install()
        assert t.is_installed(), "installation of tool %s failed" % tool_class.__name__
        log.info(".. installation of %s succeeded with installer %s" %
                 (tool_class.__name__, t.installed_method.__class__.__name__))

    def testAllToolInstallers(self):

        def iter_leaf_subclasses(aClass):
            "Iterate over subclasses at all levels that don't themselves have a subclass"
            isLeaf = True
            for subclass in sorted(aClass.__subclasses__(), key=operator.attrgetter("__name__")):
                isLeaf = False
                for leafClass in iter_leaf_subclasses(subclass):
                    if not getattr(leafClass, '_skiptest', False):
                        yield leafClass
            if isLeaf:
                yield aClass

        '''Load every tool's default chain of install methods and try them.'''
        for tool_class in iter_leaf_subclasses(tools.Tool):
            yield self.check_tool, tool_class



