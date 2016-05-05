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
import pytest

log = logging.getLogger(__name__)
util.cmd.setup_logger('INFO')


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


def all_tool_tests():
    for tool_class in iter_leaf_subclasses(tools.Tool):
        yield tool_class


@pytest.fixture(params=all_tool_tests())
def tool_class(request):
    print(request.param)
    return request.param


def test_tool_install(tool_class):
    t = tool_class()
    t.install()
    assert t.is_installed(), "installation of tool %s failed" % tool_class.__name__
    log.info(".. installation of %s succeeded with installer %s" %
             (tool_class.__name__, t.installed_method.__class__.__name__))
