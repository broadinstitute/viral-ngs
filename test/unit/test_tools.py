# Unit tests for tools/__init__.py

__author__ = "yesimon@broadinstitute.org"

import tools
from tools import *
import pytest


@pytest.fixture(params=tools.all_tool_classes())
def tool_class(request):
    print(request.param)
    return request.param


def test_tool_install(tool_class):
    t = tool_class()
    t.install()
    assert t.is_installed()
