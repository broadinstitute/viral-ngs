# Unit tests for tools/__init__.py

__author__ = "yesimon@broadinstitute.org"

import pytest
import tools
from tools import *

# Simply do nothing to override stub_conda in conftest.py
@pytest.fixture(autouse=True)
def stub_conda():
    pass

@pytest.fixture(params=tools.all_tool_classes())
def tool_class(request):
    return request.param

def test_tool_install(tool_class):
    t = tool_class()
    t.install()
    assert t.is_installed()
