# Unit tests for tools/__init__.py

__author__ = "yesimon@broadinstitute.org"

import pytest
import viral_ngs.core
from viral_ngs.core import Tool, PrexistingUnixCommand

@pytest.fixture(params=viral_ngs.core.all_tool_classes())
def tool_class(request):
    return request.param

@pytest.mark.slow
def test_tool_install(tool_class):
    t = tool_class()
    t.install()
    assert t.is_installed()
