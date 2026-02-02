# Unit tests for tool installation

__author__ = "yesimon@broadinstitute.org"

import pytest
import viral_ngs.core


# Simply do nothing to override stub_conda in conftest.py
@pytest.fixture(autouse=True)
def stub_conda():
    pass


@pytest.fixture(params=viral_ngs.core.all_tool_classes())
def tool_class(request):
    return request.param


def test_tool_install(tool_class):
    t = tool_class()
    t.install()
    assert t.is_installed()
