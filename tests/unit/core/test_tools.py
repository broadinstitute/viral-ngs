# Unit tests for tool installation

__author__ = "yesimon@broadinstitute.org"

import platform
import pytest
import viral_ngs.core

# Platform detection for x86-only tools
IS_ARM = platform.machine() in ('arm64', 'aarch64')

# Tool class names that require x86-only bioconda packages
X86_ONLY_TOOLS = {'NovoalignTool', 'MvicunaTool'}


# Simply do nothing to override stub_conda in conftest.py
@pytest.fixture(autouse=True)
def stub_conda():
    pass


@pytest.fixture(params=viral_ngs.core.all_tool_classes())
def tool_class(request):
    return request.param


def test_tool_install(tool_class):
    if IS_ARM and tool_class.__name__ in X86_ONLY_TOOLS:
        pytest.skip(f"{tool_class.__name__} requires x86-only bioconda package (not available on ARM)")
    t = tool_class()
    t.install()
    assert t.is_installed()
