# Unit tests for tool installation

__author__ = "yesimon@broadinstitute.org"

import pytest
import viral_ngs.core
from tests import HAS_NOVOALIGN, SKIP_NO_NOVOALIGN_REASON


# Simply do nothing to override stub_conda in conftest.py
@pytest.fixture(autouse=True)
def stub_conda():
    pass


@pytest.fixture(params=viral_ngs.core.all_tool_classes())
def tool_class(request):
    return request.param


def test_tool_install(tool_class):
    if tool_class.__name__ == 'NovoalignTool' and not HAS_NOVOALIGN:
        pytest.skip(SKIP_NO_NOVOALIGN_REASON)
    t = tool_class()
    t.install()
    assert t.is_installed()
