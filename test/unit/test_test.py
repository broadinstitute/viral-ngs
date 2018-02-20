"""Test the test helper code"""

from test import assert_equal_flat_dicts

import pytest

def test_cmp_flat_dicts():
    """Test comparison of flattened dicts"""

    expect_eq = (
        ({1:2}, {1:3}, [(1,)]),
        ({1:{2:3}}, {1:{2:4}}, [(1,2)])
    )
    expect_ne = (
        ({1:2}, {1:3}),
        ({1:{2:3}}, {1:{2:4}}),
        ({1:{2:3}}, {1:{2:'4'}}),
        ({1:{2:3}}, {1:{0:3}}),
    )

    for args in expect_eq:
        assert_equal_flat_dicts(*args)

    for args in expect_ne:
        with pytest.raises(AssertionError):
            assert_equal_flat_dicts(*args)
            
