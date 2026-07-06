from collections.abc import Callable
from inspect import isfunction

import pytest
from defopt import signature

from divref import main


@pytest.mark.parametrize("tool", main._tools)
def test_tools_are_defined(tool: Callable[..., None]) -> None:
    """Test that all command line tools passed to defopt are defined functions."""
    assert isfunction(tool)


@pytest.mark.parametrize("tool", main._tools)
def test_tools_have_valid_docstrings(tool: Callable[..., None]) -> None:
    """Test that all command line tools have a valid defopt docstring."""
    try:
        signature(tool)
    except TypeError:
        raise AssertionError(f"defopt could not parse docstring for {tool.__name__}") from None
