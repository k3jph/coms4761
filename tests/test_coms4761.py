"""Tests for `coms4761` module."""
from typing import Generator

import pytest

import coms4761


@pytest.fixture
def version() -> Generator[str, None, None]:
    """Sample pytest fixture."""
    yield coms4761.__version__


def test_version(version: str) -> None:
    """Sample pytest test function with the pytest fixture as an argument."""
    assert version == "0.1.0"
