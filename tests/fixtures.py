from pathlib import Path

import pytest


@pytest.fixture
def haps_path():
    return (
        Path(__file__).parent.parent / "relate" / "example" / "data" / "example.haps.gz"
    )


@pytest.fixture
def sample_path():
    return (
        Path(__file__).parent.parent
        / "relate"
        / "example"
        / "data"
        / "example.sample.gz"
    )
