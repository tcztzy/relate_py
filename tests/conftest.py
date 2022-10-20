from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def haps_path():
    return (
        Path(__file__).parent.parent / "relate" / "example" / "data" / "example.haps.gz"
    )


@pytest.fixture(scope="session")
def sample_path():
    return (
        Path(__file__).parent.parent
        / "relate"
        / "example"
        / "data"
        / "example.sample.gz"
    )


@pytest.fixture(scope="session")
def genetic_map_path():
    return (
        Path(__file__).parent.parent
        / "relate"
        / "example"
        / "data"
        / "genetic_map_GRCh37_chr1.txt"
    )


@pytest.fixture(scope="session")
def sample_ages_path():
    return (
        Path(__file__).parent.parent / "relate" / "example" / "data" / "sample_ages.txt"
    )
