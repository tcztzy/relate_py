from pathlib import Path
import base64

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


@pytest.fixture(scope="session")
def paint_bin():
    content = base64.b64decode(b'''\
AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAsiMSPrIjEj6yIxI+
siMSPrIjEj6yIxI+siMSPgEAAAAAAAAACAAAAAAAAAAt/wEA8EhCQgAAgD8AAIA/
AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAA
AAAAAAAAAACyIxI+AAAAALIjEj6yIxI+siMSPrIjEj6yIxI+siMSPgEAAAAAAAAA
CAAAAAAAAAAt/wEA05dCQgAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/
AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAAAAAAAAAAAACyIxI+siMSPgAAAACyIxI+
siMSPrIjEj6yIxI+siMSPgEAAAAAAAAACAAAAAAAAAAt/wEAMKxBQgAAgD8AAIA/
AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAA
AAAAAAAAAACyIxI+siMSPrIjEj4AAAAAsiMSPrIjEj6yIxI+siMSPgEAAAAAAAAA
CAAAAAAAAAAt/wEArGdCQgAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/
AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAAAAAAAAAAAACyIxI+siMSPrIjEj6yIxI+
AAAAALIjEj6yIxI+siMSPgEAAAAAAAAACAAAAAAAAAAt/wEADiRCQgAAgD8AAIA/
AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAA
AAAAAAAAAACyIxI+siMSPrIjEj6yIxI+siMSPgAAAACyIxI+siMSPgEAAAAAAAAA
CAAAAAAAAAAt/wEAGfxBQgAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/
AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAAAAAAAAAAAADsyxU57MsVOezLFTnsyxU5
7MsVOezLFTkAAAAA7MsVOQEAAAAAAAAACAAAAAAAAAAt/wEAHtNBQgAAgD8AAIA/
AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/AAAAAC3/AQABAAAAAAAAAAgAAAAAAAAA
AAAAAAAAAACyIxI+siMSPrIjEj6yIxI+siMSPrIjEj6yIxI+AAAAAAEAAAAAAAAA
CAAAAAAAAAAt/wEAUKxCQgAAgD8AAIA/AACAPwAAgD8AAIA/AACAPwAAgD8AAIA/''')
    return content
