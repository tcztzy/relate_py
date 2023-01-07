import logging
from pathlib import Path

import click_log

from ..utils import resource_usage
from ..io import HapsFile
from ..relatepy import make_chunks

logger = logging.getLogger(__package__)
click_log.basic_config(logger)


@resource_usage
def chunk(
    haps: Path,
    sample: Path,
    genetic_map: Path,
    output: Path,
    dist: Path | None = None,
    use_transitions: bool = True,
    memory_limit: float = 5.0,
) -> None:
    logger.debug("Parsing data.")
    make_chunks(haps, sample, genetic_map, output, dist, use_transitions, memory_limit)
