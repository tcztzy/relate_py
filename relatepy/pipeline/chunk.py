import logging
from pathlib import Path

import click_log

from ..data import RelateData
from ..utils import resource_usage

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
    output.mkdir()
    RelateData.chunk(
        haps, sample, genetic_map, dist, output, use_transitions, memory_limit
    )
