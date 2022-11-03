import logging
from pathlib import Path

import click_log

from ..utils import resource_usage
from ..io import HapsFile

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
    HapsFile(haps, sample, dist, use_transitions).make_chunks(
        output, genetic_map, min_memory=memory_limit
    )
