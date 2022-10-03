import logging
import platform
from functools import wraps
from resource import RUSAGE_SELF, getrusage
import os
from typing import Callable

import click_log

logger = logging.getLogger(__package__)
click_log.basic_config(logger)


def resource_usage(func: Callable) -> Callable:
    @wraps(func)
    def wrapper(*args, **kwargs):
        func(*args, **kwargs)
        usage = getrusage(RUSAGE_SELF)
        logger.debug(
            f"CPU Time spent: {usage.ru_utime:.6}s; "
            f"Max memory usage: {usage.ru_maxrss/(1000_000 if platform.system() == 'Darwin' else 1000)}Mb."
        )

    return func if os.getenv("RELATEPY_RESOURCE_USAGE") is None else wrapper
