import logging
import platform
from pathlib import Path
from contextlib import contextmanager
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
        result = func(*args, **kwargs)
        usage = getrusage(RUSAGE_SELF)
        logger.debug(
            f"CPU Time spent: {usage.ru_utime:.6}s; "
            f"Max memory usage: {usage.ru_maxrss/(1000_000 if platform.system() == 'Darwin' else 1000)}Mb."
        )
        return result

    return func if os.getenv("RELATEPY_RESOURCE_USAGE") is None else wrapper


@contextmanager
def chdir(original_dir, dest_dir):
    try:
        yield os.chdir(dest_dir)
    finally:
        os.chdir(original_dir)


def output_working_directory(func: Callable) -> Callable:
    @wraps(func)
    def wrapper(output: Path, **kwargs):
        if output.is_absolute():
            with chdir(Path.cwd(), output.parent):
                result = func(output=output, **kwargs)
        else:
            result = func(output=output, **kwargs)
        return result

    return wrapper
