# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import CombineSections, Options, get_options  # type: ignore


def combine_sections(
    output: Path, chunk_index: int, effective_population_size: float | None = None
):
    args = [b"relate", b"--output", bytes(output)]
    if effective_population_size is not None:
        args.extend([b"--effectiveN", f"{effective_population_size}".encode()])
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, arg in enumerate(args):
        argv[i] = arg
    options: Options = get_options(len(args), argv)
    CombineSections(options, chunk_index)
    free(argv)
