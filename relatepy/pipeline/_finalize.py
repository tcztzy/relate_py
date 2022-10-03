# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import Finalize, Options, get_options  # type: ignore


def finalize(
    output: Path, sample_ages: Path | None = None, annotation: Path | None = None
):
    args = [b"relate", b"--output", bytes(output)]
    if sample_ages is not None:
        args.extend([b"--sample_ages", bytes(sample_ages)])
    if annotation is not None:
        args.extend([b"--annot", bytes(annotation)])
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, arg in enumerate(args):
        argv[i] = arg
    options: Options = get_options(len(args), argv)
    Finalize(options)
    free(argv)
