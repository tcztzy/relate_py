# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import Options, get_options, Paint  # type: ignore

from ..utils import output_working_directory


@output_working_directory
def paint(
    *, output: Path, chunk_index: int, theta: float | None = None, rho: float | None = None
):
    args = [b"relate", b"--output", output.name.encode()]
    if theta is not None and rho is not None:
        args += [b"--painting", f"{theta},{rho}".encode()]
    elif (theta is not None and rho is None) or (theta is None and rho is not None):
        raise ValueError
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, v in enumerate(args):
        argv[i] = v
    options: Options = get_options(len(args), argv)
    Paint(options, chunk_index)
    free(argv)
