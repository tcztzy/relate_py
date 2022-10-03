# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import Options, get_options, Paint  # type: ignore


def paint(
    output: Path, chunk_id: int, theta: float | None = None, rho: float | None = None
):
    args = [b"relate", b"--output", bytes(output)]
    if theta:
        args += [b"--painting", f"{theta},{rho}".encode()]
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, v in enumerate(args):
        argv[i] = v
    options: Options = get_options(len(args), argv)
    Paint(options, chunk_id)
    free(argv)
