# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import FindEquivalentBranches, Options, get_options  # type: ignore


def find_equivalent_branches(output: Path, chunk_index: int):
    args = [b"relate", b"--output", output.name.encode()]
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, arg in enumerate(args):
        argv[i] = arg
    options: Options = get_options(len(args), argv)
    FindEquivalentBranches(options, chunk_index)
    free(argv)
