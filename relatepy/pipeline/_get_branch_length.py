# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import GetBranchLengths, Options, get_options  # type: ignore


def get_branch_length(
    output: Path,
    chunk_index: int,
    first_section: int,
    last_section: int,
    mutation_rate: float,
    effective_population_size: float,
    coal: Path | None = None,
    seed: int | None = None,
    sample_ages: Path | None = None,
):
    args = [
        b"relate",
        b"--output",
        bytes(output),
        b"--mutation_rate",
        f"{mutation_rate}".encode(),
    ]
    args.extend(
        [b"--effectiveN", f"{effective_population_size}".encode()]
        if coal is None
        else [b"--coal", bytes(coal)]
    )
    if seed is not None:
        args.extend([b"--seed", f"{seed}".encode()])
    if sample_ages is not None:
        args.extend([b"--sample_ages", bytes(sample_ages)])
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, arg in enumerate(args):
        argv[i] = arg
    options: Options = get_options(len(args), argv)
    GetBranchLengths(options, chunk_index, first_section, last_section)
    free(argv)
