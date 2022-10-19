# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython
from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.pipeline import BuildTopology, Options, get_options  # type: ignore

from ..utils import output_working_directory


@output_working_directory
def build_topology(
    *,
    output: Path,
    chunk_index: int,
    first_section: int,
    last_section: int,
    effective_population_size: float | None = None,
    theta: float = 0.001,
    rho: float = 1.0,
    seed: int | None = None,
    ancestral_state: bool = True,
    sample_ages: list[float] | None = None,
):
    args = [b"relate", b"--output", output.name.encode()]
    if effective_population_size is not None:
        args += [b"--effectiveN", f"{effective_population_size}".encode()]
    if theta:
        args += [b"--painting", f"{theta},{rho}".encode()]
    if seed is not None:
        args += [b"--seed", f"{seed}".encode()]
    if not ancestral_state:
        args.append(b"--anc_allele_unknown")
    if sample_ages is not None:
        args.extend([b"--sample_ages", bytes(sample_ages)])
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, arg in enumerate(args):
        argv[i] = arg
    options: Options = get_options(len(args), argv)
    BuildTopology(options, chunk_index, first_section, last_section)
    free(argv)
