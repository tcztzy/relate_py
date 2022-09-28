# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython

from cython.cimports.libc.stdlib import free, malloc  # type: ignore
from cython.cimports.relatepy.data import BuildTopology, CombineSections, Data, Finalize, FindEquivalentBranches, GetBranchLengths, Options, Paint, get_options  # type: ignore


@cython.cclass
class RelateData:
    data = cython.declare(Data)

    def __init__(
        self,
        *,
        number_of_sequences: int | None = None,
        number_of_alleles: int | None = None,
        hap: Path | None = None,
        pos: Path | None = None,
        dist: Path | None = None,
        rec: Path | None = None,
        rpos: Path | None = None,
        state: Path | None = None,
        param: Path | None = None,
        Ne: int = 30000,
        mu: float = 1.25e-8,
    ):
        if number_of_sequences is not None and number_of_alleles is not None:
            self.data = Data()
        elif dist is not None and param is not None:
            self.data = Data(bytes(dist), bytes(param), Ne, mu)
        elif all((hap, pos, dist, rec, rpos, state)):
            self.data = Data(
                bytes(hap),
                bytes(pos),
                bytes(dist),
                bytes(rec),
                bytes(rpos),
                bytes(state),
                Ne,
                mu,
            )
        else:
            raise ValueError

    @staticmethod
    def chunk(
        haps: Path,
        sample: Path,
        genetic_map: Path,
        dist: Path | None,
        output: Path,
        use_transitions: bool,
        memory_limit: float = 5.0,
    ):
        Data().MakeChunks(
            bytes(haps),
            bytes(sample),
            bytes(genetic_map),
            b"unspecified" if dist is None else bytes(dist),
            bytes(output),
            use_transitions,
            memory_limit,
        )

    @property
    def name(self):
        return self.data.name

    @name.setter
    def name(self, value: str):
        self.data.name = value.encode("utf-8")

    @property
    def number_of_sequences(self) -> int:
        return self.data.N

    @number_of_sequences.setter
    def number_of_sequences(self, value: int):
        self.data.N = value

    @property
    def number_of_alleles(self) -> int:
        return self.data.L

    @number_of_alleles.setter
    def number_of_alleles(self, value: int):
        self.data.L = value

    @property
    def effective_population_size(self):
        return self.data.Ne

    @effective_population_size.setter
    def effective_population_size(self, value):
        self.data.Ne = value

    @property
    def theta(self) -> float:
        return self.data.theta

    @theta.setter
    def theta(self, value: float):
        self.data.theta = value

    @property
    def ntheta(self) -> float:
        return 1 - self.theta

    @property
    def r(self) -> list:
        return self.data.r

    @r.setter
    def r(self, value: list):
        self.data.r = value


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


def build_topology(
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
    args = [b"relate", b"--output", bytes(output)]
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


def find_equivalent_branches(output: Path, chunk_index: int):
    args = [b"relate", b"--output", bytes(output)]
    argv: cython.pp_char = cython.cast(
        cython.pp_char, malloc(cython.sizeof(cython.p_char) * len(args))
    )
    for i, arg in enumerate(args):
        argv[i] = arg
    options: Options = get_options(len(args), argv)
    FindEquivalentBranches(options, chunk_index)
    free(argv)


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
