# cython: language_level=3, cpp_locals=True
from pathlib import Path

import cython

from cython.cimports.relatepy.data import Data  # type: ignore


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
