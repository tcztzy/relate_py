import os
import pathlib
from dataclasses import dataclass
from functools import cached_property

import anndata as ad
import dask.dataframe as dd
import numpy as np
import pandas as pd

from relatepy.utils import logger

lower_bound = 1e-10


def pair(p: str | pd.Series) -> str | pd.Series:
    match p:
        case "A":
            return "G"
        case "C":
            return "T"
        case "G":
            return "A"
        case "T":
            return "C"
        case s if isinstance(s, pd.Series):
            return s.replace({"A": "G", "C": "T", "G": "A", "T": "A"})
        case _:
            raise ValueError


def is_paired(a: str | pd.Series, b: str | pd.Series) -> bool | pd.Series:
    return pair(a) == b


def pack_props(row: pd.Series):
    return b"".join(
        (
            np.array([row["bp_pos"], row["dist"]], dtype="u1").tobytes(),
            row["ID"].encode().ljust(1024, b"\0"),
            row["ancestral"].encode().ljust(1024, b"\0"),
            row["alternative"].encode().ljust(1024, b"\0"),
        )
    )


class HapsFile:
    """Oxford phased haplotype file"""

    section_boundaries = ()
    window_boundaries = []

    def __init__(
        self,
        haps_path: os.PathLike,
        sample_path: os.PathLike,
        dist_path: os.PathLike | None = None,
        use_transition: bool = True,
    ) -> None:
        """
        Parameters
        ----------
        self : HapsFile
        haps_path : os.PathLike
        sample_path : os.PathLike
        """
        haps_path = pathlib.Path(haps_path)
        sample = read_sample(sample_path)
        df: dd.DataFrame = dd.read_csv(
            haps_path,
            sep=r"\s+",
            assume_missing=True,
            na_values=".",
            blocksize=None,
            names=("CHR", "ID", "bp_pos", "ancestral", "alternative", *sample.ids),
        )
        adata = ad.AnnData(
            df[sample.ids].values.T,
            dtype="u1",
            var=df.drop(columns=sample.ids).compute(),
        )
        adata.var["bp_pos"] = adata.var["bp_pos"].astype("u4")
        adata.var["ID"] = adata.var["ID"].astype(str)
        self.data = adata
        self._update_dist(dist_path)
        self.rpos = np.zeros(self.L + 1)
        self.use_transitions = use_transition

    @cached_property
    def N(self):
        """We wouldn't have 4.29 billion samples because lack of money."""
        return np.uint32(self.data.n_obs)

    @cached_property
    def L(self):
        """Human genome is 6.27~6.37 Gbp (male/female), uint32 upper limit is 4GB,
        we do not need to load all the pairs at one time, enough for most cases."""
        return np.uint32(self.data.n_vars)

    @cached_property
    def bp_pos(self):
        return np.append(self.data.var["bp_pos"], self.data.var["bp_pos"][-1:] + 1)

    @property
    def r(self):
        if "recombination_distance" not in self.data.var:
            self.data.var["recombination_distance"] = np.zeros(self.L)
        return self.data.var["recombination_distance"]

    @r.setter
    def r(self, value):
        self.data.var["recombination_distance"] = value

    @property
    def dist(self):
        return self.data.var["dist"]

    @dist.setter
    def dist(self, value):
        self.data.var["dist"] = value

    @cached_property
    def chunks(self):
        return [
            DataChunk(self, i, slice(*boundaries), self.use_transitions)
            for i, boundaries in enumerate(self.section_boundaries)
        ]

    def _update_dist(self, dist_path: os.PathLike | None):
        if dist_path is None:
            dist = np.append(np.diff(self.data.var["bp_pos"]), np.uint32(1))
        else:
            dist = pd.read_csv(dist_path, sep=r"\s+", skiprows=1, names=("bp", "dist"))[
                "dist"
            ].astype("u4")
        if (dist <= 0).any():
            raise ValueError(
                "SNPs are not sorted by bp or more than one SNP at same position."
            )
        self.data.var["dist"] = dist

    def make_chunks(
        self,
        file_out: pathlib.Path,
        filename_map: pathlib.Path,
        filename_dist: pathlib.Path | None = None,
        use_transitions: bool = True,
        min_memory: float = 5.0,
    ):
        if filename_dist is not None:
            self._update_dist(dist_path=filename_dist)
        self.use_transitions = use_transitions
        min_memory_size = min_memory * 1e9 / 4.0 - (2 * self.N**2 + 3 * self.N)
        actual_min_memory_size = 0.0
        if min_memory_size <= 0:
            raise MemoryError("Need larger memory allowance.")
        # this is the number of open files, needs to be less than 500
        windows_per_section = 500
        max_windows_per_section = 0
        overlap = 20000
        max_chunk_size = min(self.L + 1, min_memory_size // self.N)
        if min_memory >= 100:
            max_chunk_size = 2500000

        snp = 0
        window_boundaries = np.zeros(windows_per_section + 1, dtype=np.uint32)
        window_boundaries_overlap = np.zeros(windows_per_section + 1, dtype=np.uint32)
        section_boundary_start = [0]
        section_boundary_end = []
        num_windows: int = 0
        num_windows_overlap: int = 0
        chunk_size: int = 0
        chunk_index: int = 0
        while snp < self.L:
            # output chunk bed into fp_haps_chunk and chunk pos, rpos into fp_pos

            if snp > 0:
                # copy the last #overlap snps;
                assert snp - section_boundary_start[-1] >= overlap
                assert overlap <= chunk_size
                snp_section_begin: int = snp - overlap
                section_boundary_start.append(snp_section_begin)

                window_boundaries_overlap[0] = snp_section_begin
                _window_boundaries = window_boundaries[:num_windows]
                num_windows_overlap = len(_window_boundaries) + 1
                window_boundaries_overlap[1:num_windows_overlap] = _window_boundaries[
                    _window_boundaries > snp_section_begin
                ]
                assert num_windows_overlap < windows_per_section - 1

            snp_begin: int = snp
            window_memory_size = 0.0
            chunk_size = 0

            window_boundaries[0] = snp_begin
            num_windows = 1
            snps_in_window: int = 0
            it_p: int = 0

            while (
                num_windows + num_windows_overlap < windows_per_section
                and chunk_size < max_chunk_size
                and snp < self.L
            ):
                num_derived: int = self.data.X[:, it_p].sum()

                window_memory_size += num_derived * (self.N + 1)
                # 73 comes from  ((N+1+2*Node)*x + N^2 + 3*N) which
                # I am using as an approximation of memory usage
                # I am also assuming one tree in 100 SNPs on average
                if window_memory_size >= min_memory_size and snps_in_window > 10:
                    if actual_min_memory_size < window_memory_size:
                        actual_min_memory_size = window_memory_size
                    snps_in_window = 0
                    window_memory_size = 0.0
                    window_boundaries[num_windows] = snp
                    num_windows += 1

                it_p += 1
                snp += 1
                snps_in_window += 1
                chunk_size += 1
            if actual_min_memory_size < window_memory_size:
                actual_min_memory_size = window_memory_size
            mean_snps_in_window = chunk_size / num_windows
            window_boundaries[num_windows] = snp
            self.window_boundaries.append(window_boundaries[: num_windows + 1].copy())
            assert num_windows <= windows_per_section

            if max_windows_per_section < num_windows:
                max_windows_per_section = num_windows

            if mean_snps_in_window < 100:
                logger.warn(
                    f"Memory allowance should be set {100 / mean_snps_in_window} times "
                    "larger than the current setting using --memory (Default 5GB)."
                )

            section_boundary_end.append(snp)

            num_windows_in_section = num_windows + 1
            if snp_begin == 0:
                (file_out / f"parameters_c{chunk_index}.bin").write_bytes(
                    np.array(
                        [
                            self.N,
                            chunk_size,
                            num_windows_in_section,
                            *self.window_boundaries[0],
                        ],
                        dtype="u4",
                    ).tobytes()
                )
            else:
                L_chunk: int = chunk_size + overlap

                window_start: int = window_boundaries_overlap[0]
                window_boundaries_overlap[:num_windows_overlap] -= window_start

                # Output program parameters into file
                (file_out / f"parameters_c{chunk_index}.bin").write_bytes(
                    np.array(
                        [
                            self.N,
                            L_chunk,
                            num_windows_in_section,
                            *window_boundaries_overlap[:num_windows_overlap],
                            *self.window_boundaries[chunk_index] - window_start,
                        ],
                        dtype="u4",
                    ).tobytes()
                )
            chunk_index += 1

        assert section_boundary_end[-1] == self.L
        assert len(section_boundary_start) == len(section_boundary_end)
        num_chunks: int = len(section_boundary_start)
        self.section_boundaries = tuple(
            zip(section_boundary_start, section_boundary_end)
        )

        logger.warning(
            "Warning: Will use min "
            f"{2.0 * (4.0 * self.N ** 2 * (max_windows_per_section + 2.0)) / 1e9}"
            "GB of hard disc."
        )
        actual_min_memory_size += 2 * self.N**2 + 3 * self.N
        actual_min_memory_size *= 4.0 / 1e9  # to GB
        (file_out / "parameters.bin").write_bytes(
            np.array(
                [
                    self.N,
                    self.L,
                    num_chunks,
                    actual_min_memory_size,
                    *section_boundary_start,
                    *section_boundary_end,
                ],
                dtype=np.uint32,
            ).tobytes()
        )

        m_map = GeneticMapFile(filename_map)

        map_index = np.searchsorted(m_map.bp, self.bp_pos, side="right")
        # Telomere is the
        telomere = (map_index == 0) | (map_index == len(m_map.bp))
        map_pos = np.clip(map_index - 1, 0, None)
        # fmt: off
        # NOTE: black always mess up, disable it for the statement
        rpos = (
            (self.bp_pos - m_map.bp[map_pos])  # bp distance from previous position
            / np.diff(m_map.bp)[map_pos]       # bp distance
            * np.diff(m_map.gen_pos)[map_pos]  # genetic distance
            + m_map.gen_pos[map_pos]           # previous genetic position
        )
        # fmt: on
        self.rpos = np.where(telomere, m_map.gen_pos[map_pos], rpos) * 1e-2
        self.r = np.clip(np.diff(self.rpos), lower_bound, None) * 2500

        self.dump(file_out)

    def dump(self, output: pathlib.Path):
        props = self.data.var.apply(pack_props, axis=1)
        snp_bytes = props.index.astype(bytes)
        (output / "props.bin").write_bytes(b"".join(snp_bytes + props))
        for chunk in self.chunks:
            chunk.dump(output)


@dataclass
class DataChunk:
    """Will only support dump for backward compatibility"""

    data: HapsFile
    id: int
    boundaries: slice
    use_transitions: bool = True

    def __getattr__(self, attr):
        boundaries = (
            self.boundaries
            if attr != "rpos"
            else slice(self.boundaries.start, self.boundaries.stop + 1)
        )
        attr = getattr(self.data, attr)[boundaries]
        if isinstance(attr, pd.Series):
            attr = attr.values
        return attr

    @property
    def bp(self):
        return self.data.bp_pos[self.boundaries]

    @cached_property
    def size(self) -> int:
        return self.boundaries.stop - self.boundaries.start

    @property
    def state(self):
        if self.use_transitions:
            return np.ones(self.size, dtype="u4")
        else:
            return (
                ~is_paired(
                    self.data.data.var["ancestral"][self.boundaries],
                    self.data.data.var["alternative"][self.boundaries],
                )
            ).astype(int)

    @property
    def hap(self):
        return self.data.data.X[:, self.boundaries]

    def dump(self, output: pathlib.Path):
        stem = output / f"chunk_{self.id}"
        for prop in ("bp", "dist", "rpos", "r", "state"):
            value: np.ndarray = getattr(self, prop)
            content = np.uint32(len(value)).tobytes() + value.tobytes()
            stem.with_suffix("." + prop).write_bytes(content)
        stem.with_suffix(".hap").write_bytes(
            np.array([self.size, self.data.N], dtype="u8").tobytes()
            + (self.hap + 48).astype("u1").tobytes("F")
        )


class GeneticMapFile:
    bp: pd.Series
    gen_pos: pd.Series

    def __init__(self, path: pathlib.Path) -> None:
        m = pd.read_csv(path, sep=r"\s", engine="python")
        self.bp = m.iloc[:, 0]
        assert self.bp.is_monotonic_increasing
        self.gen_pos = m.iloc[:, 2]


class SampleFile:
    """Oxford sample information file

    See [1]_, [2]_ and [3]_

    Notes
    -----
    Only provide minimal support for run the model right now

    .. [1] https://www.cog-genomics.org/plink/2.0/formats#sample
    .. [2] https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html
    .. [3] https://myersgroup.github.io/relate/input_data.html
    """

    _data: pd.DataFrame

    def __init__(self, path: os.PathLike) -> None:
        self._data = pd.read_csv(path, sep=r"\s+", skiprows=[1], na_values="NA")

    @cached_property
    def ids(self):
        results = []
        for _, row in self._data.iterrows():
            for i in row[:2]:
                if i:
                    if i not in results:
                        results.append(i)
                    else:
                        results.append(i + "(1)")
        return results


def read_sample(sample_path: os.PathLike) -> SampleFile:
    """Read Oxford sample information file

    Notes
    -----
    This function only provide minimal support for run the model right now

    Parameters
    ----------
    sample_path : os.PathLike

    Returns
    -------
    SampleFile
    """
    return SampleFile(sample_path)


def read_haps(
    haps_path: os.PathLike, sample_path: os.PathLike | None = None
) -> HapsFile:
    """Read Oxford phased haplotype file

    Parameters
    ----------
    haps_path : os.PathLike
    sample_path : os.PathLike | None, optional
        default None

    Returns
    -------
    HapsFile

    Raises
    ------
    ValueError
    """
    haps_path = pathlib.Path(haps_path)
    if sample_path is None:
        if ".haps" not in haps_path.suffixes:
            raise ValueError(
                "sample file is not provided and haps file doesn't have .haps suffix, "
                "it is impossible to guess the sample file path."
            )
        sample_path = haps_path.parent / haps_path.name.replace(".haps", ".sample")
    return HapsFile(haps_path, sample_path)


def read_coal(filename: os.PathLike) -> pd.DataFrame:
    with open(filename) as f:
        groups = f.readline().split()
        epochs = [float(epoch) for epoch in f.readline().split()]
        coal = pd.read_table(
            f, sep=r"\s+", na_values="nan", names=["group1", "group2", *epochs]
        )
    if (
        len(groups) >= 1 + coal["group1"].max()
        and len(groups) >= 1 + coal["group2"].max()
    ):
        coal["group1"] = [groups[i] for i in coal["group1"]]
        coal["group2"] = [groups[i] for i in coal["group2"]]
    return coal.melt(
        id_vars=["group1", "group2"],
        var_name="epoch.start",
        value_name="haploid.coalescence.rate",
    )
