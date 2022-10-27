import os
import pathlib
from dataclasses import dataclass
from functools import cached_property
from struct import pack
from typing import Sequence

import anndata as ad
import dask.dataframe as dd
import numpy as np
import pandas as pd

from relatepy.utils import logger

lower_bound = 1e-10


@dataclass
class Alias:
    source_name: str

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return getattr(obj, self.source_name)


class HapsFile:
    """Oxford phased haplotype file"""

    def __init__(self, haps_path: os.PathLike, sample_path: os.PathLike) -> None:
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
        # adata.obs_names = sample.ids
        adata.var["bp_pos"] = adata.var["bp_pos"].astype(int)
        adata.var["dist"] = np.append(np.diff(adata.var["bp_pos"]), [1])
        self._data = adata
        self.rpos = np.zeros(self.L + 1)

    def __getattr__(self, name: str):
        return getattr(self._data, name)

    N = Alias("n_obs")
    L = Alias("n_vars")

    @cached_property
    def bp_pos(self):
        return np.append(self._data.var["bp_pos"], self._data.var["bp_pos"][-1:] + 1)

    @property
    def r(self):
        if "recombination_distance" not in self._data.var:
            self._data.var["recombination_distance"] = np.zeros(self.L)
        return self._data.var["recombination_distance"]

    @r.setter
    def r(self, value):
        self._data.var["recombination_distance"] = value

    def make_chunks(
        self,
        file_out: pathlib.Path,
        filename_map: pathlib.Path,
        filename_dist: pathlib.Path | None = None,
        use_transitions: bool = True,
        min_memory: float = 5.0,
    ):
        bp_pos = self.bp_pos
        ancestral = self._data.var["ancestral"]
        alternative = self._data.var["alternative"]
        rsid = self._data.var["ID"]
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

        p_seq: np.ndarray = np.zeros((max_chunk_size, self.N), dtype="u1")
        p_overlap: np.ndarray = np.zeros((overlap, self.N), dtype="u1")

        snp = 0
        window_boundaries = np.zeros(windows_per_section + 1, dtype=np.int32)
        window_boundaries_overlap = np.zeros(windows_per_section + 1)
        section_boundary_start = [0]
        section_boundary_end = []
        state_val: int = 1
        min_snps_in_window: int = max_chunk_size
        mean_snps_in_window: float = 0.0
        num_windows: int
        num_windows_overlap: int = 0
        overlap_in_section: int
        chunk_size: int
        chunk_index: int = 0
        window_memory_size = 0.0
        while snp < self.L:
            fp_haps_chunk = open(file_out / f"chunk_{chunk_index}.hap", "wb")
            fp_state = open(file_out / f"chunk_{chunk_index}.state", "wb")
            # output chunk bed into fp_haps_chunk and chunk pos, rpos into fp_pos

            if snp > 0:
                # copy the last #overlap snps;
                assert snp - section_boundary_start[-1] >= overlap
                overlap_in_section = overlap
                assert overlap_in_section <= chunk_size
                snp_section_begin: int = snp - overlap_in_section
                section_boundary_start.append(snp_section_begin)
                p_overlap[:overlap_in_section] = p_seq[
                    chunk_size - overlap_in_section : chunk_size
                ].copy()

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
                p_seq[it_p] = self._data.X[:, it_p]
                num_derived: int = p_seq[it_p].sum()

                window_memory_size += num_derived * (
                    self.N + 1
                )  # + 2.88*N; //2.88 = 72/25 (72 bytes per 2 nodes, 1 tree in 25 snps)
                # 73 comes from  ((N+1+2*Node)*x + N^2 + 3*N) which I am using as an approximation of memory usage
                # I am also assuming one tree in 100 SNPs on average
                if window_memory_size >= min_memory_size and snps_in_window > 10:
                    if actual_min_memory_size < window_memory_size:
                        actual_min_memory_size = window_memory_size
                    if min_snps_in_window > snps_in_window:
                        min_snps_in_window = snps_in_window
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
            if min_snps_in_window > snps_in_window:
                min_snps_in_window = snps_in_window
            mean_snps_in_window = chunk_size / num_windows
            window_boundaries[num_windows] = snp
            assert num_windows <= windows_per_section

            if num_windows > max_windows_per_section:
                max_windows_per_section = num_windows

            if mean_snps_in_window < 100:
                logger.warn(
                    f"Memory allowance should be set {100/mean_snps_in_window} times larger than the current setting using --memory (Default 5GB)."
                )

            section_boundary_end.append(snp)

            snp_tmp = section_boundary_start[-1]
            if snp_begin == 0:
                fp_haps_chunk.write(pack("QQ", chunk_size, self.N))
                num_windows_in_section = num_windows + 1
                (file_out / f"parameters_c{chunk_index}.bin").write_bytes(
                    pack(
                        "III" + "I" * num_windows_in_section,
                        self.N,
                        chunk_size,
                        num_windows_in_section,
                        *window_boundaries[:num_windows_in_section],
                    )
                )
                fp_state.write(pack("I", chunk_size))
            else:
                L_chunk: int = chunk_size + overlap_in_section
                fp_haps_chunk.write(pack("QQ", L_chunk, self.N))

                window_start: int = window_boundaries_overlap[0]
                window_boundaries_overlap[:num_windows_overlap] -= window_start

                # Output program parameters into file
                (file_out / f"parameters_c{chunk_index}.bin").write_bytes(
                    pack(
                        "III" + "I" * (num_windows + 1 + num_windows_overlap),
                        self._data.n_obs,
                        L_chunk,
                        num_windows_in_section,
                        *window_boundaries_overlap[:num_windows_overlap],
                        *(window_boundaries[: num_windows + 1] - window_start),
                    )
                )

                fp_state.write(pack("I", L_chunk))
                for p in p_overlap:
                    if not use_transitions:
                        if (
                            (ancestral[snp_tmp] == "C" and alternative[snp_tmp] == "T")
                            or (
                                ancestral[snp_tmp] == "T"
                                and alternative[snp_tmp] == "C"
                            )
                            or (
                                ancestral[snp_tmp] == "A"
                                and alternative[snp_tmp] == "G"
                            )
                            or (
                                ancestral[snp_tmp] == "G"
                                and alternative[snp_tmp] == "A"
                            )
                        ):
                            state_val = 0
                        else:
                            state_val = 1
                    fp_state.write(pack("I", state_val))
                    snp_tmp += 1
                    fp_haps_chunk.write(bytes(ord(str(i)) for i in p))
            for p in p_seq[:chunk_size]:
                if not use_transitions:
                    if (
                        (ancestral[snp_tmp] == "C" and alternative[snp_tmp] == "T")
                        or (ancestral[snp_tmp] == "T" and alternative[snp_tmp] == "C")
                        or (ancestral[snp_tmp] == "A" and alternative[snp_tmp] == "G")
                        or (ancestral[snp_tmp] == "G" and alternative[snp_tmp] == "A")
                    ):
                        state_val = 0
                    else:
                        state_val = 1
                fp_state.write(pack("I", state_val))
                snp_tmp += 1
                fp_haps_chunk.write(bytes(ord(str(i)) for i in p))
            fp_haps_chunk.close()
            fp_state.close()
            chunk_index += 1

        assert section_boundary_end[-1] == self.L
        assert len(section_boundary_start) == len(section_boundary_end)
        num_chunks: int = len(section_boundary_start)

        logger.warning(
            f"Warning: Will use min {2.0 * (4.0 * self.N ** 2 * (max_windows_per_section + 2.0)) / 1e9}GB of hard disc."
        )
        fp = open(file_out / "parameters.bin", "wb")
        actual_min_memory_size += 2 * self.N**2 + 3 * self.N
        actual_min_memory_size *= 4.0 / 1e9
        fp.write(
            pack(
                "IIId" + "I" * 2 * num_chunks,
                self.N,
                self.L,
                num_chunks,
                actual_min_memory_size,
                *section_boundary_start,
                *section_boundary_end,
            )
        )
        fp.close()
        # calculate rpos, r, pos here
        if filename_dist is None:
            dist = self._data.var["dist"]
            if (dist <= 0).any():
                raise ValueError(
                    "SNPs are not sorted by bp or more than one SNP at same position."
                )
        else:
            df = pd.read_csv(
                filename_dist, sep=r"\s+", skiprows=1, names=("mbp", "mdist")
            )
            dist = df["mdist"].values
            self._data.var["dist"] = dist
        with open(file_out / "props.bin", "wb") as fp_props:
            for snp in range(self.L):
                fp_props.write(
                    pack(
                        "III1024s1024s1024s",
                        snp,
                        int(bp_pos[snp]),
                        int(dist[snp]),
                        f"{rsid[snp]:\0<1024}".encode(),
                        f"{ancestral[snp]:\0<1024}".encode(),
                        f"{alternative[snp]:\0<1024}".encode(),
                    )
                )

        m_map = GeneticMapFile(filename_map)

        map_index = np.searchsorted(m_map.bp, bp_pos, side="right")

        if map_index[-1] < len(m_map.bp) - 2:
            # while there are still at last one position
            last_pos = map_index[-1] + 1
        else:
            last_pos = map_index[-1]
        map_diff_index = np.append(map_index, last_pos)
        map_diff = np.diff(m_map.bp[map_diff_index])
        rpos = (bp_pos - m_map.bp[map_index]) / map_diff * np.diff(
            m_map.gen_pos[map_diff_index]
        ) + m_map.gen_pos[map_index]
        cond1 = (map_diff == 0) | (bp_pos < m_map.bp[map_index])
        self.rpos = np.select([cond1, ~cond1], [m_map.gen_pos[map_index], rpos]) * 1e-2
        self.r = np.clip(np.diff(self.rpos), lower_bound, np.inf) * 2500

        for chunk in range(num_chunks):
            L_chunk = section_boundary_end[chunk] - section_boundary_start[chunk]
            L_chunk_plus_one = L_chunk + 1
            (file_out / f"chunk_{chunk}.bp").write_bytes(
                pack(
                    "I" + "i" * L_chunk,
                    L_chunk,
                    *bp_pos[
                        section_boundary_start[chunk] : section_boundary_end[chunk]
                    ],
                )
            )
            (file_out / f"chunk_{chunk}.dist").write_bytes(
                pack(
                    "I" + "i" * L_chunk,
                    L_chunk,
                    *dist[section_boundary_start[chunk] : section_boundary_end[chunk]],
                )
            )
            (file_out / f"chunk_{chunk}.rpos").write_bytes(
                pack(
                    "=I" + "d" * L_chunk_plus_one,
                    L_chunk_plus_one,
                    *self.rpos[
                        section_boundary_start[chunk] : section_boundary_end[chunk] + 1
                    ],
                )
            )
            (file_out / f"chunk_{chunk}.r").write_bytes(
                pack(
                    "I" + "d" * L_chunk,
                    L_chunk,
                    *self.r[
                        section_boundary_start[chunk] : section_boundary_end[chunk]
                    ],
                )
            )


class GeneticMapFile:
    bp: Sequence[int]
    gen_pos: Sequence[float]

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
