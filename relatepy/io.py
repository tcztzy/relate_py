import gzip
import mimetypes
import os
import pathlib
from dataclasses import dataclass
from functools import cached_property

import anndata as ad
import dask.dataframe as dd
import pandas as pd


@dataclass
class Alias:
    source_name: str

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return getattr(obj, self.source_name)

    def __set__(self, obj, value):
        setattr(obj, self.source_name, value)


def _open(path: os.PathLike):
    return (gzip.open if "gzip" in mimetypes.guess_type(path) else open)(path)


class HapsFile:
    number_of_sequences: int = 0
    number_of_SNPs: int = 0
    N: int = Alias("number_of_sequences")
    L: int = Alias("number_of_SNPs")

    def __init__(self, haps_path: os.PathLike, sample_path: os.PathLike) -> None:
        with _open(sample_path) as fp:
            for line in fp.readlines()[2:]:  # skip two header lines
                id1, id2, _ = line.split()
                self.number_of_sequences += 1 if id1 == id2 else 2
        with _open(haps_path) as fp:
            ...


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
) -> ad.AnnData:
    """Read Oxford phased haplotype file

    Parameters
    ----------
    haps_path : os.PathLike
    sample_path : os.PathLike | None, optional
        default None

    Returns
    -------
    ad.AnnData

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
    sample = read_sample(sample_path)
    df = dd.read_csv(
        haps_path,
        sep=r"\s+",
        assume_missing=True,
        na_values=".",
        names=("CHR", "ID", "POS", "ALLELE0", "ALLELE1", *sample.ids),
        blocksize=None,
    )
    adata = ad.AnnData(df[sample.ids].values.T, dtype=int)
    adata.obs_names = sample.ids
    adata.var_names = df["ID"].astype(str)
    adata.var["CHR"] = df["CHR"]
    adata.var["POS"] = df["POS"]
    adata.var["ALLELE0"] = df["ALLELE0"]
    adata.var["ALLELE1"] = df["ALLELE1"]
    return adata


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
