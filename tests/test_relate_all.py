import os
from pathlib import Path

from tempfile import TemporaryDirectory
from relatepy import all_pipeline


def test_relate():
    data_dir = (Path("relate") / "example" / "data").resolve()
    with TemporaryDirectory() as tmp:
        cwd = os.getcwd()
        os.chdir(tmp)
        all_pipeline(
            data_dir / "example.haps.gz",
            data_dir / "example.sample.gz",
            data_dir / "genetic_map_GRCh37_chr1.txt",
            Path(tmp) / "example",
            1.25e-8,
            30000,
            data_dir / "sample_ages.txt",
            seed=1,
        )
        os.chdir(cwd)

if __name__ == "__main__":
    test_relate()
