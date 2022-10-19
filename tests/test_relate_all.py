from pathlib import Path

from tempfile import TemporaryDirectory
from relatepy import all_pipeline


def test_relate():
    data_dir = (Path("relate") / "example" / "data").resolve()
    with TemporaryDirectory() as tmp:
        all_pipeline(
            haps=data_dir / "example.haps.gz",
            sample=data_dir / "example.sample.gz",
            genetic_map=data_dir / "genetic_map_GRCh37_chr1.txt",
            output=Path(tmp) / "example",
            mutation_rate=1.25e-8,
            effective_population_size=30000,
            sample_ages=data_dir / "sample_ages.txt",
            seed=1,
        )
