from pathlib import Path
from tempfile import TemporaryDirectory

from relatepy import all_pipeline
from relatepy.pipeline.paint import paint
from relatepy.io import HapsFile


def test_relate(haps_path, sample_path, genetic_map_path, sample_ages_path):
    with TemporaryDirectory() as tmp:
        all_pipeline(
            haps=haps_path,
            sample=sample_path,
            genetic_map=genetic_map_path,
            output=Path(tmp) / "example",
            mutation_rate=1.25e-8,
            effective_population_size=30000,
            sample_ages=sample_ages_path,
            seed=1,
        )


def test_paint(haps_path, sample_path, genetic_map_path, tmp_path: Path):
    (output_path := tmp_path / "output").mkdir()
    data = HapsFile(haps_path, sample_path)
    data.make_chunks(output_path, genetic_map_path)
    paint(data, 0, output_path)
