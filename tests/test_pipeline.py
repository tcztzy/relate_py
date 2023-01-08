from pathlib import Path
from tempfile import TemporaryDirectory
import struct

import numpy as np

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


def test_paint(haps_path, sample_path, genetic_map_path, paint_bin, tmp_path: Path):
    (output_path := tmp_path / "output").mkdir()
    data = HapsFile(haps_path, sample_path)
    data.make_chunks(output_path, genetic_map_path)
    paint(output=output_path, chunk_index=0)
    paint_result = output_path / "chunk_0" / "paint" / "relate_0.bin"
    assert paint_result.exists()
    content = paint_result.read_bytes()
    assert len(content) == 960
    fmt = "IIQQIf"
    size = struct.calcsize(fmt)
    start, stop, one, num_samples, boundary_begin, logscale = struct.unpack(fmt, content[:size])
    assert start == 0 and stop == 130861 and one == 1 and num_samples == 8 and boundary_begin == 0 and logscale == 0
    fmt = "f" * num_samples
    offset, size = size, struct.calcsize(fmt)
    alpha = struct.unpack(fmt, content[offset:offset+size])
    assert alpha[0] == 0 and all(map(lambda x: x == 0.14271429181098938, alpha[1:]))
    offset += size
    fmt = "QQIf"
    size = struct.calcsize(fmt)
    one, num_samples, boundary_end, logscale = struct.unpack(fmt, content[offset:offset+size])
    assert one == 1 and num_samples == 8 and boundary_end == 130861 and logscale == 48.57122802734375
    fmt = "f" * num_samples
    offset, size = offset + size, struct.calcsize(fmt)
    beta = struct.unpack(fmt, content[offset:offset+size])
    assert (np.array(beta) == 1).all()
    assert paint_bin == content
