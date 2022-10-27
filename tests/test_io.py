from pathlib import Path
from relatepy.io import read_haps
from struct import calcsize, unpack
import numpy as np


def test_haps(haps_path, sample_path, genetic_map_path, tmp_path: Path):
    data = read_haps(haps_path, sample_path)
    assert data.N == 8 and data.L == 130862
    output_dir = tmp_path / "example"
    output_dir.mkdir()
    data.make_chunks(output_dir, genetic_map_path)
    # check chunk_0.bp
    chunk_0_bp = output_dir / "chunk_0.bp"
    assert chunk_0_bp.exists()
    content = chunk_0_bp.read_bytes()
    assert len(content) == data.L * 4 + 4
    (L_chunk,) = unpack("I", content[:4])
    assert L_chunk == data.L
    assert (unpack(f"{L_chunk}I", content[4:]) == data.bp_pos[:-1]).all()
    # check chunk_0.hap
    chunk_0_hap = output_dir / "chunk_0.hap"
    assert chunk_0_hap.exists()
    content = chunk_0_hap.read_bytes()
    fmt = "QQ"
    size = calcsize(fmt)
    L_chunk, N = unpack(fmt, content[:size])
    assert L_chunk == data.L and N == data.N and len(content) == L_chunk * N + size
    assert ([int(chr(c)) for c in content[size:]] == data._data.X.T.reshape(-1)).all()
    # check chunk_0.dist
    chunk_0_dist = output_dir / "chunk_0.dist"
    assert chunk_0_dist.exists()
    content = chunk_0_dist.read_bytes()
    sizeof_uint32 = calcsize("I")
    (L,) = unpack("I", content[:sizeof_uint32])
    assert L == data.L and len(content) == sizeof_uint32 + calcsize("i") * L
    assert (unpack("i" * L, content[sizeof_uint32:]) == np.diff(data.bp_pos)).all()
    # check chunk_0.rpos
    chunk_0_rpos = output_dir / "chunk_0.rpos"
    assert chunk_0_rpos.exists()
    content = chunk_0_rpos.read_bytes()
    fmt = f"=I{data.L+1}d"
    length, *rpos = unpack(fmt, content)
    assert len(rpos) == length
    assert (data.rpos == rpos).all()
    # check chunk_0.r
    chunk_0_r = output_dir / "chunk_0.r"
    assert chunk_0_r.exists()
    content = chunk_0_r.read_bytes()
    fmt = f"=I{data.L}d"
    length, *r = unpack(fmt, content)
    assert len(r) == length
    assert (data.r == r).all()
    # check chunk_0.state
