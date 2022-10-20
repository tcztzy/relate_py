from relatepy.io import read_haps


def test_haps(haps_path, sample_path):
    data = read_haps(haps_path, sample_path)
    assert data.n_obs == data.N == 8
    assert data.n_vars == data.L == 130862
