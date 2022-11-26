import pathlib
import warnings

import numpy as np
from relatepy.io import HapsFile

LOWER_RESCALING_THRESHOLD = 1e-10
UPPER_RESCALING_THRESHOLD = 1e10


class FastPainting:
    def __init__(self, nobs: int, theta: float = 0.001) -> None:
        assert theta < 1.0
        self.theta = theta
        self.ntheta = ntheta = 1.0 - theta
        self.Nminusone = nobs - 1.0
        self.prior_theta = (theta - ntheta) / (nobs - 1)
        self.prior_ntheta = ntheta / (nobs - 1)
        self.theta_ratio = theta / ntheta - 1.0
        self.log_ntheta = np.log(ntheta)
        self.log_small = np.log(0.01)

    def paint_stepping_stones(
        self, data: HapsFile, chunk_index, k, paint_dir: pathlib.Path
    ):
        if data.data.n_obs < 100:
            warnings.warn(f"Sample number `{data.data.n_obs}` is too small.")
        window_boundaries = data.window_boundaries[chunk_index]
        num_windows = len(window_boundaries) - 1
        if window_boundaries[-1] != data.L:
            raise ValueError("")
        last_snp: np.uint32 = data.L - 1
        derived_k = np.flatnonzero(data.data.X[k] == 1)
        if 0 not in derived_k:
            derived_k = np.append(0, derived_k)
        if last_snp not in derived_k:
            derived_k = np.append(derived_k, last_snp)
        num_derived_sites = len(derived_k)
        r_sum = np.clip(
            [_.sum() for _ in np.split(data.r, derived_k[1:])], 0, -np.log(0.01)
        )
        r_prob = np.append(1 - np.exp(-r_sum), 1)
        nor_x_theta = -r_sum + self.log_ntheta

        i = np.searchsorted(derived_k, window_boundaries[1:-1])
        boundary_snp_begin = np.append(0, derived_k[i])
        boundary_snp_end = np.append(window_boundaries[1:-1], last_snp)
        alpha = np.zeros((num_windows, data.N), dtype=np.float32)
        beta = np.zeros((num_windows, data.N), dtype=np.float32)
        logscales_alpha = np.zeros(num_windows, dtype=np.float32)
        logscales_beta = np.zeros(num_windows, dtype=np.float32)

        it1_alpha = 0
        rit1_beta = -1

        it_logscale_alpha = 0
        rit_logscale_beta = -1

        # alpha_aux and beta_aux contain the alpha and beta values along the sequence.
        # I am alternating between two rows, to keep the previous and the current values
        alpha_aux = np.zeros((2, data.N), dtype=np.double)
        logscale: list = [0.0, 0.0]
        derived = data.data.X < data.data.X[k]
        it_boundary_snp_begin = 0
        alpha_sum = self.ntheta / r_prob[0] * (1 - r_prob[0])
        for i, snp in enumerate(derived_k):
            # precalculated quantities
            aux_index = i % 2
            aux_index_prev = 1 - aux_index

            # inner loop of forward algorithm
            if r_prob[i] < 1:
                r = r_prob[i] / ((1 - r_prob[i]) * self.Nminusone)
                logscale[aux_index] = logscale[aux_index_prev] + nor_x_theta[i]
                alpha_aux[aux_index] = alpha_aux[aux_index_prev] + r * alpha_sum
                alpha_aux[aux_index] *= np.where(
                    derived[:, snp], self.theta / self.ntheta, 1
                )
            else:
                r = 1
                logscale[aux_index_prev] = logscale[aux_index]
                logscale[aux_index] += np.log(self.ntheta / self.Nminusone * alpha_sum)
                alpha_aux[aux_index] = np.where(
                    derived[:, snp], self.theta / self.ntheta, 1
                )
            alpha_aux[aux_index, k] = 0.0
            alpha_sum = alpha_aux[aux_index].sum()

            # check if alpha_sums get too small, if they do, rescale
            if (
                alpha_sum < LOWER_RESCALING_THRESHOLD
                or alpha_sum > UPPER_RESCALING_THRESHOLD
            ):
                alpha_aux[aux_index] /= alpha_sum
                logscale[aux_index] += np.log(alpha_sum)
                alpha_sum = 1.0

            if it_boundary_snp_begin != len(boundary_snp_begin):
                # store first the end boundary of the current chunk and then the start boundary of the next chunk
                # copy to alpha, logscale
                while boundary_snp_begin[it_boundary_snp_begin] == snp:
                    alpha[it1_alpha] = alpha_aux[aux_index]
                    logscales_alpha[it_logscale_alpha] = logscale[aux_index]
                    it1_alpha += 1
                    it_logscale_alpha += 1
                    it_boundary_snp_begin += 1
                    if it_boundary_snp_begin == len(boundary_snp_begin):
                        break
        # Backward algorithm
        normalizing_constant = np.double(
            np.log(self.Nminusone) - num_derived_sites * self.log_ntheta
        )
        logscale = [normalizing_constant, normalizing_constant]
        beta_aux = np.zeros((2, data.N), dtype=np.double)
        beta_aux[aux_index] = 1.0
        beta_sum = np.where(derived[:, last_snp], self.theta, self.ntheta).sum() - self.ntheta
        rit_boundarySNP_end = len(boundary_snp_end) - 1
        for i, snp in enumerate(reversed(derived_k)):
            i = num_derived_sites - 1 - i
            aux_index = i % 2
            aux_index_prev = 1 - aux_index

            if r_prob[i] < 1.0:
                r = r_prob[i] / ((1.0 - r_prob[i]) * self.Nminusone)
                # inner loop of backwards algorithm
                logscale[aux_index] = logscale[aux_index_prev] + nor_x_theta[i]
                beta_aux[aux_index] = (
                    beta_aux[aux_index_prev]
                    + r * beta_sum
                    / np.where(derived[:, snp], self.theta, self.ntheta)
                )
                beta_aux[aux_index] *= np.where(
                    derived[:, snp], self.theta/self.ntheta, 1
                )
            else:
                logscale[aux_index_prev] = logscale[aux_index]
                logscale[aux_index_prev] += np.log(
                    self.ntheta / self.Nminusone * alpha_sum
                )
            beta_aux[aux_index, k] = 0
            beta_sum = (np.where(derived[:, snp], self.theta, self.ntheta) * beta_aux[aux_index]).sum()
            if (
                beta_sum < LOWER_RESCALING_THRESHOLD
                or beta_sum > UPPER_RESCALING_THRESHOLD
            ):
                beta_aux[aux_index] /= beta_sum
                logscale[aux_index] += np.log(beta_sum)
                beta_sum = 1.0

            # update beta and topology
            if rit_boundarySNP_end >= 0:
                while boundary_snp_end[rit_boundarySNP_end] == snp:
                    # store first the start boundary of the current chunk and then the end boundary of the next chunk
                    beta[rit1_beta] = beta_aux[aux_index]
                    logscales_beta[rit_logscale_beta] = logscale[aux_index]
                    rit1_beta -= 1
                    rit_logscale_beta -= 1
                    rit_boundarySNP_end -= 1
                    if rit_boundarySNP_end == -1:
                        break

        # Dump to file
        for i in range(num_windows):
            startinterval = window_boundaries[i]
            endinterval = window_boundaries[i + 1] - 1
            pfile = paint_dir / f"relate_{i}.bin"
            with pfile.open("ab") as fp:
                fp.write(
                    np.array([startinterval, endinterval], dtype=np.int32).tobytes()
                )
                dump_to_file(fp, alpha[i], boundary_snp_begin[i], logscales_alpha[i])
                dump_to_file(fp, beta[i], boundary_snp_end[i], logscales_beta[i])


def dump_to_file(f, a, b, l):
    f.write(np.array([1, len(a)], dtype=np.uint64).tobytes())
    f.write(np.int32(b).tobytes())
    f.write(np.float32(l).tobytes())
    f.write(a.tobytes())


def paint(
    data: HapsFile,
    chunk_index: int,
    output: pathlib.Path,
    theta: float | None = None,
    rho: float | None = None,
):
    if theta is not None:
        data.theta = theta
        data.ntheta = 1 - theta
    else:
        theta = 0.001
    if rho is not None:
        data.r *= rho
    else:
        rho = 1.0

    (paint_dir := output / f"chunk_{chunk_index}" / "paint").mkdir(parents=True)
    painter = FastPainting(data.data.n_obs, theta)
    for hap in range(data.N):
        painter.paint_stepping_stones(data, chunk_index, hap, paint_dir)
