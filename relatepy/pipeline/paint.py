from math import exp, log
import pathlib
from relatepy.io import HapsFile
import numpy as np


class FastPainting:
    lower_rescaling_threshold = 1e-10
    upper_rescaling_threshold = 1.0 / lower_rescaling_threshold

    def __init__(self, nobs: int, theta: float = 0.025) -> None:
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
        window_boundaries = data.window_boundaries[chunk_index]
        num_windows = len(window_boundaries) - 1
        boundary_snp_begin = np.zeros(num_windows, dtype=int)
        boundary_snp_end = np.zeros(num_windows, dtype=int)
        window_index = 1
        window_end = window_boundaries[1]
        assert window_boundaries[-1] == data.L

        it_boundary_snp_begin = 1
        it_boundary_snp_end = 0

        r_prob = np.zeros(data.L, dtype=np.double)
        nor_x_theta = np.zeros(data.L, dtype=np.double)
        derived_k = np.zeros(data.L, dtype=np.int32)
        it_r_prob = 0
        it_nor_x_theta = 0

        seq_k: bytes
        snp_next: int
        last_snp: np.uint32 = data.L - 1

        derived_k0 = derived_k[0] = 0
        num_derived_sites: int = 1
        # np.where is used to find the first SNP equal to 1, then plus one because
        # `data.r[0]` is always include, so we should start at `data.r[1]`

        snp = np.where(data.data.X[k, 1:] == 1)[0]
        if snp.size == 0:
            snp = last_snp
        else:
            snp = snp[0] + 1
        # snp: np.uint64 = min(np.where(data.data.X[k, 1:] == 1)[0] + 1, last_snp)
        # I (Tang Ziya) guess the original codes want to sum all recombination rate not
        # equal to 1 from the very beginning.
        r_prob[it_r_prob] = data.r[0] + data.r[1:snp].sum()

        if derived_k0 < window_end <= snp:
            while window_end <= snp:
                boundary_snp_end[it_boundary_snp_end] = snp
                it_boundary_snp_end += 1
                boundary_snp_begin[it_boundary_snp_begin] = derived_k0
                it_boundary_snp_begin += 1
                window_index += 1
                window_end = window_boundaries[window_index]
        it_derived_k = 1

        nor_x_theta[it_nor_x_theta] = -r_prob[it_r_prob] + self.log_ntheta
        r_prob[it_r_prob] = 1.0 - exp(-r_prob[it_r_prob])
        if r_prob[it_r_prob] > 0.99:
            r_prob[it_r_prob] = 0.99
            nor_x_theta[it_nor_x_theta] = self.log_small + self.log_ntheta

        it_r_prob += 1
        it_nor_x_theta += 1

        r_prob[it_r_prob] = data.r[snp]
        derived_k[it_derived_k] = snp

        snp += 1
        num_derived_sites += 1
        while snp < data.L:
            # skip all non-derived site
            while data.data.X[k, snp] != 1 and snp != last_snp:
                r_prob[it_r_prob] += data.r[snp]
                snp += 1

            if derived_k[it_derived_k] < window_end <= snp:
                while window_end <= snp:
                    boundary_snp_end[it_boundary_snp_end] = snp
                    it_boundary_snp_end += 1
                    boundary_snp_begin[it_boundary_snp_begin] = derived_k[it_derived_k]
                    it_boundary_snp_begin += 1
                    window_index += 1
                    window_end = window_boundaries[window_index]

            # this is a derived site
            it_derived_k += 1

            nor_x_theta[it_nor_x_theta] = -r_prob[it_r_prob] + self.log_ntheta
            r_prob[it_r_prob] = 1.0 - exp(-r_prob[it_r_prob])
            if r_prob[it_r_prob] > 0.99:
                r_prob[it_r_prob] = 0.99
                nor_x_theta[it_nor_x_theta] = self.log_small + self.log_ntheta

            it_r_prob += 1
            it_nor_x_theta += 1

            derived_k[it_derived_k] = snp
            r_prob[it_r_prob] = data.r[snp]

            snp += 1
            num_derived_sites += 1
        nor_x_theta[it_nor_x_theta] = -r_prob[it_r_prob] + self.log_ntheta
        r_prob[it_r_prob] = 1.0 - exp(-r_prob[it_r_prob])
        if r_prob[it_r_prob] > 0.99:
            r_prob[it_r_prob] = 0.99
            nor_x_theta[it_nor_x_theta] = self.log_small + self.log_ntheta

        it_r_prob += 1
        r_prob[it_r_prob] = 1  # just a technicality

        derived_k = derived_k[:num_derived_sites]
        r_prob = r_prob[:num_derived_sites+1]
        nor_x_theta = nor_x_theta[:num_derived_sites]

        boundary_snp_end[it_boundary_snp_end] = last_snp
        it_boundary_snp_end += 1
        assert it_boundary_snp_begin == len(boundary_snp_begin)
        assert it_boundary_snp_end == len(boundary_snp_end)

        it_derived_k = 1
        it_r_prob = 0
        it_nor_x_theta = 0

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
        beta_aux = np.zeros((2, data.N), dtype=np.double)
        alpha_aux_rowbegin = np.array([0, data.N], dtype=np.uint32)
        alpha_aux_rowend = np.array([data.N, data.N], dtype=np.uint32)
        beta_aux_rowbegin = np.array([0, data.N], dtype=np.uint32)
        beta_aux_rowend = np.array([data.N, data.N], dtype=np.uint32)
        aux_index = 0
        aux_index_prev = 1

        beta_sum: np.double
        logscale: list = [0.0, 0.0]

        snp = derived_k[0]

        seq_k = data.data.X[k, snp]
        derived = (data.data.X[:, snp] < seq_k).astype(np.double)
        alpha_aux[aux_index] = np.where(
            np.arange(data.data.n_obs) == k,
            0,
            derived * self.prior_theta + self.prior_ntheta,
        )

        alpha_sum = alpha_aux[aux_index].sum()

        it_boundary_snp_begin = 0
        while boundary_snp_begin[it_boundary_snp_begin] == snp:
            # store the start boundary of the current chunk

            # copy to alpha, logscale
            it2_alpha_aux = alpha_aux_rowbegin[aux_index]
            it2_alpha = 0
            while it2_alpha_aux != alpha_aux_rowend[aux_index]:
                alpha[it1_alpha, it2_alpha] = alpha_aux[aux_index, it2_alpha_aux]
                it2_alpha_aux += 1
                it2_alpha += 1
            logscales_alpha[it_logscale_alpha] = logscale[aux_index]
            it_logscale_alpha += 1
            it1_alpha += 1
            it_boundary_snp_begin += 1
            if it_boundary_snp_begin == len(boundary_snp_begin):
                break

        r_x_alpha_sum: np.double
        if r_prob[it_r_prob] < 1.0:
            r_x_alpha_sum = (
                r_prob[it_r_prob]
                / ((1.0 - r_prob[it_r_prob]) * self.Nminusone)
                * alpha_sum
            )
        else:
            r_x_alpha_sum = alpha_sum
        while it_derived_k != len(derived_k):
            # precalculated quantities
            snp = derived_k[it_derived_k]
            it2_sequence = 0
            seq_k = data.data.X[k, snp]
            aux_index = (aux_index + 1) % 2
            aux_index_prev = 1 - aux_index

            # inner loop of forward algorithm
            if r_prob[it_r_prob] < 1:
                logscale[aux_index] = logscale[aux_index_prev]
                logscale[aux_index] += nor_x_theta[it_nor_x_theta]

                it2_sequence = 0
                it2_alpha_aux_prev = alpha_aux_rowbegin[aux_index_prev]

                it2_alpha_aux = alpha_aux_rowbegin[aux_index]

                while it2_alpha_aux != alpha_aux_rowend[aux_index]:
                    alpha_aux[aux_index, it2_alpha_aux] = (
                        alpha_aux[aux_index_prev, it2_alpha_aux_prev] + r_x_alpha_sum
                    )
                    derived = np.double(seq_k > data.data.X[it2_sequence, snp])
                    alpha_aux[aux_index, it2_alpha_aux] *= (
                        derived * self.theta_ratio + 1
                    )
                    it2_alpha_aux += 1
                    it2_alpha_aux_prev += 1
                    it2_sequence += 1

                it2_alpha_aux = alpha_aux_rowbegin[aux_index]
                alpha_aux[aux_index, k] = 0.0
                alpha_sum = alpha_aux[aux_index].sum()
            else:
                logscale[aux_index_prev] = logscale[aux_index]
                logscale[aux_index] += log(self.ntheta / self.Nminusone * alpha_sum)
                it2_alpha_aux = alpha_aux_rowbegin[aux_index]
                it2_sequence = 0
                while it2_alpha_aux != alpha_aux_rowend[aux_index]:
                    derived = np.double(seq_k > data.data.X[it2_sequence, snp])
                    alpha_aux[aux_index, it2_alpha_aux] = (
                        derived * self.theta_ratio + 1.0
                    )
                    it2_alpha_aux += 1
                    it2_sequence += 1
                it2_alpha_aux = alpha_aux_rowbegin[aux_index]
                alpha_aux[aux_index, k] = 0.0
                alpha_sum = alpha_aux[aux_index].sum()

            r_x_alpha_sum = alpha_sum

            # check if alpha_sums get too small, if they do, rescale
            if (
                r_x_alpha_sum < self.lower_rescaling_threshold
                or r_x_alpha_sum > self.upper_rescaling_threshold
            ):
                tmp = r_x_alpha_sum
                alpha_aux[aux_index] /= tmp
                logscale[aux_index] += log(tmp)
                r_x_alpha_sum = 1.0
                assert logscale[aux_index] < float("inf")

            it_r_prob += 1
            if r_prob[it_r_prob] < 1.0:
                r_x_alpha_sum *= r_prob[it_r_prob] / (
                    (1 - r_prob[it_r_prob]) * self.Nminusone
                )

            if it_boundary_snp_begin != len(boundary_snp_begin):
                # store first the end boundary of the current chunk and then the start boundary of the next chunk
                # copy to alpha, logscale
                it2_alpha_aux = alpha_aux_rowbegin[aux_index]
                it2_alpha = 0
                while it2_alpha_aux != alpha_aux_rowend[aux_index]:
                    alpha[it1_alpha, it2_alpha] = alpha_aux[aux_index, it2_alpha_aux]
                    it2_alpha_aux += 1
                    it2_alpha += 1
                logscales_alpha[it_logscale_alpha] = logscale[aux_index]
                it1_alpha += 1
                it_logscale_alpha += 1

                it_boundary_snp_begin += 1
                if it_boundary_snp_begin == len(boundary_snp_begin) - 1:
                    break

            it_derived_k += 1
            it_nor_x_theta += 1

        assert it_boundary_snp_begin == len(boundary_snp_begin) - 1
        assert it1_alpha == alpha.shape[0] - 1
        assert it_logscale_alpha == len(logscales_alpha) - 1

        # Backward algorithm

        normalizing_constant = np.double(
            log(self.Nminusone) - num_derived_sites * self.log_ntheta
        )

        # SNP L-1
        it_derived_k -= 1
        assert last_snp == derived_k[it_derived_k]

        logscale[aux_index] = normalizing_constant
        logscale[aux_index_prev] = normalizing_constant
        beta_sum = 0.0

        it2_sequence = last_snp
        seq_k = data.data.X[k, it2_sequence]

        beta_aux_rowbegin[aux_index, :] = 1.0
        it2_beta_aux = beta_aux_rowbegin[aux_index]
        while it2_beta_aux != beta_aux_rowend[aux_index]:
            if seq_k > data.data.X[it2_sequence, last_snp]:
                beta_sum += self.theta
            else:
                beta_sum += self.ntheta
            it2_sequence += 1
            it2_beta_aux += 1
        beta_sum -= self.ntheta

        rit_boundarySNP_end = len(boundary_snp_end) - 1
        while boundary_snp_end[rit_boundarySNP_end] == last_snp:
            # copy beta
            it2_beta = 0
            it2_beta_aux = 0
            while it2_beta_aux != beta_aux_rowend[aux_index]:
                beta[rit1_beta, it2_beta] = beta_aux[aux_index, it2_beta_aux]
                it2_beta += 1
                it2_beta_aux += 1
            logscales_beta[rit_logscale_beta] = logscale[aux_index]
            rit1_beta -= 1
            rit_logscale_beta -= 1
            rit_boundarySNP_end -= 1
            if rit_boundarySNP_end == 0:
                break

        # SNP < L-1
        if r_prob[it_r_prob] < 1.0:
            r_x_beta_sum = (
                r_prob[it_r_prob]
                / ((1.0 - r_prob[it_r_prob]) * self.Nminusone)
                * beta_sum
            )
        else:
            r_x_alpha_sum = beta_sum
        snp_next = last_snp

        while snp > 0:
            it_derived_k -= 1
            snp = derived_k[it_derived_k]
            aux_index = (aux_index + 1) % 2
            aux_index_prev = 1 - aux_index

            if r_prob[it_r_prob] < 1.0:

                # inner loop of backwards algorithm
                logscale[aux_index] = logscale[aux_index_prev]
                logscale[aux_index] += nor_x_theta[it_nor_x_theta]
                beta_sum = 0
                beta_sum_oneminustheta = r_x_beta_sum / self.ntheta
                beta_sum_theta = r_x_beta_sum / self.theta - beta_sum_oneminustheta
                it2_sequence = 0

                it2_beta_aux_next = beta_aux_rowbegin[aux_index_prev]
                it2_beta_aux = beta_aux_rowbegin[aux_index]
                while it2_beta_aux != beta_aux_rowend[aux_index]:
                    derived = np.double(seq_k > data.data.X[it2_sequence, snp_next])
                    beta_aux[aux_index, it2_beta_aux] = (
                        beta_aux[aux_index_prev, it2_beta_aux_next]
                        + derived * beta_sum_theta
                        + beta_sum_oneminustheta
                    )
                    beta_aux[aux_index, it2_beta_aux] *= (
                        derived * self.theta_ratio + 1.0
                    )
                    it2_sequence += 1
                    it2_beta_aux += 1
                    it2_beta_aux_next += 1

                it2_sequence = 0
                seq_k = data.data.X[k, snp]
                it2_beta_aux = beta_aux_rowbegin[aux_index]
                beta_aux[aux_index, k] = 0
                while it2_beta_aux != beta_aux_rowend[aux_index]:
                    if seq_k > data.data.X[it2_sequence, snp]:
                        beta_sum += self.theta * beta_aux[aux_index, it2_beta_aux]
                    else:
                        beta_sum += self.ntheta * beta_aux[aux_index, it2_beta_aux]
                    it2_sequence += 1
                    it2_beta_aux += 1
            else:
                logscale[aux_index_prev] = logscale[aux_index]
                logscale[aux_index_prev] += log(
                    self.ntheta / self.Nminusone * alpha_sum
                )
                beta_sum = 0.0

                it2_beta_aux = beta_aux_rowbegin[aux_index]
                while it2_beta_aux != beta_aux_rowend[aux_index]:
                    if seq_k > data.data.X[it2_sequence, snp]:
                        beta_sum += self.theta * beta_aux[aux_index, it2_beta_aux]
                    else:
                        beta_sum += self.ntheta * beta_aux[aux_index, it2_beta_aux]
                    it2_sequence += 1
                    it2_beta_aux += 1

            r_x_beta_sum = beta_sum

            if (
                r_x_beta_sum < self.lower_rescaling_threshold
                or r_x_beta_sum > self.upper_rescaling_threshold
            ):
                tmp = r_x_beta_sum
                it2_beta_aux = beta_aux_rowbegin[aux_index]
                while it2_beta_aux < beta_aux_rowend[aux_index]:
                    beta_aux[aux_index, it2_beta_aux] /= tmp
                    it2_beta_aux += 1
                logscale[aux_index] += log(tmp)
                r_x_beta_sum = 1.0
                assert logscale[aux_index] < float("inf")

            it_r_prob -= 1
            if r_prob[it_r_prob] < 1.0:
                r_x_beta_sum *= r_prob[it_r_prob] / (
                    (1.0 - r_prob[it_r_prob]) * self.Nminusone
                )

            # update beta and topology
            if rit_boundarySNP_end != 0:
                while boundary_snp_end[rit_boundarySNP_end] == snp:
                    # store first the start boundary of the current chunk and then the end boundary of the next chunk

                    it2_beta = 0
                    it2_beta_aux = beta_aux_rowbegin[aux_index]
                    while it2_beta_aux != beta_aux_rowend[aux_index]:
                        beta[rit1_beta, it2_beta] = beta_aux[aux_index, it2_beta_aux]
                        it2_beta += 1
                        it2_beta_aux += 1

                    logscales_beta[rit_logscale_beta] = logscale[aux_index]
                    rit1_beta -= 1
                    rit_logscale_beta -= 1
                    rit_boundarySNP_end -= 1
                    if rit_boundarySNP_end == 0:
                        break

            snp_next = snp
            it_nor_x_theta -= 1  # I want this to be pointing at snp_next

        assert derived_k[it_derived_k] == 0
        assert rit_boundarySNP_end == 0

        # Dump to file
        for i in range(num_windows):
            startinterval = window_boundaries[i]
            endinterval = window_boundaries[i + 1] - 1
            pfile = paint_dir / f"relate_{i}.bin"
            with pfile.open("wb") as fp:
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
        theta = 0.025
    if rho is not None:
        data.r *= rho
    else:
        rho = 1.0

    (paint_dir := output / f"chunk_{chunk_index}" / "paint").mkdir(parents=True)
    painter = FastPainting(data.data.n_obs, theta)
    for hap in range(data.N):
        painter.paint_stepping_stones(data, chunk_index, hap, paint_dir)
