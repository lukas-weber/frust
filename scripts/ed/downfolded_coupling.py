import math
import numpy as np
import numba
from numba import prange
import itertools


@numba.jit(nopython=True)
def disp_op(n, m, g) -> float:
    s = 0
    for k in range(min(n, m) + 1):
        s += (1j * g) ** (n + m - 2 * k) * np.exp(
            0.5 * math.lgamma(1 + n)
            + 0.5 * math.lgamma(1 + m)
            - math.lgamma(1 + k)
            - math.lgamma(1 + n - k)
            - math.lgamma(1 + m - k)
        )
    return s * np.exp(-(g**2) / 2)


def elem(ns, ms, omegas, gs) -> float:
    s = 0
    lmaxs = 2 * np.minimum(ns, ms) + 20

    for ls in itertools.product(*(range(lmax) for lmax in lmaxs.astype(int))):
        disp_ops = 1
        for n, m, l, g in zip(ns, ms, ls, gs):
            disp_ops *= (
                disp_op(n, l, g)
                * np.conj(disp_op(l, m, g))
                * (1 - 2 * ((n - m) // 2 % 2))
            )
        s += disp_ops.real * (
            1 / (1 + omegas @ (ls - ms)) + 1 / (1 + omegas @ (ls - ns))
        )

    return 0.5 * s


def matrix(max_photons, omegas, gs) -> np.array:
    photon_dim = np.prod(max_photons)
    res = np.zeros([photon_dim, photon_dim])

    for i, ns in enumerate(itertools.product(*(range(m) for m in max_photons))):
        for j, ms in enumerate(itertools.product(*(range(m) for m in max_photons))):
            res[i, j] = elem(np.array(ns), np.array(ms), np.array(omegas), np.array(gs))
    return res
