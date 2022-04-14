import math
import numpy as np
import numba
from numba import prange


@numba.jit(nopython=True)
def disp_op(n, m, g):
    s = 0
    for k in range(min(n, m) + 1):
        s += (1j * g) ** (n + m - 2 * k) * np.exp(
            0.5 * math.lgamma(1 + n)
            + 0.5 * math.lgamma(1 + m)
            - math.lgamma(1 + k)
            - math.lgamma(1 + n - k)
            - math.lgamma(1 + m - k)
        )
    return s


@numba.jit(nopython=True, parallel=True)
def elem(n, m, omega, g):
    s = 0

    for l in prange(2 * min(n, m) + 20):
        s += (disp_op(n, l, g) * np.conj(disp_op(l, m, g))).real * (
            1 / (1 + omega * (l - m)) + 1 / (1 + omega * (l - n))
        )
    return 0.5 * (1 - 2 * ((n - m) // 2 % 2)) * np.exp(-(g**2)) * s


def matrix(max_bosons, omega, g):
    return np.array(
        [[elem(n, m, omega, g) for n in range(max_bosons)] for m in range(max_bosons)]
    )
