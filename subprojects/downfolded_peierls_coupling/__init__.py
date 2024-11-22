import cppyy
from pathlib import Path as _Path
import numpy as _np

root = _Path(__file__).resolve().parent

cppyy.include(str(root / "downfolded_peierls_coupling.h"))
cppyy.load_library(str(next(root.glob("downfolded_peierls_coupling.*.so"))))

from cppyy.gbl import downfolded_coupling as _matrix


def matrix(omegas, gs, max_photons):
    j = _np.array(_matrix([(omega, g, max_photons) for omega, g in zip(omegas, gs)]))
    n = round(len(j) ** 0.5)
    return j.reshape([n, n])
