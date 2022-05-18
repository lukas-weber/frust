import cppyy as _cppyy
from pathlib import Path as _Path
import numpy as _np

root = _Path(__file__).resolve().parent

_cppyy.include(str(root / "downfolded_peierls_coupling.h"))
_cppyy.load_library(str(next(root.glob("downfolded_peierls_coupling.*.so"))))

from cppyy.gbl import downfolded_peierls_coupling as _dpc

def elem(m, n, omegas, gs, max_photons):
    return _dpc.elem([_dpc.mode_params(omega, g, max_photons) for omega, g in zip(omegas, gs)])

def matrix(omegas, gs, max_photons):
    j = _np.array(_dpc.matrix([_dpc.mode_params(omega, g, max_photons) for omega, g in zip(omegas, gs)]))
    n = round(len(j) ** 0.5)
    return j.reshape([n, n])
