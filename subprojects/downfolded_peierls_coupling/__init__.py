import cppyy as _cppyy
from pathlib import Path as _Path
import numpy as _np

root = _Path(__file__).resolve().parent

_cppyy.include(str(root / "downfolded_peierls_coupling.h"))
_cppyy.load_library(str(next(root.glob("downfolded_peierls_coupling.*.so"))))

from cppyy.gbl import downfolded_peierls_coupling as _dpc


class Generator:
    def __init__(self, omegas, gs, max_photons):
        self._impl = _dpc.generator(
            [_dpc.mode_params(omega, g, max_photons) for omega, g in zip(omegas, gs)]
        )

    def elem(self, m, n):
        return self._impl.elem(m, n)

    def matrix(self):
        res = _np.array(self._impl.matrix())
        n = round(len(res) ** 0.5)
        return res.reshape([n, n])
