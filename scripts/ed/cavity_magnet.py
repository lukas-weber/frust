from . import hamiltonian
import numpy as np
import math
import scipy.sparse as sps


class Model:
    def __init__(self, model_data):
        self.model_data = model_data
        self.boson_dimension = np.prod([m.max_bosons for m in model_data.modes])
        self.spin_dimension = np.prod([s.spin_dim for s in model_data.sites])
        self.N = len(model_data.sites)
        self.lifter = hamiltonian.SpinLifter([s.spin_dim for s in model_data.sites])
        self.boson_lifter = hamiltonian.Lifter([m.max_bosons for m in model_data.modes])

    def _downfolded_coupling(self, mode):
        def disp_op(n, m, g):
            return sum(
                (1j * g) ** (n + m - 2 * k)
                * np.exp(
                    0.5 * math.lgamma(1 + n)
                    + 0.5 * math.lgamma(1 + m)
                    - math.lgamma(1 + k)
                    - math.lgamma(1 + n - k)
                    - math.lgamma(1 + m - k)
                )
                for k in range(min(n, m) + 1)
            )

        def elem(n, m, omega, g):
            return (
                0.5
                * np.exp(-(g**2))
                * sum(
                    (disp_op(n, l, g) * np.conj(disp_op(l, m, g))).real
                    * (1 / (1 + omega * (l - m)) + 1 / (1 + omega * (l - n)))
                    for l in range(2 * min(n, m) + 20)
                )
            )

        return np.array(
            [
                [elem(n, m, mode.omega, mode.coupling) for n in range(mode.max_bosons)]
                for m in range(mode.max_bosons)
            ]
        )

    def hamiltonian(self):
        boson_numbers = []
        dim = 1

        for m in self.model_data.modes:
            boson_numbers.append(
                sps.kron(
                    sps.identity(dim),
                    sps.kron(
                        sps.diags(np.arange(m.max_bosons)),
                        sps.identity(self.boson_dimension / m.max_bosons / dim),
                    ),
                )
            )
            dim *= m.max_bosons

        boson_identity = sps.identity(self.boson_dimension)
        spin_identity = sps.identity(self.spin_dimension)

        dim = self.boson_dimension * self.spin_dimension
        H = sps.dok_matrix((dim, dim))

        for b in self.model_data.bonds:
            assert len(self.model_data.modes) == 1
            H += b.J * sps.kron(
                self._downfolded_coupling(self.model_data.modes[0]),
                self.lifter.heisen_bond(b.i, b.j),
            )

        # for i, s in enumerate(self.model_data.sites):
        #    for j, m in enumerate(self.model_data.modes):
        #        H += m.coupling * sps.kron(boson_numbers[j], self.lifter.Sz(i))

        for j, m in enumerate(self.model_data.modes):
            H += m.omega * sps.kron(boson_numbers[j], spin_identity)
        return H

    def observables(self, params, Ts, E, psi):
        ens = hamiltonian.Ensemble(E, psi, Ts)

        obs = {}
        obs["Energy"] = ens.diag_mean(E) / self.N
        obs["SpecificHeat"] = (
            Ts ** (-2) * (ens.diag_mean(E**2) - ens.diag_mean(E) ** 2) / self.N
        )

        return obs
