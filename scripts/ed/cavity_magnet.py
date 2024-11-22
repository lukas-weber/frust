from . import hamiltonian
import numpy as np
import math
import scipy.sparse as sps
from . import downfolded_coupling
from . import model_common


class Model(model_common.Magnet):
    def __init__(self, model_data):
        self.model_data = model_data
        self.photon_dimension = np.prod([m.max_photons for m in model_data.modes])
        self.spin_dimension = np.prod([s.spin_dim for s in model_data.sites])
        self.N = len(model_data.sites)
        self.lifter = hamiltonian.SpinLifter([s.spin_dim for s in model_data.sites])
        self.photon_lifter = hamiltonian.Lifter(
            [m.max_photons for m in model_data.modes]
        )

    def photon_number_ops(self):
        photon_numbers = []

        for i, m in enumerate(self.model_data.modes):
            photon_numbers.append(
                self.photon_lifter.lift_op(i, sps.diags(np.arange(m.max_photons)))
            )

        return photon_numbers

    def hamiltonian(self):
        photon_numbers = self.photon_number_ops()

        photon_identity = sps.identity(self.photon_dimension)
        spin_identity = sps.identity(self.spin_dimension)

        dim = self.photon_dimension * self.spin_dimension
        H = sps.dok_matrix((dim, dim))

        for b in self.model_data.bonds:
            assert len(self.model_data.modes) == 1
            max_photons = self.model_data.modes[0].max_photons
            omega = self.model_data.modes[0].omega
            coupling = b.mode_couplings[0]
            H += b.J * sps.kron(
                downfolded_coupling.matrix(max_photons, omega, coupling),
                (
                    self.lifter.heisen_bond(b.i, b.j)
                    - 0.25 * sps.eye(self.spin_dimension)
                ),
            )

        # for i, s in enumerate(self.model_data.sites):
        #    for j, m in enumerate(self.model_data.modes):
        #        H += m.coupling * sps.kron(photon_numbers[j], self.lifter.Sz(i))

        for j, m in enumerate(self.model_data.modes):
            H += m.omega * sps.kron(photon_numbers[j], spin_identity)
        return H

    def observables(self, params, Ts, E, psi):
        ens = hamiltonian.Ensemble(E, psi, Ts)

        obs = {}
        obs["Energy"] = ens.diag_mean(E) / self.N
        obs["SpecificHeat"] = (
            Ts ** (-2) * (ens.diag_mean(E**2) - ens.diag_mean(E) ** 2) / self.N
        )

        for name, op in self.mag_sign_cfgs.items():
            if name in params["measure"]:
                M = sps.kron(
                    sps.eye(self.photon_dimension),
                    self.signed_magnetization(lambda i: self.lifter.Sz(i), *op[1]),
                )
                obs.update(
                    {
                        op[0] + k: v
                        for k, v in self.observables_magnetization(ens, M).items()
                    }
                )

                photon_numbers = [
                    sps.kron(op, sps.eye(self.spin_dimension))
                    for op in self.photon_number_ops()
                ]

                obs["PhotonNum"] = np.array(
                    [
                        ens.mean(photon_numbers[m])
                        for m in range(len(self.model_data.modes))
                    ]
                ).T
                obs["PhotonNum2"] = np.array(
                    [
                        ens.mean(photon_numbers[m] @ photon_numbers[m])
                        for m in range(len(self.model_data.modes))
                    ]
                ).T
                obs["PhotonNumVar"] = obs["PhotonNum2"] - obs["PhotonNum"] ** 2

        return obs
