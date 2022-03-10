from . import hamiltonian
import numpy as np
import scipy.sparse as sps


class Model:
    def __init__(self, model_data):
        self.model_data = model_data
        self.Nfull = len(self.model_data.sites)
        half2full = sum(
            [site.nspinhalfs * [i] for i, site in enumerate(self.model_data.sites)], []
        )
        self.full2half = [
            [i for i, x in enumerate(half2full) if x == idx]
            for idx in range(self.Nfull)
        ]
        self.N = len(half2full)
        self.lifter = hamiltonian.SpinLifter([2] * self.N)

    def chirality_operators(self):
        j = np.exp(2j * np.pi / 3)

        uL = 1 / 3**0.5 * np.matrix([0, 1, j, 0, j * j, 0, 0, 0])
        uR = 1 / 3**0.5 * np.matrix([0, 1, j * j, 0, j, 0, 0, 0])
        dL = 1 / 3**0.5 * np.matrix([0, 0, 0, j * j, 0, j, 1, 0])
        dR = 1 / 3**0.5 * np.matrix([0, 0, 0, j, 0, j * j, 1, 0])

        us = 1 / 2**0.5 * np.matrix([0, 0, 1, 0, -1, 0, 0, 0])
        ut = 1 / 6**0.5 * np.matrix([0, -2, 1, 0, 1, 0, 0, 0])
        ds = 1 / 2**0.5 * np.matrix([0, 0, 0, -1, 0, 1, 0, 0])
        dt = 1 / 6**0.5 * np.matrix([0, 0, 0, 1, 0, 1, -2, 0])

        tauz = uL.H @ uL + dL.H @ dL - uR.H @ uR - dR.H @ dR
        tauy = -1j * (uL.H @ uR + dL.H @ dR - uR.H @ uL - dR.H @ dL)

        def lift_tri(i, op):
            return sps.kron(
                sps.kron(sps.identity(8**i), op),
                sps.identity(8 ** (self.Nfull - i - 1)),
            )

        nem_diag = ut.H @ ut - us.H @ us + dt.H @ dt - ds.H @ ds
        nem_off = 3 / 4 * (ut.H @ us + us.H @ ut + dt.H @ ds + ds.H @ dt)

        nem_diag_corr = []
        tauzcorr = []
        tauycorr = []

        nem_off_corr = []
        nem_cross_corr = []
        for i in range(self.Nfull):
            tauzcorr.append(lift_tri(0, tauz) @ lift_tri(i, tauz))
            tauycorr.append(lift_tri(0, tauy) @ lift_tri(i, tauy))
            nem_diag_corr.append(lift_tri(0, nem_diag) @ lift_tri(i, nem_diag))
            nem_off_corr.append(lift_tri(0, nem_off) @ lift_tri(i, nem_off.H))
            nem_cross_corr.append(
                1j
                * (
                    lift_tri(0, nem_diag) @ lift_tri(i, nem_off.H)
                    - lift_tri(0, nem_off) @ lift_tri(i, nem_diag.H)
                )
            )

        ops = {}
        ops["TauZ"] = tauzcorr
        ops["TauY"] = tauycorr
        ops["NematicityDiagCorr"] = nem_diag_corr
        ops["NematicityOffCorr"] = nem_off_corr
        ops["NematicityCrossCorr"] = nem_cross_corr

        C = (
            self.lifter.heisen_bond(0, 1)
            + j * self.lifter.heisen_bond(1, 2)
            + j**2 * self.lifter.heisen_bond(2, 0)
        )

        ops["NematicityCorr"] = [
            C
            @ (
                self.lifter.heisen_bond(3 * i + 0, 3 * i + 1)
                + j * self.lifter.heisen_bond(3 * i + 1, 3 * i + 2)
                + j**2 * self.lifter.heisen_bond(3 * i + 2, 3 * i + 0)
            ).H
            for i in range(0, self.Nfull)
        ]
        return ops

    def signed_mag(self, sx, sy, suc):
        dim = 2**self.N
        M = sps.dok_matrix((dim, dim))
        for x in range(self.model_data.Lx):
            for y in range(self.model_data.Ly):
                for uc in range(self.model_data.uc_site_count):
                    i = (
                        self.model_data.uc_site_count * (self.model_data.Lx * y + x)
                        + uc
                    )
                    sign = (
                        sx**x
                        * sy**y
                        * (self.model_data.sites[i].sublattice_sign if suc < 0 else 1)
                    )

                    for idx in self.full2half[i]:
                        M += sign * self.lifter.Sz(idx)
        return M / self.N

    def hamiltonian(self):
        dim = 2**self.N
        Id = sps.identity(dim)

        def onsite_term(Jin, site, h):
            res = sps.dok_matrix((dim, dim))
            idx = 0
            for i in self.full2half[site]:
                res += self.lifter.Sz(i) * h
                for j in self.full2half[site]:
                    if j < i:
                        res += Jin[idx] * self.lifter.heisen_bond(i, j)
                        idx += 1
            return res

        H = sps.dok_matrix((dim, dim))
        for b in self.model_data.bonds:
            for spini, i in enumerate(self.full2half[b.i]):
                for spinj, j in enumerate(self.full2half[b.j]):
                    H += b.J[
                        spini * len(self.full2half[b.j]) + spinj
                    ] * self.lifter.heisen_bond(i, j)

        for idx, s in enumerate(self.model_data.sites):
            H += onsite_term(s.Jin, idx, s.h)
        return H

    def l_operator(n):
        lifter = hamiltonian.SpinLifter([2] * n)
        S2 = sum(
            lifter.Sa(a, i) * lifter.Sa(a, j)
            for a in range(3)
            for i in range(n)
            for j in range(n)
        ).real.todense()

        w, v = np.linalg.eigh(S2)
        assert np.allclose(v @ np.diag(w) @ v.T, S2)

        return v @ np.diag(np.sqrt(w + 0.25) - 0.5) @ v.T

    def j_operators(self):
        dim = 2**self.N
        l_opers = {}
        for nspinhalf in set(s.nspinhalfs for s in self.model_data.sites):
            l_opers[nspinhalf] = hamiltonian.l_operator(nspinhalf)

        ops = {}

        for name, q in {
            "J": 0,
            "Jq1": 2 * np.pi / self.model_data.Lx,
            "Jq2": 4 * np.pi / self.model_data.Lx,
        }.items():
            op = sps.dok_matrix((dim, dim), dtype=np.complex)
            pos = 0
            for i, s in enumerate(self.model_data.sites):
                x = (i // self.model_data.uc_site_count) % self.model_data.Lx
                op += np.exp(1j * x * q) * sps.kron(
                    sps.kron(sps.identity(2**pos), l_opers[s.nspinhalfs]),
                    sps.identity(2 ** (self.N - pos - s.nspinhalfs)),
                )
                pos += s.nspinhalfs
            ops[name] = op / len(self.model_data.sites)

        if all(len(s) == 3 for s in self.full2half):
            ops["JDim"] = self.lifter.heisen_bond(0, 1) + 0.75 * Id

        return ops

    def observables(self, params, Ts, E, psi):
        N = sum(s.nspinhalfs for s in self.model_data.sites)

        ens = hamiltonian.Ensemble(E, psi, Ts)

        def mag_obs(prefix, M):
            M2 = M @ M
            M4 = M2 @ M2

            obs = {}
            obs[prefix + "Mag"] = ens.mean(M)
            obs[prefix + "Mag2"] = ens.mean(M2)
            obs[prefix + "Mag4"] = ens.mean(M4)

            obs[prefix + "BinderRatio"] = (
                obs[prefix + "Mag2"] ** 2 / obs[prefix + "Mag4"]
            )

            obs[prefix + "MagChi"] = ens.chi(M, self.N)

            print('observable prefix "{}"'.format(prefix))

            return obs

        def j_obs():
            ops = self.j_operators()
            J = ops["J"]
            Jq1 = ops["Jq1"]
            Jq2 = ops["Jq2"]

            obs = {}
            obs["J"] = ens.mean(J)
            obs["J2"] = ens.mean(J @ J)
            obs["JVar"] = ens.mean(J @ J) - obs["J"] ** 2

            obs["JStruc1"] = ens.mean(
                np.real(Jq1) @ np.real(Jq1) + np.imag(Jq1) @ np.imag(Jq1)
            )
            obs["JStruc2"] = ens.mean(
                np.real(Jq2) @ np.real(Jq2) + np.imag(Jq2) @ np.imag(Jq2)
            )
            r = obs["JStruc1"] / obs["JStruc2"]
            obs["JCorrLen"] = (
                self.model_data.Lx
                / 2
                / np.pi
                * np.sqrt(np.maximum((r - 1) / (4 - r), 0))
            )

            if "JDim" in ops.keys():
                obs["JDim"] = ens.mean(ops["JDim"])

            return obs

        obs = {}
        obs["Energy"] = ens.diag_mean(E) / self.N
        obs["SpecificHeat"] = (
            Ts ** (-2) * (ens.diag_mean(E**2) - ens.diag_mean(E) ** 2) / self.N
        )

        mag_ops = {
            "mag": ("", self.signed_mag(1, 1, 1)),
            "sxmag": ("StagX", self.signed_mag(-1, 1, 1)),
            "symag": ("StagY", self.signed_mag(1, -1, 1)),
            "sxsymag": ("StagXStagY", self.signed_mag(-1, -1, 1)),
            "sxsucmag": ("StagXStagUC", self.signed_mag(-1, 1, -1)),
        }

        for name, op in mag_ops.items():
            if name in params["measure"]:
                obs.update(mag_obs(op[0], op[1]))

        if "chirality" in params["measure"]:
            chirality_ops = self.chirality_operators()
            for name, corrs in chirality_ops.items():
                obs[name] = np.real(np.array([ens.mean(corr) for corr in corrs]).T)

            obs["NematicityAltCorr"] = (
                9 / 16 * obs["NematicityDiagCorr"] + obs["NematicityOffCorr"]
            )

        if "j" in params["measure"]:
            obs.update(j_obs())

        for obsname in ["TauZ", "TauY"] + [
            "Nematicity{}Corr".format(n) for n in ["", "Diag", "Off", "Cross"]
        ]:
            if obsname in obs.keys():
                obs[obsname.strip("Corr") + "Struc"] = np.sum(obs[obsname], axis=1)

        return obs
