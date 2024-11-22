import scipy.sparse as sps


class Magnet:
    mag_sign_cfgs = {
        "mag": ("", (1, 1, 1)),
        "sucmag": ("StagUC", (1, 1, -1)),
        "sxmag": ("StagX", (-1, 1, 1)),
        "symag": ("StagY", (1, -1, 1)),
        "sxsymag": ("StagXStagY", (-1, -1, 1)),
        "sxsucmag": ("StagXStagUC", (-1, 1, -1)),
    }

    def signed_magnetization(self, spin_op, sx, sy, suc):
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

                    M += sign * spin_op(i)
        return M / self.N

    def observables_magnetization(self, ens, M):
        M2 = M @ M
        M4 = M2 @ M2

        obs = {}
        obs["Mag"] = ens.mean(M)
        obs["Mag2"] = ens.mean(M2)
        obs["Mag4"] = ens.mean(M4)

        obs["BinderRatio"] = obs["Mag2"] ** 2 / obs["Mag4"]

        obs["MagChi"] = ens.chi(M, self.N)
        return obs
