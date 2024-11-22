import numpy as np
import scipy.sparse as sps
import scipy.linalg as spl


def spinop(component, spin_dim):
    S = (spin_dim - 1) * 0.5
    m = S - np.arange(spin_dim)
    if component < 2:
        splus = sps.diags([np.sqrt(S * (S + 1) - m[1:] * (m[1:] + 1))], offsets=[1])
        if component == 0:
            return (splus + splus.H) * 0.5
        return -(splus - splus.H) * 0.5j
    if component == 2:
        return sps.diags([m], [0])


σz2 = np.array([[1, 0], [0, -1]])
σy2 = np.array([[0, -1j], [1j, 0]])
σx2 = np.array([[0, 1], [1, 0]])


class Lifter:
    def __init__(self, site_dims):
        self.site_dims = site_dims
        self.cum_site_dims = np.append([1], np.cumprod(site_dims))

    def lift_op(self, pos, op):
        assert all(d == self.site_dims[pos] for d in op.shape)
        return sps.kron(
            sps.kron(sps.identity(self.cum_site_dims[pos]), op),
            sps.identity(self.cum_site_dims[-1] / self.cum_site_dims[pos + 1]),
        )


class SpinLifter(Lifter):
    def __init__(self, site_dims):
        Lifter.__init__(self, site_dims)

    def Sa(self, a, pos):
        return self.lift_op(pos, spinop(a, self.site_dims[pos]))

    def Sz(self, pos):
        return self.Sa(2, pos)

    def heisen_bond(self, i, j):
        return sum(self.Sa(a, i) @ self.Sa(a, j) for a in range(3)).real


class Ensemble:
    def __init__(self, E, psi, Ts):
        self.Ts = Ts
        Enorm = E - E.min()
        self.ρ = np.exp(-Enorm[None, :] / Ts[:, None])
        Z = np.sum(self.ρ, axis=1)
        self.ρ /= Z[:, None]
        self.psi = psi
        self.E = E

    def mean(self, A):
        nAn = np.zeros(self.psi.shape[1], dtype=np.complex128)
        for i in range(self.psi.shape[1]):
            nAn[i] = self.psi[:, i].conj() @ A @ self.psi[:, i]

        res = np.sum(nAn[None, :] * self.ρ, axis=1)
        if np.allclose(np.imag(res), 0):
            res = np.real(res)
        return res

    def diag_mean(self, A):
        return np.sum(A[None, :] * self.ρ, axis=1)

    def integcorrfunc(self, A):  # ∫^β_0 dτ ⟨A(τ)A⟩ integrated correlation function
        result = np.zeros_like(self.Ts)

        Anm = np.dot(self.psi.conj().T, np.dot(np.array(A.todense()), self.psi))
        for i, t in enumerate(self.Ts):
            deltaE = self.E[:, None] - self.E[None, :]

            deltaE[np.isclose(deltaE, 0, atol=1e-6)] = t * 2

            result[i] = np.sum(2 * self.ρ[i, None, :] / deltaE * np.abs(Anm) ** 2)
        return result

    def chi(self, M, N):
        return self.integcorrfunc(M) * N
