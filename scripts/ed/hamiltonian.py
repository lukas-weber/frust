import numpy as np
import scipy.sparse as sps
import scipy.linalg as spl

σz2 = np.array([[1,0],[0,-1]])
σy2 = np.array([[0,-1j],[1j,0]])
σx2 = np.array([[0,1],[1,0]])

def lift_op(pos, op, N):
    return sps.kron(sps.kron(sps.identity(2**(pos)), op), sps.identity(2**(N-pos-1)))

def Sx(pos, N):
    return lift_op(pos, σx2/2, N)
def Sy(pos, N):
    return lift_op(pos, σy2/2, N)
def Sz(pos, N):
    return lift_op(pos, σz2/2, N)
def heisen_bond(i, j, N):
    return (Sx(i,N)@Sx(j,N) + Sy(i,N)@Sy(j,N) + Sz(i,N)@Sz(j,N)).real

def l_operator(n):
    S2 = sum(Sx(i,n)*Sx(j,n)+Sy(i,n)*Sy(j,n)+Sz(i,n)*Sz(j,n) for i in range(n) for j in range(n)).real.todense()

    w, v = np.linalg.eigh(S2)
    assert np.allclose(v@np.diag(w)@v.T,S2)

    return v@np.diag(np.sqrt(w+0.25)-0.5)@v.T
    
class Ensemble:
    def __init__(self, E, psi, Ts):
        self.Ts = Ts
        Enorm = E - E.min()
        self.ρ = np.exp(-Enorm[None,:]/Ts[:,None])
        Z = np.sum(self.ρ, axis=1)
        self.ρ /= Z[:,None]
        self.psi = psi
        self.E = E

    def mean(self, A):
        nAn = np.zeros(self.psi.shape[1], dtype=np.complex128)
        for i in range(self.psi.shape[1]):
            nAn[i] = self.psi[:,i].conj()@ A @self.psi[:,i]

        res = np.sum(nAn[None,:]*self.ρ, axis=1)
        if np.allclose(np.imag(res), 0):
            res = np.real(res)
        return res

    def diag_mean(self, A):
        return np.sum(A[None,:]*self.ρ,axis=1)

    def integcorrfunc(self, A): # ∫^β_0 dτ ⟨A(τ)A⟩ integrated correlation function
        result = np.zeros_like(self.Ts)
        
        Anm = np.dot(self.psi.conj().T,np.dot(np.array(A.todense()),self.psi))
        for i, t in enumerate(self.Ts):
            deltaE = self.E[:,None]-self.E[None,:]
            
            deltaE[np.isclose(deltaE,0,atol=1e-6)] = t*2
            
            result[i] = np.sum(2*self.ρ[i,None,:]/deltaE*np.abs(Anm)**2)
        return result

    def chi(self, M, N):
        return self.integcorrfunc(M)*N
