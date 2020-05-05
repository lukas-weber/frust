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

def l_operator(n):
    S2 = sum(Sx(i,n)*Sx(j,n)+Sy(i,n)*Sy(j,n)+Sz(i,n)*Sz(j,n) for i in range(n) for j in range(n)).real.todense()

    w, v = np.linalg.eigh(S2)
    assert np.allclose(v@np.diag(w)@v.T,S2)

    return v@np.diag(np.sqrt(w+0.25)+0.5)@v.T
    

def construct(lat):
    Nmult = len(lat.sites)
    half2mult = sum([site.nspinhalfs*[i] for i, site in enumerate(lat.sites)], [])
    mult2half = [[i for i, x in enumerate(half2mult) if x == idx] for idx in range(Nmult)]

    N = len(half2mult)

    Id = sps.identity(2**N)

    def H_heisen_bond(i, j):
        return Sx(i,N)@Sx(j,N) + Sy(i,N)@Sy(j,N) + Sz(i,N)@Sz(j,N)

    def onsite_term(site):
        res = sps.dok_matrix((2**(N), 2**(N)))
        for i in mult2half[idx]:
            for j in mult2half[idx]:
                if i < j:
                    res += H_heisen_bond(i, j)
        return res

    H = sps.dok_matrix((2**N, 2**N))
    for b in lat.bonds:
        for spini, i in enumerate(mult2half[b.i]):
            for spinj, j in enumerate(mult2half[b.j]):
                H += b.J[spini*len(mult2half[b.j])+spinj] * H_heisen_bond(i, j)

    for idx, s in enumerate(lat.sites):
        H += s.Jin * onsite_term(idx)


    l_opers = {}
    for nspinhalf in set(s.nspinhalfs for s in lat.sites):
        l_opers[nspinhalf] = l_operator(nspinhalf)
    
    obs_ops = {}
    obs_ops['M'] = sum(Sz(i, N) for i in range(N))/N
    obs_ops['sM'] = sum(lat.sites[half2mult[i]].sublattice * Sz(i,N) for i in range(N))/N
    
    return H, obs_ops
