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

    return v@np.diag(np.sqrt(w+0.25)-0.5)@v.T
    

def chirality(Nfull, obs_ops):
    j = np.exp(2j*np.pi/3)

    uL = 1/3**0.5*np.matrix([0,1,j,0,j*j,0,0,0])
    uR = 1/3**0.5*np.matrix([0,1,j*j,0,j,0,0,0])
    dL = 1/3**0.5*np.matrix([0,0,0,j*j,0,j,1,0])
    dR = 1/3**0.5*np.matrix([0,0,0,j,0,j*j,1,0])

    us = 1/2**0.5 * np.matrix([0,0,1,0,-1,0,0,0])
    ut = 1/6**0.5 * np.matrix([0,-2,1,0,1,0,0,0])
    ds = 1/2**0.5 * np.matrix([0,0,0,-1,0,1,0,0])
    dt = 1/6**0.5 * np.matrix([0,0,0,1,0,1,-2,0])

    tauz = uL.H@uL + dL.H@dL - uR.H@uR - dR.H@dR
    #tauz = -1j*(uL.H@uR + dL.H@dR - uR.H@uL - dR.H@dL)
    #tauz = us.H@ut + ds.H@dt + ut.H@us + dt.H@ds
    #tauz = us.H@ut + ds.H@dt - ut.H@us - dt.H@ds
    tauy = -1j*(uL.H@uR + dL.H@dR - uR.H@uL - dR.H@dL)

    
    def lift_tri(i, op):
        return sps.kron(sps.kron(sps.identity(8**i), op), sps.identity(8**(Nfull-i-1)))

    nem_diag = ut.H@ut-us.H@us + dt.H@dt - ds.H@ds
    nem_off = 3/4*(ut.H@us  + us.H@ut + dt.H@ds + ds.H@dt)

    
    nem_diag_corr = []
    tauzcorr = []
    tauycorr = []

    nem_off_corr = []
    nem_cross_corr = []
    for i in range(Nfull):
        tauzcorr.append(lift_tri(0,tauz)@lift_tri(i,tauz))
        tauycorr.append(lift_tri(0,tauy)@lift_tri(i,tauy))
        nem_diag_corr.append(lift_tri(0,nem_diag)@lift_tri(i,nem_diag))
        nem_off_corr.append(lift_tri(0,nem_off)@lift_tri(i,nem_off.H))
        nem_cross_corr.append(1j*(lift_tri(0,nem_diag)@lift_tri(i,nem_off.H)-lift_tri(0,nem_off)@lift_tri(i,nem_diag.H)))
    #meantauz = -sum(lift_tauz(i)@lift_tauz(j) for i in range(Nfull) for j in range(Nfull) if i != j)/Nfull**2
    #
    obs_ops['chirality_tauz'] = tauzcorr
    obs_ops['chirality_tauy'] = tauycorr
    obs_ops['NematicityDiagCorr'] = nem_diag_corr
    obs_ops['NematicityOffCorr'] = nem_off_corr
    obs_ops['NematicityCrossCorr'] = nem_cross_corr


def construct(lat):
    Nfull = len(lat.sites)
    half2full = sum([site.nspinhalfs*[i] for i, site in enumerate(lat.sites)], [])
    full2half = [[i for i, x in enumerate(half2full) if x == idx] for idx in range(Nfull)]

    N = len(half2full)

    Id = sps.identity(2**N)

    def H_heisen_bond(i, j, N=N):
        return (Sx(i,N)@Sx(j,N) + Sy(i,N)@Sy(j,N) + Sz(i,N)@Sz(j,N)).real

    def onsite_term(Jin, site,h):
        res = sps.dok_matrix((2**(N), 2**(N)))
        idx = 0
        for i in full2half[site]:
            res += Sz(i,N) * h
            for j in full2half[site]:
                if j < i:
                    res += Jin[idx]*H_heisen_bond(i, j)
                    idx += 1
        return res

    H = sps.dok_matrix((2**N, 2**N))
    for b in lat.bonds:
        for spini, i in enumerate(full2half[b.i]):
            for spinj, j in enumerate(full2half[b.j]):
                H += b.J[spini*len(full2half[b.j])+spinj] * H_heisen_bond(i, j)

    for idx, s in enumerate(lat.sites):
        H += onsite_term(s.Jin, idx,s.h)

    l_opers = {}
    for nspinhalf in set(s.nspinhalfs for s in lat.sites):
        l_opers[nspinhalf] = l_operator(nspinhalf)
    
    obs_ops = {}
    obs_ops['M'] = sum(Sz(i, N) for i in range(N))/N

    w = np.exp(1j*2*np.pi/3)


    for name, q in {'J':0, 'Jq1':2*np.pi/lat.Lx, 'Jq2':4*np.pi/lat.Lx}.items():
        op = sps.dok_matrix((2**N,2**N),dtype=np.complex)
        pos = 0
        for i, s in enumerate(lat.sites):
            x = (i//lat.uc_spin_count)%lat.Lx
            op += np.exp(1j*x*q)*sps.kron(sps.kron(sps.identity(2**pos), l_opers[s.nspinhalfs]),sps.identity(2**(N-pos-s.nspinhalfs)))
            pos += s.nspinhalfs
        obs_ops[name] = op/len(lat.sites) 

    def signed_mag(sx, sy, suc):
        M = sps.dok_matrix((2**N,2**N))
        for x in range(lat.Lx):
            for y in range(lat.Ly):
                for uc in range(lat.uc_spin_count):
                    i = lat.uc_spin_count*(lat.Lx*y+x)+uc
                    sign = sx**x * sy**y * (lat.uc_signs[uc] if suc < 0 else 1)

                    for idx in full2half[i]:
                        M += sign*Sz(idx,N)
        return M/N
    obs_ops['sxM'] = signed_mag(-1,1,1)
    obs_ops['syM'] = signed_mag(1,-1,1)
    obs_ops['sxsyM'] = signed_mag(-1,-1,1)
    obs_ops['sxsucM'] = signed_mag(-1,1,-1)

    if all(len(s)==3 for s in full2half):
        chirality(Nfull, obs_ops)
    
        C = H_heisen_bond(0,1) + w * H_heisen_bond(1,2) + w**2 * H_heisen_bond(2,0)

        obs_ops['NematicityCorr'] = [C @ (H_heisen_bond(3*i+0,3*i+1) +
                w * H_heisen_bond(3*i+1,3*i+2) + w**2 * H_heisen_bond(3*i+2,3*i+0)).H for i in range(0,Nfull)]
        obs_ops['C'] = C

        #idx=0
        #nem_diff = 9/16*obs_ops['NematicityDiagCorr'][idx] + obs_ops['NematicityOffCorr'][idx] + 3/4*obs_ops['NematicityCrossCorr'][idx] - obs_ops['NematicityCorr'][idx]
        #nem_diff[np.abs(nem_diff) < 1e-10] = 0
        #nem_diff = nem_diff.todense()

        #print(nem_diff[np.nonzero(nem_diff)])


        obs_ops['JDim'] = H_heisen_bond(0,1) + 0.75*Id
        
    return H, obs_ops
