import numpy as np
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import lattice
import argparse
import json
from collections import defaultdict
import ed_hamiltonian
from loadleveller import jobfile
import os

parser = argparse.ArgumentParser(description='ED code for the frust project. Automatically calculates the results for a chosen loadleveller jobfile. If there is more than one task in the jobfile, those tasks must only vary by temperature.')

parser.add_argument('jobfile', type=str, help='Jobscript describing the model.')
parser.add_argument('-o', '--outfile', type=str, help='.json file to save results to')

args = parser.parse_args()

def extract_jobdata(job):
    summarized_tasks = []
    Ts = []
    for task, params in job.tasks.items():

        T = params.pop('T')
        found = False
        for previous in summarized_tasks:
            if params == previous[1]:
                previous[2].append(T)
                found = True
        if not found:
            summarized_tasks.append((task, params, [T]))
        

    return summarized_tasks

def solve_task(H, num_states=0):
    if num_states == 0:
        H = np.array(H.todense())
        E, psi = np.linalg.eigh(H)
        print('Full diagonalization done.')
    else:
        E, psi = spsl.eigsh(H, num_states, which='SA')
        Emax = E.max()
        psi = psi[:,E<E.max()-1e-3]
        E = E[E<E.max()-1e-3]
        print('Sparse diagonalization done ({}/{} states).'.format(len(E),num_states))

    Egap = E[E>E.min()+1e-8].min()-E.min()
    print('{} energies: {:.3g}..{:.3g}, gap={:.3g}'.format(len(E), E.min(), E.max(), Egap))

    return E, psi

def calc_observables(params, Ts, E, psi, obs_ops):
    N = sum(s.nspinhalfs for s in lat.sites)
    na = np.newaxis
    Enorm = E - E.min()
    ρ = np.exp(-Enorm[na,:]/Ts[:,na])
    Z = np.sum(ρ, axis=1)
    ρ /= Z[:,na]

    def mean(A):
        nAn = np.zeros(psi.shape[1], dtype=np.complex128)
        for i in range(psi.shape[1]):
            nAn[i] = psi[:,i].conj()@ A @psi[:,i]

        res = np.sum(nAn[na,:]*ρ, axis=1)
        if np.allclose(np.imag(res), 0):
            res = np.real(res)
        return res

    def integcorrfunc(A): # ∫^β_0 dτ ⟨A(τ)A⟩ integrated correlation function
        result = np.zeros_like(Ts)
        
        Anm = np.dot(psi.conj().T,np.dot(np.array(A.todense()),psi))
        for i, t in enumerate(Ts):
            deltaE = E[:,na]-E[na,:]
            
            deltaE[np.isclose(deltaE,0,atol=1e-6)] = t*2
            
            result[i] = np.sum(2*ρ[i,na,:]/deltaE*np.abs(Anm)**2)
        return result

    def chi(M, N):
        return integcorrfunc(M)*N

    # def structure_fac(q, spin, positions, sset):
    #     #struc = sum(np.cos(np.dot(q, positions[i]-positions[0]))*sign(0,i)*np.dot(spin(0),spin(i)) for i in range(N))
    #     struc = sum(np.exp(1j*np.dot(q, positions[i]))*spin(i) for i in sset)

    #     return mean(np.dot(struc, struc.conj().T))/len(sset)

    def mag_obs(prefix, M):
        M2 = M @ M
        M4 = M2 @ M2

        obs = {}
        obs[prefix+'Mag'] = mean(M)
        obs[prefix+'Mag2'] = mean(M2)
        obs[prefix+'Mag4'] = mean(M4)
        
        obs[prefix+'BinderRatio'] = obs[prefix+'Mag2']**2/obs[prefix+'Mag4']

        obs[prefix+'MagChi'] = chi(M, N)
                
        #Lx = job.tasks[taskname]['Lx']
        #obs[prefix+'StrucFac2'] = structure_fac([2*np.pi/Lx, 0], spin, positions, sset)
        #obs[prefix+'CorrLenF'] = np.sqrt(np.maximum(len(sset)*obs[prefix+'Mag2']/obs[prefix+'StrucFac2'],1)-1)/2/np.pi

        print('observable prefix "{}"'.format(prefix))

        return obs
    def j_obs():
        J = obs_ops['J']
        Jq1 = obs_ops['Jq1']
        Jq2 = obs_ops['Jq2']

        obs = {}
        obs['J'] = mean(J)
        obs['J2'] = mean(J@J)
        obs['JVar'] = mean(J@J)-obs['J']**2

        obs['JStruc1'] = mean(np.real(Jq1)@np.real(Jq1)+np.imag(Jq1)@np.imag(Jq1))
        obs['JStruc2'] = mean(np.real(Jq2)@np.real(Jq2)+np.imag(Jq2)@np.imag(Jq2))
        r = obs['JStruc1']/obs['JStruc2']
        obs['JCorrLen'] = lat.Lx/2/np.pi*np.sqrt(np.maximum((r-1)/(4-r),0))

        return obs


    obs = {}
    obs['Z'] = Z
    obs['Energy'] = np.sum(E[na,:]*ρ,axis=1)/N
    obs['SpecificHeat'] = Ts**(-2)*(np.sum(E[na,:]**2*ρ,axis=1) - np.sum(E[na,:]*ρ,axis=1)**2)/N

    if 'mag' in params['measure']:
        obs.update(mag_obs('', obs_ops['M']))
    if 'sxmag' in params['measure']:
        obs.update(mag_obs('StagX', obs_ops['sxM']))
    if 'symag' in params['measure']:
        obs.update(mag_obs('StagY', obs_ops['syM']))
    if 'sxsymag' in params['measure']:
        obs.update(mag_obs('StagXStagY', obs_ops['sxsyM']))
    if 'sxsucmag' in params['measure']:
        obs.update(mag_obs('StagXStagUC', obs_ops['sxsucM']))

    if 'chirality' in params['measure']:
        obs['TauZ'] = np.real(np.array([mean(taucorr) for taucorr in obs_ops['chirality_tauz']]).T)
        obs['TauY'] = np.real(np.array([mean(taucorr) for taucorr in obs_ops['chirality_tauy']]).T)
        obs['NematicityCorr'] = np.real(np.array([mean(corr) for corr in obs_ops['NematicityCorr']]).T)
        obs['NematicityDiagCorr'] = np.real(np.array([mean(corr) for corr in obs_ops['NematicityDiagCorr']]).T)
        obs['NematicityOffCorr'] = np.real(np.array([mean(corr) for corr in obs_ops['NematicityOffCorr']]).T)
        #obs['NematicityCrossCorr'] = np.array([mean(corr) for corr in obs_ops['NematicityCrossCorr']]).T

        obs['NematicityAltCorr'] = 9/16*obs['NematicityDiagCorr'] + obs['NematicityOffCorr']
        
    if 'j' in params['measure']:
        obs.update(j_obs())
        if 'JDim' in obs_ops.keys():
            obs['JDim'] = mean(obs_ops['JDim'])
    for obsname in ['TauZ', 'TauY'] + ['Nematicity{}Corr'.format(n) for n in ['', 'Diag', 'Off', 'Cross']]:
        if obsname in obs.keys():
            obs[obsname.strip('Corr') + 'Struc'] = np.sum(obs[obsname],axis=1)

    return obs

def to_list(a):
    if type(a) == np.ndarray:
        return a.tolist()
    return a

jobdir = os.path.dirname(args.jobfile)
if jobdir != '':
    os.chdir(jobdir)

job = jobfile.JobFile('./'+os.path.basename(args.jobfile))

output = []

summarized_tasks = extract_jobdata(job)
for i, (taskname, params, Ts) in enumerate(summarized_tasks):
    print('task {}/{}...'.format(i+1, len(summarized_tasks)))    
    lat = lattice.load(job, taskname, force_overwrite=True)

    H, obs_ops = ed_hamiltonian.construct(lat)
    print('dimension = 2**{:d} = {}'.format(int(np.log(H.shape[0])/np.log(2)), H.shape[0]))
    
    E, psi = solve_task(H, params.get('ed_num_states', 0))

    obs = calc_observables(params, np.array(Ts), E, psi, obs_ops)
    print('observables calculated')

    for i, T in enumerate(Ts):
        task_observables = {obsname: to_list(obsvalue[i]) for obsname, obsvalue in obs.items()}
        output.append({'parameters': {'T': T, **params}, 'results': task_observables})

if args.outfile:
    with open(args.outfile, 'w') as outfile:
        json.dump(output, outfile, indent=1)
else:
    print(json.dumps(output, indent=4))
    

