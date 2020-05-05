import numpy as np
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import lattice
import argparse
import ed_hamiltonian
from loadleveller import jobfile
import os

parser = argparse.ArgumentParser(description='ED code for the frust project. Automatically calculates the results for a chosen loadleveller jobfile. If there is more than one task in the jobfile, those tasks must only vary by temperature.')

parser.add_argument('jobfile', type=str, help='Jobscript describing the model.')
parser.add_argument('-o', '--outfile', type=str, help='.npz file to save results to')

args = parser.parse_args()

def extract_jobdata(job):
    tasknames = []
    Ts = []
    for task, params in job.tasks.items():
        if params['ed_compare']:
            tasknames.append(task)
            Ts.append(params['T'])
    if len(tasknames) == 0:
        raise Exception('jobfile did not contain any tasks with ed_compare=True')

    return tasknames[0], np.array(Ts)

jobdir = os.path.dirname(args.jobfile)
if jobdir != '':
    os.chdir(jobdir)

job = jobfile.JobFile('./'+os.path.basename(args.jobfile))

taskname, Ts = extract_jobdata(job)

num_states = job.tasks[taskname].get('ed_num_states', 0)
lat = lattice.load(job, taskname)

N = sum(s.nspinhalfs for s in lat.sites)
print('dimension = 2**{} = {}'.format(N, 2**N))
H, obs_ops = ed_hamiltonian.construct(lat)

print('H constructed.')
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

na = np.newaxis
Egap = E[E>E.min()+1e-8].min()-E.min()
print('{} energies: {:.3g}..{:.3g}, gap={:.3g}'.format(len(E), E.min(), E.max(), Egap))

def calc_observables(E, psi):
    Enorm = E-E.min()
    ρ = np.exp(-Enorm[na,:]/Ts[:,na])
    Z = np.sum(ρ, axis=1)
    ρ /= Z[:,na]

    def mean(A):
        nAn = np.zeros(psi.shape[1])
        for i in range(psi.shape[1]):
            nAn[i] = np.real(psi[:,i].conj()@ A @psi[:,i])

        return np.sum(nAn[na,:]*ρ, axis=1)
        

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
        M2 = np.dot(M,M)
        M4 = np.dot(M2,M2)

        obs = {}
        obs[prefix+'Mag'] = mean(M)
        obs[prefix+'Mag2'] = mean(M2)
        obs[prefix+'Mag4'] = mean(M4)

        if num_states == 0:
            obs[prefix+'MagChi'] = chi(M, N)
                
        #obs[prefix+'BinderRatio'] = obs[prefix+'Mag2']**2/obs[prefix+'Mag4']
        
        #Lx = job.tasks[taskname]['Lx']
        #obs[prefix+'StrucFac2'] = structure_fac([2*np.pi/Lx, 0], spin, positions, sset)
        #obs[prefix+'CorrLenF'] = np.sqrt(np.maximum(len(sset)*obs[prefix+'Mag2']/obs[prefix+'StrucFac2'],1)-1)/2/np.pi

        print('observable prefix "{}"'.format(prefix))

        return obs
    def j_obs():
        J = obs_ops['J']

        obs = {}
        obs['J'] = mean(J)
        obs['J2'] = mean(J@J)
        obs['JVar'] = mean(J@J)-obs['J']**2

        return obs


    obs = {}
    obs['Energy'] = np.sum(E[na,:]*ρ,axis=1)/N
    obs['SpecificHeat'] = Ts**(-2)*(np.sum(E[na,:]**2*ρ,axis=1) - np.sum(E[na,:]*ρ,axis=1)**2)/N
    obs.update(mag_obs('', obs_ops['M']))
    obs.update(mag_obs('Stag', obs_ops['sM']))
    if 'J' in obs_ops.keys():
        obs.update(j_obs())

    return obs

obs = calc_observables(E, psi)

if args.outfile:
    np.savez_compressed(args.outfile, **obs)
    print('Saved to "{}"'.format(args.outfile))
else:
    for k, v in obs.items():
        print('{}: {}'.format(k,v))
        


