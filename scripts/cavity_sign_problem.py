import ed.downfolded_coupling as coupling
import matplotlib.pyplot as plt
import numpy as np

gs = np.linspace(0.01, 1, 80)
omegas = np.linspace(0, 1.5, 81)
nphot = 10
data = np.zeros([len(gs), len(omegas)])

for i, g in enumerate(gs):
    for j, omega in enumerate(omegas):
        mat = coupling.matrix(nphot, omega, g)
        data[i, j] = np.all(mat >= -1e-9)

momegas, mgs = np.meshgrid(omegas, gs)
plt.pcolor(mgs, momegas, data)
plt.colorbar(label="sign")
plt.ylabel("$\Omega/U$")
plt.xlabel("$g$")
plt.tight_layout()
plt.show()
