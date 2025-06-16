#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, subprocess
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

# Load data
nactxt = np.loadtxt('NATXT', dtype=float)
eigtxt = np.loadtxt('EIGTXT', dtype=float)

nbasis = eigtxt.shape[1]
nactxt = nactxt.reshape((-1, nbasis, nbasis))

hbar = 0.6582119281559802
potim = 1.0

# NAC in unit of meV
a_nac = np.average(np.abs(nactxt), axis=0) * hbar / (2 * potim) * 1000
a_eig = np.average(eigtxt, axis=0)

# Save `a_nac` to .dat files
np.savetxt('nac_data.dat', a_nac, fmt='%.6f', delimiter='\t')

# Plotting and saving the image
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches((4, 3))
plt.subplots_adjust(left=0.10, right=0.95, bottom=0.10, top=0.95, wspace=0.20, hspace=0.20)
ax.set_aspect(1.0)

# Plot the heatmap
img = ax.imshow(a_nac, origin='lower', interpolation='none', cmap='viridis')
cbar = plt.colorbar(img, pad=0.01)
cbar.ax.tick_params(axis='y', which='both', labelsize='medium')
cbar.ax.text(1.15, 1.0, 'meV', fontsize='medium', ha='left', va='top', transform=cbar.ax.transAxes)

# Add annotations
for i in range(a_nac.shape[0]):
    for j in range(a_nac.shape[1]):
        ax.text(j, i, f"{a_nac[i, j]:.1f}", ha="center", va="center", color="white", fontsize=6)

ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout(pad=0.2)
plt.savefig('nac.png', dpi=300)

# Open saved image
if os.name == "nt":  # Windows
    subprocess.call(["explorer", 'nac.png'])
elif "microsoft" in os.uname().release.lower():
    subprocess.call(["explorer.exe", 'nac.png'])
else:  # Linux/macOS
    subprocess.call(["xdg-open", 'nac.png']) 
