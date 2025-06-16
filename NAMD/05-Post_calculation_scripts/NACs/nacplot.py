#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

nactxt = np.loadtxt('NATXT', dtype=float)
eigtxt = np.loadtxt('EIGTXT', dtype=float)

nbasis = eigtxt.shape[1]
nactxt = nactxt.reshape((-1, nbasis, nbasis))

hbar = 0.6582119281559802
potim = 1.0

# NAC in unit of meV
a_nac = np.average(np.abs(nactxt), axis=0) * hbar / (2 * potim) * 1000
a_eig = np.average(eigtxt, axis=0) 

# plot part
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches((4, 3))
plt.subplots_adjust(left=0.10, right=0.95,
                    bottom=0.10, top=0.95,
                    wspace=0.20, hspace=0.20)
ax.set_aspect(1.0)

# bd = ax.get_position().bounds
# ax_cbar = fig.add_axes([bd[0] + bd[2], bd[1], 0.05, bd[-2]])

# x, y = np.meshgrid(a_eig, a_eig)
# img = ax.scatter(x, y, c=a_nac,
#                  s=3,
#                  cmap='Reds',
#                  # lw=0.5, edgecolors='k',
#                  marker='s')

img = ax.imshow(a_nac, origin='lower', interpolation='none',
                # cmap='jet'
                )

cbar = plt.colorbar(img, pad=0.01)
cbar.ax.tick_params(axis='y', which='both', labelsize='medium')
cbar.ax.text(1.15, 1.0, 'meV',
            fontsize='medium',
            ha='left',
            va='top', transform=cbar.ax.transAxes)

# ax.axis([a_eig.min(), a_eig.max(), a_eig.min(), a_eig.max()])
# ax.set_xlim(a_eig.min(), a_eig.max())
# ax.set_ylim(a_eig.min(), a_eig.max())

# ax.set_xlabel('Energy (eV)', labelpad=5)
# ax.set_ylabel('Energy (eV)', labelpad=5)

ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout(pad=0.2)
plt.savefig('nac.png', dpi=300)
