#!/usr/bin/env python3
############################################################
import os, re
import numpy as np
from glob import glob

import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['axes.unicode_minus'] = False

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
############################################################
def WeightFromPro(infile='PROCAR', whichAtom=None, spd=None):
    """
    Contribution of selected atoms to the each KS orbital
    """

    print(infile) 
    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub('[^0-9]', ' ', FileContents[1]).split()]

    if spd:
        Weights = np.asarray([line.split()[1:-1] for line in FileContents
                              if not re.search('[a-zA-Z]', line)], dtype=float)
        Weights = np.sum(Weights[:,spd], axis=1)
    else:
        Weights = np.asarray([line.split()[-1] for line in FileContents
                              if not re.search('[a-zA-Z]', line)], dtype=float)
    
    nspin = Weights.shape[0] / (nkpts * nbands * nions)
    Weights.resize(nspin, nkpts, nbands, nions)

    Energies = np.asarray([line.split()[-4] for line in FileContents
                            if 'occ.' in line], dtype=float)
    Energies.resize(nspin, nkpts, nbands)
    
    if whichAtom is None:
        return Energies, np.sum(Weights, axis=-1)
    else:
        # whichAtom = [xx - 1 for xx in whichAtom]
        return Energies, np.sum(Weights[:,:,:,whichAtom], axis=-1)

def parallel_wht(runDirs, whichAtoms, nproc=None):
    '''
    calculate localization of some designated in parallel.
    '''
    import multiprocessing
    nproc = multiprocessing.cpu_count() if nproc is None else nproc
    pool = multiprocessing.Pool(processes=nproc)

    results = []
    for rd in runDirs:
        res = pool.apply_async(WeightFromPro, (rd + '/PROCAR', whichAtoms, None,))
        results.append(res)

    enr = []
    wht = []
    for ii in range(len(results)):
        tmp_enr, tmp_wht = results[ii].get()
        enr.append(tmp_enr)
        wht.append(tmp_wht)

    return np.array(enr), np.array(wht)

############################################################
# calculate spatial localization
############################################################
nsw     = 2000
dt      = 1.0
nproc   = 32
prefix  = '../run'
runDirs = [prefix + '/{:04d}'.format(ii + 1) for ii in range(nsw)]
# which spin
whichS  = 0
# which k-point
whichK  = 0 
# which atoms
whichA = np.arange(84)      # Indices 0 to 83 for the first 84 atoms (MSN)
whichB  = np.arange(84, 124)  # Indices 84 to 123 for the remaining atoms (NTO)
Alabel  = 'MoSi$_2$N$_4$'
Blabel  = 'Ni$_3$TeO$_6$'

if os.path.isfile('all_wht.npy'):
    Wht = np.load('all_wht.npy')
    Enr = np.load('all_en.npy')
else:
    # for gamma point version, no-spin
    Enr, Wht = parallel_wht(runDirs, whichA, nproc=nproc)
    Enr = Enr[:, whichS,whichK, :]
    Wht = Wht[:, whichS,whichK, :]

    # Enr, Wht1 = parallel_wht(runDirs, whichA, nproc=nproc)
    # Enr, Wht2 = parallel_wht(runDirs, whichB, nproc=nproc)
    # Enr = Enr[:, whichS,whichK, :]
    # Wht1 = Wht1[:, whichS,whichK, :]
    # Wht2 = Wht2[:, whichS,whichK, :]
    # Wht = Wht1 / (Wht1 + Wht2)

    np.save('all_wht.npy', Wht)
    np.save('all_en.npy', Enr)

############################################################
fig = plt.figure()
fig.set_size_inches(4.8, 3.0)

########################################
ax      = plt.subplot()
nband   = Enr.shape[1]
T, dump = np.mgrid[0:nsw:dt, 0:nband]
sFac    = 8
'''
############################################################
# METHOD 1.
############################################################
# use scatter to plot the band 
# ax.scatter( T, Enr, s=Wht / Wht.max() * sFac, color='red', lw=0.0, zorder=1)
# for ib in range(nband):
#     ax.plot(T[:,ib], Enr[:,ib], lw=0.5, color='k', alpha=0.5)

############################################################
# METHOD 2.
############################################################
# use colored scatter to plot the band 
img = ax.scatter(T, Enr, s=1.0, c=Wht, lw=0.0, zorder=1,
                 vmin=Wht.min(),
                 vmax=Wht.max(),
                 cmap='jet_r')
# for ib in range(nband):
#     ax.plot(T[:,ib], Enr[:,ib], lw=0.5, color='k', alpha=0.5)

divider = make_axes_locatable(ax)
ax_cbar = divider.append_axes('right', size='5%', pad=0.02)
cbar = plt.colorbar(img, cax=ax_cbar,
                    orientation='vertical')
cbar.set_ticks([Wht.min(), Wht.max()])
cbar.set_ticklabels([Alabel, Blabel])

############################################################
# METHOD 3.
############################################################
# # use color strip to plot the band
#
# LW    = 1.0
# DELTA = 0.3
# norm  = mpl.colors.Normalize(vmin=Wht.min(),
#                              vmax=Wht.max())
# # create a ScalarMappable and initialize a data structure
# s_m   = mpl.cm.ScalarMappable(cmap='jet_r', norm=norm)
# s_m.set_array([Wht])
#
# x     = np.arange(0, nsw, dt)
# # for iband in range(nband):
# for iband in range(100, 110):
#     print('Processing band: {:4d}...'.format(iband))
#     y = Enr[:,iband]
#     z = Wht[:,iband]
#
#     ax.plot(x, y,
#             lw=LW + 2 * DELTA,
#             color='gray', zorder=1)
#
#     points = np.array([x, y]).T.reshape(-1, 1, 2)
#     segments = np.concatenate([points[:-1], points[1:]], axis=1)
#     lc = LineCollection(segments,
#                         colors=[s_m.to_rgba(ww) for ww in (z[1:] + z[:-1])/2.]
#                         )
#     # lc.set_array((z[1:] + z[:-1]) / 2)
#     lc.set_linewidth(LW)
#     ax.add_collection(lc)
#
#     divider = make_axes_locatable(ax)
#     ax_cbar = divider.append_axes('right', size='5%', pad=0.02)
#     cbar = plt.colorbar(s_m, cax=ax_cbar,
#                         # ticks=[Wht.min(), Wht.max()],
#                         orientation='vertical')
#     cbar.set_ticks([Wht.min(), Wht.max()])
#     cbar.set_ticklabels([Alabel, Blabel])

ax.set_xlim(0, nsw)
ax.set_ylim(-3.0, 3.0)

ax.set_xlabel('Time [fs]',   fontsize='small', labelpad=5)
ax.set_ylabel('Energy [eV]', fontsize='small', labelpad=8)
ax.tick_params(which='both', labelsize='x-small')

########################################
plt.tight_layout(pad=0.2)
plt.savefig('ksen_wht.png', dpi=360)

############################################################
'''
############################################################
# load FSSH result files
############################################################

########################################
bmin     = 403
bmax     = 408
namdTime = 1000
potim    = 1.0
Nt       = namdTime - 2
inpFiles = glob('./SHPROP.*')
########################################

if not os.path.isfile('weight_sh.dat'):
    if (len(inpFiles) > 0):
        iniTimes = [int(F.split('.')[-1]) for F in inpFiles]

        psi_a    = np.array([np.loadtxt(F)[:,2:] for F in inpFiles])
        chgocc   = Wht[:,bmin:bmax]

        # chgocc   = Wht[:,1:,bmin:bmax]
        # chgocc[:,0,:] /= (chgocc[:,0,:] + chgocc[:,1,:])
        # chgocc[:,1,:] /= (chgocc[:,0,:] + chgocc[:,1,:])

        dcdt     = np.diff(chgocc, axis=0) / potim
        dpdt     = np.diff(psi_a,  axis=1) / potim

        rho      = []
        for start, ci in zip(iniTimes, psi_a):
            end  = start + namdTime
            P    = chgocc[start:end,:] * ci[:,:]
            # P    = chgocc[start:end,0,:] * ci[:,:]

            rho.append(np.sum(P, axis=1))
        else:
            rho = np.sum(rho, axis=0) / len(iniTimes)

        np.savetxt('weight_sh.dat', rho)

        ########################################
        Ns = len(inpFiles)

        NA = np.zeros((Ns,Nt))
        AD = np.zeros((Ns,Nt))
        ET = np.zeros((Ns,Nt, bmax-bmin))
        ########################################

        for i, j in enumerate(iniTimes):
            for k in range(Nt):
                NA[i,k] = np.sum(dpdt[i,k,:]  * chgocc[j-1+k,:])
                AD[i,k] = np.sum(psi_a[i,k,:] *   dcdt[j-1+k,:])
                # ET[i,k] = np.sum(psi_a[i,k,:] * chgocc[j-1+k,:])
                ET[i,k,:] = psi_a[i,k,:] * chgocc[j-1+k,:]

        # Average ET over all simulations
        et_avg = np.average(ET, axis=0)

        # Extract NTO and MSN bands from ET
        et_nto = et_avg[:, -2:]  # Last two bands (407, 408)
        et_msn = et_avg[:, :-2]  # Remaining bands (403-406)

        # Save the separate files
        np.save('et_nto.npy', et_nto)
        np.save('et_msn.npy', et_msn)

        na = np.average(NA, axis=0)
        ad = np.average(AD, axis=0)
        # et = np.average(ET, axis=0)

        np.save('na.npy', na)
        np.save('ad.npy', ad)
        # np.save('et.npy', et)

else:
    rho= np.loadtxt('weight_sh.dat')
    et_nto = np.load('et_nto.npy')
    et_msn = np.load('et_msn.npy')
    na = np.load('na.npy')
    ad = np.load('ad.npy')

############################################################
if os.path.isfile('et_nto.npy') and os.path.isfile('et_msn.npy') and os.path.isfile('na.npy') \
                            and os.path.isfile('ad.npy'):
    
    # Load the previously saved data
    et_nto = np.load('et_nto.npy')
    et_msn = np.load('et_msn.npy')
    na = np.load('na.npy')
    ad = np.load('ad.npy')

    # Save the data to a .dat file
    output_data = np.column_stack((np.arange(len(et_nto)), et_nto, et_msn, np.cumsum(na), np.cumsum(ad)))
    header = "Time [fs]    ET_NTO    ET_MSN    Cumulative NA    Cumulative AD"
    np.savetxt('hole_localization.dat', output_data, header=header, fmt='%12.6f')

    plt.clf()

    fig = plt.figure()
    fig.set_size_inches(3.6, 3.0)
    ax = plt.subplot()

    # Plot only one line from each 2D array
    ax.plot(et_nto[:, 1],       ls='-', lw=1.0, color='purple',    label='VBM@NTO')
    ax.plot(et_msn[:, 2],       ls='-', lw=1.0, color='darkgreen', label='VBM@MSN')
    ax.plot(np.cumsum(na),      ls='-', lw=1.0, color='#e377c2',   label='NA')
    ax.plot(np.cumsum(ad),      ls='-', lw=1.0, color='navy',      label='AD')

    ax.legend(loc='upper right')

    ax.set_xlim(0, namdTime)
    ax.set_ylim(-0.4, 0.8)

    ax.set_xlabel('Time [fs]', fontsize='small', labelpad=5)
    ax.set_ylabel('Hole Localization', fontsize='small', labelpad=5)

    plt.tight_layout(pad=0.2)
    plt.savefig('kspat_hole.png', dpi=300)
