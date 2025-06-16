#!/usr/bin/env python3.6

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
#from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# def gaussian(t, c, tau):
#     return np.exp(-0.5 * (t / tau)**2)

def gaussian(x,c):
    result = np.exp(-x**2/(2*c**2))
    return result

def dephase(Et, dt=1.0):
    r'''
    Calculate the autocorrelation function (ACF), dephasing function, and FT of
    ACF.

    The dephasing function was calculated according to the following formula:

    G(t) = (1 / hbar**2) \int_0^{t_1} dt_1 \int_0^{t_2} dt_2 <E(t_2)E(0)>
    D(t) = exp(-G(t))

    Fourier Transform (FT) of the normalized ACF gives the phonon influence
    spectrum, also known as the phonon spectral density.

    I(\omega) \propto | FT(Ct / Ct[0]) |**2
    '''

    from scipy.fftpack import fft

    hbar = 0.6582119513926019  # eV fs

    Et = np.asarray(Et)
    Et -= np.average(Et)

    # Autocorrelation Function (ACF) of Et
    Ct = np.correlate(Et, Et, 'full')[Et.size:] / Et.size
    
    # Cumulative integration of the ACF
    Gt = cumulative_trapezoid(cumulative_trapezoid(Ct, dx=dt, initial=0), dx=dt, initial=0)
    Gt /= hbar**2
    
    # Dephasing function
    Dt = np.exp(-Gt)
    
    # FT of normalized ACF
    Iw = np.abs(fft(Ct / Ct[0]))**2

    return Ct, Dt, Gt, Iw


energy = np.loadtxt('EIGTXT')
nbasis = energy.shape[1]
matrix = np.zeros((nbasis, nbasis), dtype=float)

dt = 1.0 # fs

############################################################
# Initialize accumulators for averaging
T_all = None
Et_all = None
Ct_all = None
Dt_all = None
Gt_all = None
Dt_fit_all = None
Iw_all = None

count = 0  # To keep track of how many pairs are accumulated

############################################################
# Loop through each pair of basis states
# Loop through each pair of basis states
for ii in range(nbasis):
    for jj in range(ii):
        Et = energy[:, ii] - energy[:, jj]
        T = np.arange(Et.size-1) * dt

        # Perform dephasing calculations
        Ct, Dt, Gt, Iw = dephase(Et)

        # Fit dephasing function using Gaussian
        popt, pcov = curve_fit(gaussian, T, Dt, maxfev=10000)
        Dt_fit = gaussian(T, *popt)
        matrix[ii, jj] = popt[0]
        matrix[jj, ii] = matrix[ii, jj]

        # Accumulate data for averaging
        if T_all is None:
            # Initialize accumulators
            T_all = T
            Et_all = Et[:-1]
            Ct_all = Ct
            Dt_all = Dt
            Gt_all = Gt  # Initialize Gt_all accumulator
            Dt_fit_all = Dt_fit
            Iw_all = Iw
        else:
            # Accumulate by summing
            Et_all += Et[:-1]
            Ct_all += Ct
            Dt_all += Dt
            Gt_all += Gt  # Accumulate Gt
            Dt_fit_all += Dt_fit
            Iw_all += Iw
        
        count += 1  # Increment pair count

############################################################
# Calculate means by dividing the accumulators by the number of pairs
Et_avg = Et_all / count
Ct_avg = Ct_all / count
Dt_avg = Dt_all / count
Gt_avg = Gt_all / count  # Calculate Gt average
Dt_fit_avg = Dt_fit_all / count
Iw_avg = Iw_all / count

############################################################
# Save the averaged data to a single .dat file
with open('dephasing_averaged_data.dat', 'w') as f:
    f.write("# Time [fs]    Et_avg [eV]    Ct_avg    Dt_avg    Gt_avg    Dt_fit_avg    Iw_avg\n")
    for t, e, ct, dt, gt, dtf, iw in zip(T_all, Et_avg, Ct_avg, Dt_avg, Gt_avg, Dt_fit_avg, Iw_avg):
        f.write(f"{t:12.4f} {e:12.6f} {ct:12.6f} {dt:12.6f} {gt:12.6f} {dtf:12.6f} {iw:12.6f}\n")

############################################################
# Save DEPHTIME matrix to file
np.savetxt('DEPHTIME', matrix, fmt='%10.4f')

############################################################
fig, axes = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.10, top=0.95, wspace=0.10, hspace=0.20)
fig.set_size_inches((4.5, 9))

inp = [Et_avg, Ct_avg, Dt_avg, Iw_avg]
ylabs = [r'$\Delta$E (eV)', 'ACF', 'Dephasing Function', 'Spectral Density']
clrs = ['k', 'r', 'b', 'g']

for ii in range(4):
    ax = axes.flatten()[ii]
    dat = inp[ii]
    ax.minorticks_on()
    N = min(T_all.size, dat.size)
    ax.plot(T_all[:N], dat[:N], ls='-', lw=1.0, markevery=3, mew=0.0, ms=3, color=clrs[ii])
    ax.set_ylabel(ylabs[ii], labelpad=5, fontsize='medium')

if ii == 3:
    ax.set_xlabel('Time [fs]', labelpad=10, fontsize='medium')
    ax.set_xlim(0, 500)
    ax.tick_params(which='both', labelsize='medium')

plt.tight_layout()
plt.savefig('averaged_decoherence.png', dpi=300)