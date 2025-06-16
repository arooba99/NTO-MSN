#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Load SHPROP files
shprop_files = glob('SHPROP.*')
n_steps = 1000  # Adjust based on your NAMDTIME
dt = 1.0  # Timestep (fs) from POTIM

# Initialize arrays
pop_cbm = []  # Population of CBM
pop_vbm = []   # Population of VBM
pop_vbm1 = []  # Population of VBM-1

for file in shprop_files:
    data = np.loadtxt(file)
    pop_cbm.append(data[:n_steps, 5])  # Column 6 = CBM
    pop_vbm.append(data[:n_steps, 4])   # Column 5 = VBM
    pop_vbm1.append(data[:n_steps, 3])  # Column 4 = VBM-1

# Average across trajectories
avg_pop_cbm = np.mean(pop_cbm, axis=0)
avg_pop_vbm = np.mean(pop_vbm, axis=0)
avg_pop_vbm1 = np.mean(pop_vbm1, axis=0)
time = np.arange(n_steps) * dt  # Time axis in fs

# Define exponential decay function for CBM
def exponential_decay(t, tau, amplitude):
    return amplitude * np.exp(-t / tau)

# Define exponential growth function for VBM and VBM-1
def exponential_growth(t, tau, amplitude):
    return amplitude * (1 - np.exp(-t / tau))

# Fit CBM population decay
popt_cbm, pcov_cbm = curve_fit(exponential_decay, time, avg_pop_cbm, p0=[100, 1.0])
tau_cbm_fs = popt_cbm[0]
tau_cbm_ps = tau_cbm_fs * 1e-3
print(f"Relaxation time (CBM Population Decay): {tau_cbm_ps:.3f} ps")

# Fit VBM population growth
popt_vbm, pcov_vbm = curve_fit(exponential_growth, time, avg_pop_vbm, p0=[100, 1.0])
tau_vbm_fs = popt_vbm[0]
tau_vbm_ps = tau_vbm_fs * 1e-3

# Fit VBM-1 population growth
popt_vbm1, pcov_vbm1 = curve_fit(exponential_growth, time, avg_pop_vbm1, p0=[100, 1.0])
tau_vbm1_fs = popt_vbm1[0]
tau_vbm1_ps = tau_vbm1_fs * 1e-3

# Save CBM, VBM, and VBM-1 population data
np.savetxt('cbm_decay.dat', np.column_stack((time, avg_pop_cbm, exponential_decay(time, *popt_cbm))), header='Time(fs) Population Fitted_Population')

# Normalize the population
norm_pop_cbm = avg_pop_cbm / avg_pop_cbm[0]
norm_pop_vbm = avg_pop_vbm / np.max(avg_pop_vbm)
norm_pop_vbm1 = avg_pop_vbm1 / np.max(avg_pop_vbm1)

# Define normalized decay and growth functions
def exp_decay_fixed(t, tau):
    return np.exp(-t / tau)

def exp_growth_fixed(t, tau):
    return (1 - np.exp(-t / tau))

# Fit normalized CBM decay
tau_norm_cbm, _ = curve_fit(exp_decay_fixed, time, norm_pop_cbm, p0=[100])
tau_norm_cbm_ps = tau_norm_cbm[0] * 1e-3
print(f"Normalized Relaxation time (CBM Decay): {tau_norm_cbm_ps:.3f} ps")

# Fit normalized VBM growth
tau_norm_vbm, _ = curve_fit(exp_growth_fixed, time, norm_pop_vbm, p0=[100])
tau_norm_vbm_ps = tau_norm_vbm[0] * 1e-3

# Fit normalized VBM-1 growth
tau_norm_vbm1, _ = curve_fit(exp_growth_fixed, time, norm_pop_vbm1, p0=[100])
tau_norm_vbm1_ps = tau_norm_vbm1[0] * 1e-3

# Save normalized CBM, VBM, and VBM-1 population data
np.savetxt('cbm_normalized_decay.dat', np.column_stack((time, norm_pop_cbm, exp_decay_fixed(time, tau_norm_cbm[0]))), header='Time(fs) Normalized_Population Fitted_Population')


# Combined plot for CBM decay and VBM, VBM-1 growth
plt.figure(figsize=(8, 5))
plt.plot(time, avg_pop_cbm, 'b-', label='CBM Population')
plt.plot(time, exponential_decay(time, *popt_cbm), 'r--', label=f'CBM Fit: τ = {tau_cbm_ps:.3f} ps')
plt.plot(time, avg_pop_vbm, 'g-', label='VBM Population')
plt.plot(time, exponential_growth(time, *popt_vbm), 'm--')
plt.plot(time, avg_pop_vbm1, 'c-', label='VBM-1 Population')
plt.plot(time, exponential_growth(time, *popt_vbm1), 'y--')
plt.xlabel('Time (fs)')
plt.ylabel('Population')
plt.legend()
plt.tight_layout()
plt.savefig('direct_recomb.png', dpi=300)

# Combined plot for normalized CBM decay and VBM growth
plt.figure(figsize=(8, 5))
plt.plot(time, norm_pop_cbm, 'b-', label='Normalized CBM Population')
plt.plot(time, exp_decay_fixed(time, tau_norm_cbm[0]), 'r--', label=f'CBM Fit: τ = {tau_norm_cbm_ps:.3f} ps')
plt.plot(time, norm_pop_vbm, 'g-', label='Normalized VBM Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_vbm[0]), 'm--')
plt.plot(time, norm_pop_vbm1, 'c-', label='Normalized VBM-1 Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_vbm1[0]), 'y--')
plt.xlabel('Time (fs)')
plt.ylabel('Normalized Population')
plt.legend()
plt.tight_layout()
plt.savefig('direct_recomb_norm.png', dpi=300)