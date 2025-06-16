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
pop_cbm1 = []  # Population of CBM+1
pop_vbm = []   # Population of VBM

for file in shprop_files:
    data = np.loadtxt(file)
    pop_cbm.append(data[:n_steps, 3])  # Column 4 = CBM
    pop_cbm1.append(data[:n_steps, 4]) # Column 5 = CBM+1
    pop_vbm.append(data[:n_steps, 2])   # Column 3 = VBM

# Average across trajectories
avg_pop_cbm = np.mean(pop_cbm, axis=0)
avg_pop_cbm1 = np.mean(pop_cbm1, axis=0)
avg_pop_vbm = np.mean(pop_vbm, axis=0)
time = np.arange(n_steps) * dt  # Time axis in fs

# Define exponential decay function for VBM
def exponential_decay(t, tau, amplitude):
    return amplitude * np.exp(-t / tau)

# Define exponential growth function for CBM and CBM+1
def exponential_growth(t, tau, amplitude):
    return amplitude * (1 - np.exp(-t / tau))

# Fit VBM population decay
popt_vbm, pcov_vbm = curve_fit(exponential_decay, time, avg_pop_vbm, p0=[100, 1.0])
tau_vbm_fs = popt_vbm[0]
tau_vbm_ps = tau_vbm_fs * 1e-3
print(f"Relaxation time (VBM Population Decay): {tau_vbm_ps:.3f} ps")

# Fit CBM population growth
popt_cbm, pcov_cbm = curve_fit(exponential_growth, time, avg_pop_cbm, p0=[100, 1.0])
tau_cbm_fs = popt_cbm[0]
tau_cbm_ps = tau_cbm_fs * 1e-3

# Fit CBM+1 population growth
popt_cbm1, pcov_cbm1 = curve_fit(exponential_growth, time, avg_pop_cbm1, p0=[100, 1.0])
tau_cbm1_fs = popt_cbm1[0]
tau_cbm1_ps = tau_cbm1_fs * 1e-3

# Save VBM, CBM, and CBM+1 population data
np.savetxt('vbm_decay.dat', np.column_stack((time, avg_pop_vbm, exponential_decay(time, *popt_vbm))), header='Time(fs) Population Fitted_Population')

# Normalize the population
norm_pop_cbm = avg_pop_cbm / np.max(avg_pop_cbm)
norm_pop_cbm1 = avg_pop_cbm1 / np.max(avg_pop_cbm1)
norm_pop_vbm = avg_pop_vbm / avg_pop_vbm[0]

# Define normalized decay and growth functions
def exp_decay_fixed(t, tau):
    return np.exp(-t / tau)

def exp_growth_fixed(t, tau):
    return (1 - np.exp(-t / tau))

# Fit normalized VBM decay
tau_norm_vbm, _ = curve_fit(exp_decay_fixed, time, norm_pop_vbm, p0=[100])
tau_norm_vbm_ps = tau_norm_vbm[0] * 1e-3
print(f"Normalized Relaxation time (VBM Decay): {tau_norm_vbm_ps:.3f} ps")

# Fit normalized CBM growth
tau_norm_cbm, _ = curve_fit(exp_growth_fixed, time, norm_pop_cbm, p0=[100])
tau_norm_cbm_ps = tau_norm_cbm[0] * 1e-3

# Fit normalized CBM+1 growth
tau_norm_cbm1, _ = curve_fit(exp_growth_fixed, time, norm_pop_cbm1, p0=[100])
tau_norm_cbm1_ps = tau_norm_cbm1[0] * 1e-3

# Save normalized VBM, CBM, and CBM+1 population data
np.savetxt('vbm_normalized_decay.dat', np.column_stack((time, norm_pop_vbm, exp_decay_fixed(time, tau_norm_vbm[0]))), header='Time(fs) Normalized_Population Fitted_Population')
np.savetxt('cbm_normalized_growth.dat', np.column_stack((time, norm_pop_cbm, exp_growth_fixed(time, tau_norm_cbm[0]))), header='Time(fs) Normalized_Population Fitted_Population')
np.savetxt('cbm1_normalized_growth.dat', np.column_stack((time, norm_pop_cbm1, exp_growth_fixed(time, tau_norm_cbm1[0]))), header='Time(fs) Normalized_Population Fitted_Population')

# Combined plot for VBM decay and CBM, CBM+1 growth
plt.figure(figsize=(8, 5))
plt.plot(time, avg_pop_vbm, 'b-', label='VBM Population')
plt.plot(time, exponential_decay(time, *popt_vbm), 'r--', label=f'VBM Fit: τ = {tau_vbm_ps:.3f} ps')
plt.plot(time, avg_pop_cbm, 'g-', label='CBM Population')
plt.plot(time, exponential_growth(time, *popt_cbm), 'm--')
plt.plot(time, avg_pop_cbm1, 'c-', label='CBM+1 Population')
plt.plot(time, exponential_growth(time, *popt_cbm1), 'y--')
plt.xlabel('Time (fs)')
plt.ylabel('Population')
plt.legend()
plt.tight_layout()
plt.savefig('reverse_recomb.png', dpi=300)

# Combined plot for normalized CBM decay and VBM growth
plt.figure(figsize=(8, 5))
plt.plot(time, norm_pop_vbm, 'b-', label='Normalized VBM Population')
plt.plot(time, exp_decay_fixed(time, tau_norm_vbm[0]), 'r--', label=f'VBM Fit: τ = {tau_norm_vbm_ps:.3f} ps')
plt.plot(time, norm_pop_cbm, 'g-', label='Normalized CBM Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_cbm[0]), 'm--')
plt.plot(time, norm_pop_cbm1, 'c-', label='Normalized CBM+1 Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_cbm1[0]), 'y--')
plt.xlabel('Time (fs)')
plt.ylabel('Normalized Population')
plt.legend()
plt.tight_layout()
plt.savefig('reverse_recomb_norm.png', dpi=300)