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
pop_cbm2 = []  # Population of CBM+2
pop_cbm1 = []  # Population of CBM+1
pop_cbm = []   # Population of CBM

for file in shprop_files:
    data = np.loadtxt(file)
    pop_cbm2.append(data[:n_steps, 4])  # Column 5 = CBM+2
    pop_cbm1.append(data[:n_steps, 3])  # Column 4 = CBM+1
    pop_cbm.append(data[:n_steps, 2])   # Column 3 = CBM

# Average across trajectories
avg_pop_cbm2 = np.mean(pop_cbm2, axis=0)
avg_pop_cbm1 = np.mean(pop_cbm1, axis=0)
avg_pop_cbm = np.mean(pop_cbm, axis=0)
time = np.arange(n_steps) * dt  # Time axis in fs

# Define exponential decay function for CBM+2
def exponential_decay(t, tau, amplitude):
    return amplitude * np.exp(-t / tau)

# Define exponential growth function for CBM and CBM+1
def exponential_growth(t, tau, amplitude):
    return amplitude * (1 - np.exp(-t / tau))

# Fit CBM+2 population decay
popt_cbm2, pcov_cbm2 = curve_fit(exponential_decay, time, avg_pop_cbm2, p0=[100, 1.0])
tau_cbm2_fs = popt_cbm2[0]
tau_cbm2_ps = tau_cbm2_fs * 1e-3
# print(f"Relaxation time (CBM+2 Population Decay): {tau_cbm2_ps:.3f} ps")

# Fit CBM population growth
popt_cbm, pcov_cbm = curve_fit(exponential_growth, time, avg_pop_cbm, p0=[100, 1.0])
tau_cbm_fs = popt_cbm[0]
tau_cbm_ps = tau_cbm_fs * 1e-3
# print(f"Relaxation time (CBM Population Growth): {tau_cbm_ps:.3f} ps")

# Fit CBM+1 population growth
popt_cbm1, pcov_cbm1 = curve_fit(exponential_growth, time, avg_pop_cbm1, p0=[100, 1.0])
tau_cbm1_fs = popt_cbm1[0]
tau_cbm1_ps = tau_cbm1_fs * 1e-3
# print(f"Relaxation time (CBM+1 Population Growth): {tau_cbm1_ps:.3f} ps")

# Save CBM+2, CBM, and CBM+1 population data
np.savetxt('cbm2_decay.dat', np.column_stack((time, avg_pop_cbm2, exponential_decay(time, *popt_cbm2))), header='Time(fs) Population Fitted_Population')
np.savetxt('cbm_growth.dat', np.column_stack((time, avg_pop_cbm, exponential_growth(time, *popt_cbm))), header='Time(fs) Population Fitted_Population')
np.savetxt('cbm1_growth.dat', np.column_stack((time, avg_pop_cbm1, exponential_growth(time, *popt_cbm1))), header='Time(fs) Population Fitted_Population')

# Normalize the population
norm_pop_cbm2 = avg_pop_cbm2 / np.max(avg_pop_cbm2)
norm_pop_cbm = avg_pop_cbm / np.max(avg_pop_cbm)
norm_pop_cbm1 = avg_pop_cbm1 / np.max(avg_pop_cbm1)

# Define normalized decay and growth functions
def exp_decay_fixed(t, tau):
    return np.exp(-t / tau)

def exp_growth_fixed(t, tau):
    return (1 - np.exp(-t / tau))

# Fit normalized CBM+2 decay
tau_norm_cbm2, _ = curve_fit(exp_decay_fixed, time, norm_pop_cbm2, p0=[100])
tau_norm_cbm2_ps = tau_norm_cbm2[0] * 1e-3
print(f"Normalized Relaxation time (CBM+2 Decay): {tau_norm_cbm2_ps:.3f} ps")

# Fit normalized CBM growth
tau_norm_cbm, _ = curve_fit(exp_growth_fixed, time, norm_pop_cbm, p0=[100])
tau_norm_cbm_ps = tau_norm_cbm[0] * 1e-3
print(f"Normalized Relaxation time (CBM Growth): {tau_norm_cbm_ps:.3f} ps")

# Fit normalized CBM+1 growth
tau_norm_cbm1, _ = curve_fit(exp_growth_fixed, time, norm_pop_cbm1, p0=[100])
tau_norm_cbm1_ps = tau_norm_cbm1[0] * 1e-3
print(f"Normalized Relaxation time (CBM+1 Growth): {tau_norm_cbm1_ps:.3f} ps")

# Save normalized CBM+2, CBM, and CBM+1 population data
np.savetxt('cbm2_normalized_decay.dat', np.column_stack((time, norm_pop_cbm2, exp_decay_fixed(time, tau_norm_cbm2[0]))), header='Time(fs) Normalized_Population Fitted_Population')
np.savetxt('cbm_normalized_growth.dat', np.column_stack((time, norm_pop_cbm, exp_growth_fixed(time, tau_norm_cbm[0]))), header='Time(fs) Normalized_Population Fitted_Population')
np.savetxt('cbm1_normalized_growth.dat', np.column_stack((time, norm_pop_cbm1, exp_growth_fixed(time, tau_norm_cbm1[0]))), header='Time(fs) Normalized_Population Fitted_Population')

# Combined plot for CBM+2 decay and CBM, CBM+1 growth
plt.figure(figsize=(8, 5))
plt.plot(time, avg_pop_cbm2, 'b-', label='CBM+2 Population')
plt.plot(time, exponential_decay(time, *popt_cbm2), 'r--', label=f'CBM+2 Fit: τ = {tau_cbm2_ps:.3f} ps')
plt.plot(time, avg_pop_cbm, 'g-', label='CBM Population')
plt.plot(time, exponential_growth(time, *popt_cbm), 'm--', label=f'CBM Fit: τ = {tau_cbm_ps:.3f} ps')
plt.plot(time, avg_pop_cbm1, 'c-', label='CBM+1 Population')
plt.plot(time, exponential_growth(time, *popt_cbm1), 'y--', label=f'CBM+1 Fit: τ = {tau_cbm1_ps:.3f} ps')
plt.xlabel('Time (fs)')
plt.ylabel('Population')
plt.legend()
plt.tight_layout()
plt.savefig('electron_population.png', dpi=300)

# Combined plot for normalized CBM+2 decay and CBM growth
plt.figure(figsize=(8, 5))
plt.plot(time, norm_pop_cbm2, 'b-', label='Normalized CBM+2 Population')
plt.plot(time, exp_decay_fixed(time, tau_norm_cbm2[0]), 'r--', label=f'CBM+2 Fit: τ = {tau_norm_cbm2_ps:.3f} ps')
plt.plot(time, norm_pop_cbm, 'g-', label='Normalized CBM Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_cbm[0]), 'm--', label=f'CBM Fit: τ = {tau_norm_cbm_ps:.3f} ps')
plt.plot(time, norm_pop_cbm1, 'c-', label='Normalized CBM+1 Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_cbm1[0]), 'y--', label=f'CBM+1 Fit: τ = {tau_norm_cbm1_ps:.3f} ps')
plt.xlabel('Time (fs)')
plt.ylabel('Normalized Population')
plt.legend()
plt.tight_layout()
plt.savefig('electron_population_norm.png', dpi=300)