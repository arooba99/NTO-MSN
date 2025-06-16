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
pop_vbm2 = []  # Population of VBM-2
pop_vbm1 = []  # Population of VBM-1
pop_vbm = []   # Population of VBM

for file in shprop_files:
    data = np.loadtxt(file)
    pop_vbm2.append(data[:n_steps, 4])  # Column 5 = VBM-2
    pop_vbm1.append(data[:n_steps, 5])  # Column 6 = VBM-1
    pop_vbm.append(data[:n_steps, 6])   # Column 7 = VBM

# Average across trajectories
avg_pop_vbm2 = np.mean(pop_vbm2, axis=0)
avg_pop_vbm1 = np.mean(pop_vbm1, axis=0)
avg_pop_vbm = np.mean(pop_vbm, axis=0)
time = np.arange(n_steps) * dt  # Time axis in fs

# Define exponential decay and growth functions
def exponential_decay(t, tau, amplitude):
    return amplitude * np.exp(-t / tau)

def exponential_growth(t, tau, amplitude):
    return amplitude * (1 - np.exp(-t / tau))

# Fit VBM-2 population decay
popt_vbm2, _ = curve_fit(exponential_decay, time, avg_pop_vbm2, p0=[100, 1.0])

# Fit VBM and VBM-1 population growth
popt_vbm, _ = curve_fit(exponential_growth, time, avg_pop_vbm, p0=[100, 1.0])
popt_vbm1, _ = curve_fit(exponential_growth, time, avg_pop_vbm1, p0=[100, 1.0])

# Normalize the population
norm_pop_vbm2 = avg_pop_vbm2 / avg_pop_vbm2[0]
norm_pop_vbm = avg_pop_vbm / np.max(avg_pop_vbm)
norm_pop_vbm1 = avg_pop_vbm1 / np.max(avg_pop_vbm1)

# Define normalized decay and growth functions
def exp_decay_fixed(t, tau):
    return np.exp(-t / tau)

def exp_growth_fixed(t, tau):
    return (1 - np.exp(-t / tau))

# Fit normalized VBM-2 decay
tau_norm_vbm2, _ = curve_fit(exp_decay_fixed, time, norm_pop_vbm2, p0=[0.1])
tau_norm_vbm2_ps = tau_norm_vbm2[0] * 1e-3
print(f"Normalized Relaxation time (VBM-2 Decay): {tau_norm_vbm2_ps:.3f} ps")

# Fit normalized VBM and VBM-1 growth
tau_norm_vbm, _ = curve_fit(exp_growth_fixed, time, norm_pop_vbm, p0=[0.1])
tau_norm_vbm_ps = tau_norm_vbm[0] * 1e-3
print(f"Normalized Relaxation time (VBM Growth): {tau_norm_vbm_ps:.3f} ps")

tau_norm_vbm1, _ = curve_fit(exp_growth_fixed, time, norm_pop_vbm1, p0=[0.1])
tau_norm_vbm1_ps = tau_norm_vbm1[0] * 1e-3
print(f"Normalized Relaxation time (VBM-1 Growth): {tau_norm_vbm1_ps:.3f} ps")

# Save results
np.savetxt('vbm2_decay.dat', np.column_stack((time, avg_pop_vbm2, exponential_decay(time, *popt_vbm2))), header='Time(fs) Population Fitted_Population')
np.savetxt('vbm_growth.dat', np.column_stack((time, avg_pop_vbm, exponential_growth(time, *popt_vbm))), header='Time(fs) Population Fitted_Population')
np.savetxt('vbm1_growth.dat', np.column_stack((time, avg_pop_vbm1, exponential_growth(time, *popt_vbm1))), header='Time(fs) Population Fitted_Population')

# Save normalized data
np.savetxt('vbm2_normalized_decay.dat', np.column_stack((time, norm_pop_vbm2, exp_decay_fixed(time, tau_norm_vbm2[0]))), header='Time(fs) Normalized_Population Fitted_Population')
np.savetxt('vbm_normalized_growth.dat', np.column_stack((time, norm_pop_vbm, exp_growth_fixed(time, tau_norm_vbm[0]))), header='Time(fs) Normalized_Population Fitted_Population')
np.savetxt('vbm1_normalized_growth.dat', np.column_stack((time, norm_pop_vbm1, exp_growth_fixed(time, tau_norm_vbm1[0]))), header='Time(fs) Normalized_Population Fitted_Population')

# Plot for VBM-2, VBM, and VBM-1
plt.figure(figsize=(8, 5))
plt.plot(time, avg_pop_vbm2, 'b-', label='VBM-2 Population')
plt.plot(time, exponential_decay(time, *popt_vbm2), 'r--', label='VBM-2 Fit')
plt.plot(time, avg_pop_vbm, 'g-', label='VBM Population')
plt.plot(time, exponential_growth(time, *popt_vbm), 'm--', label='VBM Fit')
plt.plot(time, avg_pop_vbm1, 'c-', label='VBM-1 Population')
plt.plot(time, exponential_growth(time, *popt_vbm1), 'y--', label='VBM-1 Fit')
plt.xlabel('Time (fs)')
plt.ylabel('Population')
plt.legend()
plt.tight_layout()
plt.savefig('hole_population.png', dpi=300)

# Combined plot for normalized VBM-2 decay and VBM growth
plt.figure(figsize=(8, 5))
plt.plot(time, norm_pop_vbm2, 'b-', label='Normalized VBM-2 Population')
plt.plot(time, exp_decay_fixed(time, tau_norm_vbm2[0]), 'r--', label=f'VBM-2 Fit: τ = {tau_norm_vbm2_ps:.3f} ps')
plt.plot(time, norm_pop_vbm, 'g-', label='Normalized VBM Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_vbm[0]), 'm--', label=f'VBM Fit: τ = {tau_norm_vbm_ps:.3f} ps')
plt.plot(time, norm_pop_vbm1, 'c-', label='Normalized VBM-1 Population')
plt.plot(time, exp_growth_fixed(time, tau_norm_vbm1[0]), 'y--', label=f'VBM-1 Fit: τ = {tau_norm_vbm1_ps:.3f} ps')
plt.xlabel('Time (fs)')
plt.ylabel('Normalized Population')
plt.legend()
plt.tight_layout()
plt.savefig('hole_population_norm.png', dpi=300)