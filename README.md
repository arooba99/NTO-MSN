# Ni₃TeO₆/MoSi₂N₄ Heterostructure for Photocatalytic Water Splitting

This repository contains data, structures, scripts, and supporting calculations for "Ultrafast Charge Transfer and Delayed Recombination in type-II 
Ni₃TeO₆/MoSi₂N₄ Heterostructure: A Time-Domain Ab Initio Study".

## Repository Structure

### 1. `Structure/`
- POSCAR/CIF files for:
  - Ni₃TeO₆ (NTO)
  - MoSi₂N₄ (MSN)
  - NTO/MSN heterostructure

---

### 2. `Stability/`
- **Dynamic Stability:**
  - Phonon dispersion plots (PHONOPY)
- **Mechanical Stability:**
  - Stiffness tensor (in N/m)
- **Thermodynamic Stability:**
  - AIMD input/output files and representative trajectory

---

### 3. `Electronic_Properties/`
- Band structures and DOS plots (GGA+U and HSE)
- Scripts for band alignment and electric field calculation
- Work function
- Charge density difference (CDD) plot

---

### 4. `NAMD/`
- Scripts and selected outputs for:
  - Time-dependent Kohn-Sham energies
  - Carrier populations
  - Carrier Localization

---

### 5. `Gibbs_HER/`
- POSCARs and CONTCARs for pristine and defect systems for Gibbs free energy calculations for HER

---

### 6. `STH_Efficiency/`
- Script for solar-to-hydrogen (STH) efficiency estimation

---

## Notes
- All calculations were performed using VASP and analyzed using Python-based post-processing tools.
- Gibbs free energy calculations were performed using VASP with vibrational corrections.
- GGA+U functional was adopted for all calculations (ground-state properties, excited-state dynamics and vibrational analysis)
- NAMD simulations were performed using the Hefei-NAMD package on 2x2x1 supercells with ~2000 SCF steps.

---
