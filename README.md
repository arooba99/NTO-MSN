# Ni₃TeO₆/MoSi₂N₄ Heterostructure for Photocatalytic Water Splitting

This repository contains data, structures, scripts, and supporting calculations for "Ultrafast Charge Transfer and Delayed Recombination in type-II 
Ni₃TeO₆/MoSi₂N₄ Heterostructure: A Time-Domain Ab Initio Study".

## Repository Structure

### 1. `Structure/`
- Input and output files for
  - Ni₃TeO₆ (NTO)
  - MoSi₂N₄ (MSN)
  - NTO/MSN heterostructure

---

### 2. `Stability/`
- **Mechanical Stability**
  - Stiffness tensor (in N/m)
- **Thermodynamic Stability (AIMD)**
  - AIMD input/output files and representative trajectory
- **AIMD_M3GNet**
  - Stability test using M3GNet pretrained model

---

### 3. `Electronic_Properties/`
- Files are provided for comparative analysis between GGA+U, HSE06 and HSE06+U.
- Projected Band structures and DOS data files
- Scripts for band alignment and electric field calculation
- Work function
- Charge density difference (CDD)
- Bader charge analysis

---

### 4. `NAMD/`
- **01-Supercell_optimization**
  - INCAR file indicating input parameters used
- **02-AIMD_NVT_Equilibration_run**
  - Example INCAR and KPOINTS files along with post-processing bash script (dynamics.sh)
- **03-AIMD_NVE_Production_run**
  - Example INCAR and KPOINTS files
- **04-NAMD_run**
  - xdat2pos.py script to obtain input POSCAR files from last NSCF steps (from NVE run) for SCF calculations
  - Example INCAR and KPOINTS files for scf-based runs
- **05-Post_calculation_scripts**
  - Example scripts and outputs for
    - Time-dependent Kohn-Sham energies
    - Carrier lifetimes
    - Carrier Localization
    - Nonadiabatic Couplings (NACs)
    - Dephasing function, Autocorrelation function (ACF), and Spectral density 

---

### 5. `Gibbs_HER/`
- Input and output files for selected adsorption sites in pristine and vacancy-defective systems for Gibbs free energy calculations for HER

---

### 6. `STH_Efficiency/`
- Script for solar-to-hydrogen (STH) efficiency estimation using GGA+U, HSE06 and HSE06+U functionals.

---

## Notes
- All ground-state calculations were performed using VASP and analyzed using vaspkit and python-based post-processing scripts.
- GGA+U functional with vdW correction (D3) was adopted for all calculations (ground-state properties, excited-state dynamics and vibrational analysis)
- Since HSE06 and HSE06+U are computationally intensive functionals for this magnetic system, GGA+U was adopted for spin-polarized calculations here.
- Excited-state NAMD simulations were performed using the Hefei-NAMD package on 2x2x1 supercells with 2000 SCF steps.

---
