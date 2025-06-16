import warnings
from pymatgen.core import Structure
from m3gnet.models import MolecularDynamics

for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")

structure = Structure.from_file("POSCAR.vasp")

md = MolecularDynamics(
    atoms=structure,
    temperature=300,  # 300 K
    ensemble="nvt",  # NVT ensemble
    timestep=0.2,  # 0.2 fs,
    trajectory="nvt_run.traj",  # save trajectory to .traj
    logfile="nvt_run.log",  # log file for MD
    loginterval=5,  # record every 1 fs
)

md.run(steps=35000)  # 0.2 fs * 35000 = 7000 fs = 7ps

from ase.io import read

# Read the last frame from the trajectory file
final_structure = read("nvt_run.traj", index=-1)

# Save final structure in VASP POSCAR format
final_structure.write("CONTCAR")


import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("nvt_run.log", sep='\s+', comment="#")

# Create figure and 3 subplots arranged vertically (3 rows, 1 column)
fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# Plot Temperature vs Time
axs[0].plot(df["Time[ps]"], df["T[K]"], color='tab:red')
axs[0].set_ylabel("Temperature (K)", fontsize=16)
axs[0].grid(True)
axs[0].tick_params(axis='both', labelsize=14)

# Plot Potential Energy vs Time
axs[1].plot(df["Time[ps]"], df["Epot[eV]"], color='tab:blue')
axs[1].set_ylabel("Potential Energy (eV)", fontsize=16)
axs[1].grid(True)
axs[1].tick_params(axis='both', labelsize=14)

# Plot Total Energy vs Time
axs[2].plot(df["Time[ps]"], df["Etot[eV]"], color='tab:green')
axs[2].set_xlabel("Time (ps)", fontsize=16)
axs[2].set_ylabel("Total Energy (eV)", fontsize=16)
axs[2].grid(True)
axs[2].tick_params(axis='both', labelsize=14)

plt.tight_layout()
plt.savefig("aimd.png", dpi=300)