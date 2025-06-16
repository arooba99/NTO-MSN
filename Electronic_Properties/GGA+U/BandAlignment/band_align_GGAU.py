import matplotlib.pyplot as plt

# === Heterostructure ===
VBM = 0.3628          # Valence Band Maximum (in eV)
CBM = 1.5121           # Conduction Band Minimum (in eV)
E_vacuum = 5.594       # Vacuum level (in eV)

# VBM = 0.2029          # Valence Band Maximum (in eV)
# CBM = 1.5073           # Conduction Band Minimum (in eV)
# E_vacuum = 5.521       # Vacuum level (in eV)

# === NTO ===
# VBM = -3.9425          # Valence Band Maximum (in eV)
# CBM = -2.0787           # Conduction Band Minimum (in eV)
# E_vacuum = 1.808       # Vacuum level (in eV)

# === MSN ===
# VBM = 1.1745          # Valence Band Maximum (in eV)
# CBM = 3.1303           # Conduction Band Minimum (in eV)
# E_vacuum = 6.865       # Vacuum level (in eV)

# === Shifted to Vacuum Level ===
VBM_vac = VBM - E_vacuum
CBM_vac = CBM - E_vacuum

# Water redox levels (in eV vs vacuum)
water_reduction = -4.44  # H⁺/H₂
water_oxidation = -5.67  # O₂/H₂O

# === Plotting ===
plt.figure(figsize=(5, 6))
plt.bar(1, VBM_vac - CBM_vac, bottom=CBM_vac, width=0.4, color='skyblue', edgecolor='k', label='Band Gap')

# Redox and vacuum level lines
plt.axhline(y=water_reduction, color='red', linestyle='--', label='Water Reduction (-4.44 eV)')
plt.axhline(y=water_oxidation, color='green', linestyle='--', label='Water Oxidation (-5.67 eV)')
plt.axhline(y=0, color='gray', linestyle=':')

# Formatting
plt.xlim(0.5, 1.5)
plt.ylim(-3, -7)
plt.ylabel('Energy (eV vs Vacuum)', fontsize=12)
plt.xticks([])
plt.xlabel('')
plt.gca().invert_yaxis()  # Vacuum at top
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('band_alignment_GGA.png', dpi=300)